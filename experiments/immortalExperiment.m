rng shuffle

% Load network and distance matrices and glomeruli
load('fullNetworkImmortal.mat') % Load your network data
load('distanceImmortal.mat') % Load you distance data
load('glomeruliImmortal.mat') % Load your glomeruli data
load('mitralCellsImmortal.mat')

% Number of mitral cells and granule cells
[mitralNum,granuleNum] = size(network);

% Set up the glomeruli (just an array from 1 to the total number of
% glomeruli)
glomeruli = unique(glomArray);
glomNum = length(glomeruli);

% Intrinsic MC and GC parameters
mParam = zeros(9, mitralNum);
gParam = zeros(8, granuleNum);

% Establish the intrinsic MC parameters
for i = 1:mitralNum
    while mParam(1,i) >= mParam(2,i)
        mParam(1,i) = normrnd(-58,58/10); % mVr (resting potential)
        mParam(2,i) = normrnd(-49,49/10); % mVt (threshold potential)
    end
end
mParam(3,:) = normrnd(0.02,0.02/10,1,mitralNum); %a (Izhikevich parameter a)
mParam(4,:) = normrnd(12,12/10,1,mitralNum); %b_m (Izhikevich parameter b)
mParam(5,:) = normrnd(2.5,2.5/10,1,mitralNum); %k_m (Izhikevich parameter k)

for i = 1:mitralNum
    while mParam(6,i) >= mParam(2,i)
        mParam(6,i) = normrnd(-65,65/10); % c_m (Izhikevich parameter c)
    end
end
mParam(7,:) = normrnd(13,13/10,1,mitralNum); %d_m (Izhikevich parameter d)  
mParam(8,:) = normrnd(191,191/10,1,mitralNum); % cap_m (capacitance)
mParam(9,:) = normrnd(0.006,0.006/10,mitralNum,1);

% Set the length constant
L = normrnd(675,675/10,mitralNum,1);
Lmat = repmat(L,[1,granuleNum]); 

% MC synaptic parameters 
mGABA = normrnd(0.13,0.13/10,mitralNum,granuleNum).* exp(-distance./Lmat); % GABA adjusted for distance
tauG_m = normrnd(18,18/10,mitralNum, granuleNum); % tau_GABA


% Establish the intrinsic granule cell parameters
for i = 1:granuleNum
    while gParam(1,i) >= gParam(2,i) || gParam(2,i) - gParam(1,i) < 20
        gParam(1,i) = normrnd(-71,71/10); % gVr (resting potential)
        gParam(2,i) = normrnd(-39,39/10); % gVt (threshold potential)
    end
end
gParam(3,:) = normrnd(0.01,0.01/10,1,granuleNum); %a_g (Izhikevich parameter a)

for i = 1:granuleNum
    % Establish Izhikevich parameters b and k based on rheobase and input resistance
    rheobase = 0;
    inres = -1;
    while (rheobase < 10 || rheobase > 70 || inres < 0.25 || inres > 1.5) || b > 0
        b = normrnd(-2/15,2/10); 
        k = normrnd(1/15,1/10);
        rheobase = (b+k*(-gParam(1,i)+gParam(2,i)))^2/(4*k);
        inres = 1/(b - k*(gParam(1,i)-gParam(2,i)));
    end
    gParam(4,i) = b; %b_g (Izhekevich parameter b)
    gParam(5,i) = k; %k_g (Izhikevich parameter k)
end

for i = 1:granuleNum
    % Establish Izhikevich parameter c and ensure it is not above the
    % threshold potential
    while gParam(6,i) >= gParam(2,i)
        gParam(6,i) = normrnd(-75,75/10); %c_g (Izhikevich parameter c)
    end
end

gParam(7,:) = normrnd(1.2,1.2/10,1,granuleNum); %d_g (Izhikevich parameter d)
gParam(8,:) = normrnd(48,48/10,1,granuleNum); %cap_g (capacitance)

% GC synaptic parameters
tauA_g = normrnd(5.5,5.5/10,mitralNum,granuleNum); % tau_AMPA
tauNr_g = normrnd(10,10/10,mitralNum,granuleNum); % tau_NMDA_rise
tauNd_g = normrnd(80,80/10,mitralNum,granuleNum); % tau_NMDA_decay
gAMPA = normrnd(0.73,0.73/10,mitralNum, granuleNum); % g_AMPA
gNMDA = normrnd(0.84,0.84/10,mitralNum, granuleNum); % g_NMDA



%Find MCs belonging to each glomeruli and assign a respiratory input between
%0 and maximum input (for when odor is not present)
mResp_Amp = zeros(mitralNum,1);
mResp_Phase = zeros(mitralNum,1);
for i = 1:length(glomeruli)
    R = find(glomArray == glomeruli(i));
    mean_phase = rand*2*pi;

    % input strength
    meanInput = rand*150;
    mResp_Amp(R) = normrnd(meanInput,meanInput/5,length(R),1);
    
    % phase of input
    mResp_Phase(R) = normrnd(mean_phase,pi/4,length(R),1);
end

% number of odors you want to make
numOdors = 10;

% rounds of neurogenesis
numGeneses = 10;

% cell array to hold the spike trains (like a python list)
spikeTrainRecord = cell(numGeneses,numOdors);
gSpikes = cell(numGeneses,numOdors);
feedbackRecord = cell(1, numOdors);

% set the respiratory rates
fnorm = 2/1000; % 2 Hz breathing rate
fodor = 6/1000; % 6 Hz sniff rate during odor presentation

% file to write to keep track of how many trials have passed if you're
% impatient like me
fileID = fopen('trialResultsI.txt','w');

% Input currents and phases for odor presentation
odorAmpMatrix = zeros(mitralNum,numOdors);
odorPhaseMatrix = zeros(mitralNum, numOdors);

% glomeruli for your odors, where each row is a different odor
odorIndicesMatrix = zeros(numOdors, length(glomeruli));

% for setting granule cell
rmax = 300;

% Generate your odors
for i = 1:numOdors
    J = randperm(glomNum);
    odorGlomNum = 4 + randi(6);    
    for j = 1:odorGlomNum
        odorIndicesMatrix(i,J(j)) = 1;
    end
    [odorAmpMatrix(:,i), odorPhaseMatrix(:,i)] = odorGenerator(odorIndicesMatrix(i,:), glomArray, glomeruli, mitralNum);
end



% Generates a feedback pattern to GCs (ignore for now)
FB_frac = 0.25; % fraction of GCs targeted
FBnum = round(FB_frac*granuleNum); 
gFeedback1 = zeros(granuleNum,1);   
GC_order = randperm(granuleNum);
FB_cells = GC_order(1:FBnum); % select GCs 
gFeedback1(FB_cells) = rand(length(FB_cells),1) * 100; % set feedback strength to each GC

beginOdor = 0/fnorm; % adjust numerator for how much respiraation (i.e. no odor breaths) you want before odor input
endOdor = beginOdor + 2/fodor; % adjust numerator for how long you want the sniffing to go on - I recommend at least 2
beginFB = 0; % when do you want feedback to begin - ignore for now

addFrac = 0.1;
totalAddNum = mitralNum * 15 - granuleNum;
targetNum = round(addFrac * totalAddNum);

networkList = cell(1,numGeneses);

for gene = 1:numGeneses   
    networkList{gene} = network;
    gTotal = zeros(granuleNum,1);
    for odor = 1:numOdors
        % simulate without feedback
        [spikeTrainRecord{gene, odor}, ~, ~, gSpikes{gene,odor}]  = simulator(network, mParam, gParam,...
            mResp_Amp, mResp_Phase, fnorm, odorAmpMatrix(:,odor), odorPhaseMatrix(:,odor), fodor, zeros(granuleNum,1), beginOdor, beginFB, endOdor,...
            mGABA, tauG_m,...
            tauA_g, tauNr_g, tauNd_g, gAMPA, gNMDA);
        

        % simulate with feedback
    %     [spikeTrainRecord{2, trial}, ~, ~, ~]  = simulator(network, mParam, gParam,...
    %         mResp_Amp, mResp_Phase, fnorm, mOdor_Amp, mOdor_Phase, fodor, gFeedback1, beginOdor, beginFB, endOdor,...
    %         mGABA, tauG_m, mAMPA, tauA_m, mNMDA, tauNr_m, tauNd_m,...
    %         tauA_g, tauNr_g, tauNd_g, gAMPA, gNMDA,...
    %         tauAx_g, tauNrx_g, tauNdx_g, gAMPAx, gNMDAx);
        
        formatSpec = 'Trial is %d,%d\n';
        fprintf(fileID, formatSpec, gene, odor);
    end
    
    network = [network zeros(mitralNum,targetNum)];
    distance = [distance ones(mitralNum,targetNum)*-1];
    addedNum = 1;
    
    while addedNum <= targetNum
    
        % generate a new GC and assign properties
        newGranule = granule();
        newGranule = newGranule.assignProperties(rmax);

        % generate an empty matrix representing potential connections between
        % the MCs of the network and the new GC
        tempNet = zeros(mitralNum,1);
        tempDistance = ones(mitralNum,1)*-1;

        % for each MC, determine the probability of synapse with the new GC
        for j = 1:mitralNum

            % set the MC
            mitralCell = mitralArray(j);

            % calculate the total pre-existing synapses 
            deadSpace = sum(network(j,:));

            % determine the location of the center of the GC cone at the height
            % of the MC
            [g_loc_x, g_loc_y] = newGranule.calculateLocation(mitralCell.z);

            % calculate the radius of the GC cone at the height of the MC
            g_radius = newGranule.calculateRadius(mitralCell.z);

            % calculate the distance between the GC cone and MC at the height
            % of the MC
            s = norm([mitralCell.x, mitralCell.y] - [g_loc_x, g_loc_y]);

            % if the MC and GC are in range of one another, calculate the
            % average number of synapses
            if s < mitralCell.radius + g_radius && mitralCell.z > newGranule.z0 && mitralCell.z < newGranule.zmax

                % the average number of synapses for the MC/GC pair
                lambda = synProb(newGranule, mitralCell, deadSpace, g_radius, s);

                % calculate the probability of synapse assuming a Poisson
                % distribution
                prob = 1 - exp(-lambda);

                % sample the given probability and assign synaptic distance
                if rand < prob
                    tempNet(j,1) = 1;
                    while tempDistance(j,1) == -1
                       r = rand*mitralCell.radius;
                       theta = rand * 2 * pi;
                       x_r = r*cos(theta);
                       y_r = r*sin(theta);
                       if norm([x_r,y_r] - [s, 0]) < g_radius
                           tempDistance(j,1) = r;
                       end
                    end       
                end
            end
        end
    
        % calculate the total number of synapses the GC has made and
        % incorporate the new GC if it is greater than 0
        totalSynapses = sum(tempNet(:,1));

        if totalSynapses > 0
            % if the total number of synapses is greater than the number of
            % available spines, select a random subset of MCs equal to the
            % number of available spines to be incorporated into the network
            if totalSynapses > newGranule.availableSpines
               synmitrals = find(tempNet);
               synmitrals = synmitrals(randperm(length(synmitrals)));
               tempNet(:) = 0;
               tempNet(synmitrals(1:newGranule.availableSpines)) = 1;
               tempDistance(synmitrals(newGranule.availableSpines+1:totalSynapses)) = -1;
            end

            % assign the new connections to the appropriate entry in the
            % network for the GC
            network(:, addedNum) = tempNet;
            distance(:,addedNum) = tempDistance;

            % update the GC array
            addedNum = addedNum + 1;
        end
    end
    
    gParam = [gParam zeros(8,targetNum)];
    tauA_g = [tauA_g normrnd(5.5,5.5/10,mitralNum,targetNum)]; % tau_AMPA
    tauNr_g = [tauNr_g normrnd(10,10/10,mitralNum,targetNum)]; % tau_NMDA_rise
    tauNd_g = [tauNd_g normrnd(80,80/10,mitralNum,targetNum)]; % tau_NMDA_decay
    gAMPA = [gAMPA normrnd(0.73,0.73/10,mitralNum, targetNum)]; % g_AMPA
    gNMDA = [gNMDA normrnd(0.84,0.84/10,mitralNum, targetNum)]; % g_NMDA
    
    LmatT = repmat(L,[1,targetNum]); 
    mGABA = [mGABA normrnd(0.13,0.13/10,mitralNum,targetNum).* exp(-distance(:,granuleNum+1:granuleNum+targetNum)./LmatT)]; 
    tauG_m = [tauG_m normrnd(18,18/10,mitralNum, targetNum)];
    

        % Establish the intrinsic granule cell parameters
    for i = granuleNum+1:granuleNum+targetNum
        while gParam(1,i) >= gParam(2,i) || gParam(2,i) - gParam(1,i) < 20
            gParam(1,i) = normrnd(-71,71/10); % gVr (resting potential)
            gParam(2,i) = normrnd(-39,39/10); % gVt (threshold potential)
        end
    end
    gParam(3,granuleNum+1:granuleNum+targetNum) = normrnd(0.01,0.01/10,1,targetNum); %a_g (Izhikevich parameter a)

    for i = granuleNum+1:granuleNum+targetNum
        % Establish Izhikevich parameters b and k based on rheobase and input resistance
        rheobase = 0;
        inres = -1;
        while (rheobase < 10 || rheobase > 70 || inres < 0.25 || inres > 1.5) || b > 0
            b = normrnd(-2/15,2/10); 
            k = normrnd(1/15,1/10);
            rheobase = (b+k*(-gParam(1,i)+gParam(2,i)))^2/(4*k);
            inres = 1/(b - k*(gParam(1,i)-gParam(2,i)));
        end
        gParam(4,i) = b; %b_g (Izhekevich parameter b)
        gParam(5,i) = k; %k_g (Izhikevich parameter k)
    end

    for i = granuleNum+1:granuleNum+targetNum
        % Establish Izhikevich parameter c and ensure it is not above the
        % threshold potential
        while gParam(6,i) >= gParam(2,i)
            gParam(6,i) = normrnd(-75,75/10); %c_g (Izhikevich parameter c)
        end
    end

    gParam(7,granuleNum+1:granuleNum+targetNum) = normrnd(1.2,1.2/10,1,targetNum); %d_g (Izhikevich parameter d)
    gParam(8,granuleNum+1:granuleNum+targetNum) = normrnd(48,48/10,1,targetNum); %cap_g (capacitance)
    
    granuleNum = granuleNum + targetNum; %update total number of GCs in the network
    
    fname = sprintf('immortalResults.mat');
    save(fname,'spikeTrainRecord','odorAmpMatrix', 'odorPhaseMatrix', 'gSpikes','networkList','-v7.3');
end

    

    



%Generates the odor for the given number of glomeruli
function [mOdor_Amp, mOdor_Phase] = odorGenerator(odorIndices, glomArray, glomeruli, mitralNum)

    mOdor_Amp = zeros(mitralNum,1);
    mOdor_Phase = zeros(mitralNum,1);
    %Find MCs belonging to each glomeruli and assign an odor input between
    %0 and maximum inputs 
    for i = 1:length(glomeruli)
        mean_phase = 2*pi*rand;
        R = find(glomArray == i);
        if odorIndices(i) == 1

            % input strength
            meanInput = 150+450*rand; % you can make this a vector if you want to have different strengths for each glomerulus
            mOdor_Amp(R) = normrnd(meanInput,meanInput/5,length(R),1);

            % phase of input
            mOdor_Phase(R) = normrnd(mean_phase,pi/4,length(R),1);
        else
            % input strength
            meanInput = rand*150;
            mOdor_Amp(R) = normrnd(meanInput,meanInput/5,length(R),1);

            % phase of input
            mOdor_Phase(R) = normrnd(mean_phase,pi/4,length(R),1);
        end
    end
end




function [mSpikeTrain, mVolt, gVolt, gSpikes] = simulator(network, mParam, gParam,...
    mResp_Amp, mResp_Phase, fnorm, mOdor_Amp, mOdor_Phase, fodor, gFeedback, beginOdor, beginFB, endOdor,...
    mGABA, tauG_m,...
    tauA_g, tauNr_g, tauNd_g, gAMPA, gNMDA)

    % time step is 0.1 ms
    TS = 0.1; 
    
    % total record time proceeds to end of odor
    recordTime = endOdor;
    
    % set the time vector
    tspan = 0:TS:recordTime;   
    
    % Number of mitral cells and granule cells
    [mitralNum,granuleNum] = size(network);
    
    % set number of GCs to record for voltage traces
    granuleNumRecord = 1000;
    
    % unpack MC and GC parameters
    mVr = mParam(1,:); mVt = mParam(2,:);
    a_m = mParam(3,:); b_m = mParam(4,:); k_m = mParam(5,:); c_m = mParam(6,:); d_m = mParam(7,:);
    Cm = mParam(8,:);  k = mParam(9,:);
    k = repmat(k,[1,granuleNum]);       
    
    gVr = gParam(1,:); gVt = gParam(2,:);
    a_g = gParam(3,:); b_g = gParam(4,:); k_g = gParam(5,:); c_g = gParam(6,:); d_g = gParam(7,:);
    Cg = gParam(8,:); 
    
    % gating variable jump constant
    W = 0.5;
    
    % magnesium concentration
    mag = 1;

    % cells start at their resting potentials with 0 recoveery current
    mV = ones(1, mitralNum) .* mVr;
    mU = zeros(1,mitralNum);
    gV = ones(1, granuleNum) .* gVr;
    gU = zeros(1,granuleNum);
    
    % set variables to record voltages
    mVolt = zeros(mitralNum,length(tspan));
    gVolt = zeros(granuleNumRecord, length(tspan));
    
    % set synaptic gating variables 
    mIntegG = zeros(mitralNum, granuleNum);
    gIntegN = zeros(mitralNum, granuleNum);
    gIntegA = zeros(mitralNum, granuleNum);
    gX = zeros(mitralNum, granuleNum);
    
    % record MC and GC spikes
    mSpikeTrain = zeros(mitralNum,length(tspan));
    gSpikes = zeros(granuleNum,1);
    fileID = fopen('results7.txt','w');

    for t=1:length(tspan)     

        % find all MC and GCs which have spiked
        mFired = find(mV >= 30);     
        gFired = find(gV >= 25);  
        
        % record MC and GC spike times
        if t > 1667
            gSpikes(gFired) = gSpikes(gFired) + 1;
        end
        mSpikeTrain(mFired,t) = 1;
        
        % record voltage
        mVolt(:,t) = mV;
        
        % record GC voltage
        gVolt(:,t) = gV(1:granuleNumRecord);   

        % if in odor presentation, give odor input, otherwise give regular respiratory input 
        if tspan(t) > beginOdor && tspan(t) <= endOdor          
           mI = ((mOdor_Amp/2 + mOdor_Amp/4 .* (sin(2*pi*fodor*tspan(t) - mOdor_Phase)+1)));
        else
           mI = ((mResp_Amp/3 + mResp_Amp/3 .* (sin(2*pi*fnorm*tspan(t) - mResp_Phase)+1)));  
        end
        
        if tspan(t) > beginFB 
            gI = gFeedback;
        else
            gI = zeros(granuleNum,1);
        end

       
        % All synapses for firing MCs have NMDA and AMPA receptors updated
        gIntegA(mFired,:) = gIntegA(mFired,:) + network(mFired,:)*W.*(1-gIntegA(mFired,:));
        gX(mFired,:) = gX(mFired,:) + network(mFired,:)*W.*(1-gX(mFired,:));
        
        % All synapses for indirectly connected MCs have GABA receptors
        % updated
        for i = 1:length(mFired)
            gRel = find(network(mFired(i),:));
            mIntegG(:,gRel) = mIntegG(:,gRel) + k(:,gRel).*network(:,gRel)*W.*(1-mIntegG(:,gRel));
        end
        
        % All synapses for firing GCs have GABA receptors updated
        mIntegG(:,gFired) = mIntegG(:,gFired) + network(:,gFired)*W.*(1-mIntegG(:,gFired));
        
        % calculate external and synaptic currents 
        [mIreturn, gIreturn, mIntegG, gIntegA,gIntegN, gX] = current(mV, gV, mag,...
        mIntegG, gIntegA, gIntegN, gX,...
        mGABA, tauG_m, ...
        gAMPA, tauA_g, gNMDA, tauNr_g, tauNd_g, TS);
    
        mI = mI'+mIreturn;
        gI = gI'+gIreturn;
        
        % update Izhikevich values
        [mV, mU, gV, gU] = izhikevich(mV, mU, gV, gU, mI, gI, mFired, gFired, TS,...
                a_m, b_m, c_m, d_m, k_m, a_g, b_g, c_g, d_g, k_g, mVr, mVt, gVr, gVt, Cm, Cg);  
        
        formatSpec = 'Time is %1.2f\n';
        fprintf(fileID, formatSpec, tspan(t));
    end
end

% updates the currents based on the gating variables
function [mI, gI, mIG, gIA,gIN, Xg] = current(mV, gV, Mg,...
        mIntegG, gIntegA, gIntegN, gX,...
        mGABA, tauG_m, ...
        gAMPA, tauA_g, gNMDA, tauNr_g, tauNd_g, TS)
    
      
        % inhibitory reversal potential
        Ei = -70; 
        
        % excitatory reversal potential
        Ee = 0;
        
        % alpha values for synaptic and external NMDA receptors
        alphaG = 0.1;

        % calculate currents from gating variables
        Ig = -(mV-Ei) .* sum(mGABA .* mIntegG,2)';
        mI = Ig; 
        Ja = -(gV-Ee) .* sum(gAMPA.*gIntegA);
        Jn = -(gV-Ee) .* (sum(gNMDA.*gIntegN).* 1./(1+exp(-0.062*gV)*Mg/3.57));
        gI = Ja + Jn;       
        
        % update gating variables
        mIG = mIntegG ./ exp(TS./tauG_m);
        gIA = gIntegA ./ exp(TS./tauA_g);
        gIN = gIntegN + TS*(-gIntegN./tauNd_g + alphaG*gX.*(1-gIntegN));
        Xg = gX ./ exp(TS./tauNr_g); 
end     
        
        
% computes the Izhikevich equations
function [mitV, mitU, graV, graU] =  izhikevich(mV, mU, gV, gU, mI, gI, mFired, gFired, TS,...
                a_m, b_m, c_m, d_m, k_m, a_g, b_g, c_g, d_g, k_g, mVr, mVt, gVr, gVt, Cm, Cg)
        
        % adjust voltages and recovery currents for cells which have fired
        gV(gFired) = c_g(gFired);
        gU(gFired) = gU(gFired) + d_g(gFired);
        mV(mFired) = c_m(mFired);
        mU(mFired) = mU(mFired) + d_m(mFired);

        %Updating values for all mitral and granule cells based on
        %Izhikevich functions
        mV = mV + TS * (k_m./Cm .* (mV-mVr).*(mV-mVt) - mU./Cm + mI./Cm + 0.25*randn);
        mU = mU + TS * a_m.*(b_m .* (mV-mVr)-mU);

        gV = gV + TS * (k_g./Cg .* (gV-gVr).*(gV-gVt) - gU./Cg + gI./Cg + 0.25*randn);
        gU = gU + TS * a_g.*(b_g .* (gV-gVr)-gU);
     
        mitV = mV; mitU = mU;  
        graV = gV; graU = gU; 
end













