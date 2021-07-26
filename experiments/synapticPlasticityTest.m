load('simData/gAMPA.mat');
load('simData/gNMDA.mat');
load('simData/gSpikes.mat');
load('simData/mSpikes.mat');
load('data69/fullNetworkControl.mat');
load('simData/mSpikeTrain');

mFired = mSpikeTrain;

mFired = mFired(:, 501);

%disp(gAMPA(15, find(network(15, :))));

[ga, gn] = synapticPlasticity(mFired, gAMPA, gNMDA, network, mSpikes, gSpikes);

%disp(ga(15, find(network(15, :))));

function [gAMPA, gNMDA] = synapticPlasticity(mFired, gAMPA, gNMDA,... 
network, mSpikes, gSpikes)

    % find the MCs and GCs corresponding to the synapsesa
    MCs = find(mFired);
    for i = 1:length(MCs)
        mc = MCs(i);
        GCs = find(network(mc, :));
        for m = 1:length(GCs)
            gc = GCs(m);
            sum1 = 0;
            sum2 = 0;
            for j = 1:length(mSpikes{mc})
                for k = 1:length(gSpikes{gc})
                    diff = (gSpikes{gc}(k) - mSpikes{mc}(j))/1000;

                    A1 = gAMPA(mc, gc) * 0.05;
                    A2 = gNMDA(mc, gc) * 0.05;
                    % arbitrary params to 'slow down' hebbian updates
                    tau = 0.1;
                    if diff > 0
                        W = @(x) exp(-x/tau);
                    else
                        W = @(x) -exp(x/tau);
                    end
                    sum1 = sum1 + A1 * W(diff);
                    sum2 = sum2 + A2 * W(diff);
                end
            end
            gAMPA(mc, gc) = gAMPA(mc, gc) + sum1;
            gNMDA(mc, gc) = gNMDA(mc, gc) + sum2;
        end
    end

end