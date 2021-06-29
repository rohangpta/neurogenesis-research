% calculate the average number of synapses between a given MC and GC

function synapse = synProb(gr,mi,deadSpace, g_radius, s)

    % coefficient for converting dendrite length to volume of interaction
    q = 2.32;
    
    % average volume of a spine (Woolf 1991)
    spineVolume = 0.58;
    
    % calculate the spine density at the height of the MC
    density = gr.calculateDensity(mi.z); 

    % these functions represent the simplification of the inside
    % integral in the Methods section of the paper, when g(theta) and
    % g'(theta) are the limits of integration, respectively
    gfun = @(theta) mi.alpha ./(2.*pi).*atan(mi.k .* (s.*cos(theta)-sqrt(g_radius.^2-s.^2.*sin(theta).^2))-tan(mi.m));
    gfunprime = @(theta) mi.alpha ./(2.*pi).*atan(mi.k .* (s.*cos(theta)+sqrt(g_radius.^2-s.^2*sin(theta).^2))-tan(mi.m));

    % Below are the calculations for the expected number of synapses
    % depending on the distance between the MC and GC and the relative size
    % of the MC and GC radii
    
    % for MC radius greater than twice the GC radius
    if mi.radius > 2*g_radius

        if s < (g_radius + mi.radius) && s >= sqrt(g_radius^2 + mi.radius^2)
            mu = acos((mi.radius^2 + s^2 - g_radius^2) / (2*mi.radius*s)); 
            len = 2 *(mu *mi.alpha/(2*pi)* atan(mi.k*(mi.radius)-tan(mi.m)) - integral(gfun, 0, mu));
            
            % compute number of expected synapses adjusted for pre-existing
            % synapse volume (same for all conditions)
            synapse = density...
                     *(q*pi*len -spineVolume*deadSpace*(len/mi.totlength));


        elseif s >= (mi.radius - g_radius) && s < sqrt(g_radius^2 + mi.radius^2)
            mu = acos((mi.radius^2 + s^2 - g_radius^2) / (2*mi.radius*s));
            gamma = asin(g_radius/s);  
            len = 2 * (mu *mi.alpha/(2*pi)* atan(mi.k*(mi.radius)-tan(mi.m))...
                      - integral(gfun, 0, gamma) + integral(gfunprime, mu, gamma));
            synapse = density...
                      *(q*pi*len-spineVolume*deadSpace*(len/mi.totlength));

        elseif s >= g_radius && s < mi.radius - g_radius
            gamma = asin(g_radius/s);
            len = 2 * (integral(gfunprime,0,gamma) - integral(gfun,0,gamma));
            synapse = density...
                      *(q*pi*len -spineVolume*deadSpace*(len/mi.totlength)); 

        elseif s < g_radius
            len = (integral(gfunprime,0,2*pi)+mi.m*mi.alpha);
            synapse = density...
                      * (q*pi*len-spineVolume*deadSpace*(len/mi.totlength)); 

        else
            synapse = 0;
        end    

    % for MC radius greater than the GC radius but less than twice the GC radius    
    elseif mi.radius > g_radius

        if s < (g_radius + mi.radius) && s >= sqrt(g_radius^2 + mi.radius^2)
            mu = acos((mi.radius^2 + s^2 - g_radius^2) / (2*mi.radius*s)); 
            len = 2 *(mu *mi.alpha/(2*pi)* atan(mi.k*(mi.radius)-tan(mi.m)) - integral(gfun, 0, mu));
            synapse = density...
                     *(q*pi*len-spineVolume*deadSpace*(len/mi.totlength)); 

        elseif s >= g_radius && s < sqrt(g_radius^2 + mi.radius^2)
            mu = acos((mi.radius^2 + s^2 - g_radius^2) / (2*mi.radius*s));
            gamma = asin(g_radius/s);  
            len = 2 * (mu *mi.alpha/(2*pi)* atan(mi.k*(mi.radius)-tan(mi.m))...
                      - integral(gfun, 0, gamma) + integral(gfunprime, mu, gamma));
            synapse = density...
                      *(q*pi*len-spineVolume*deadSpace*(len/mi.totlength));

        elseif s >= (mi.radius - g_radius) && s < g_radius
            mu = acos((mi.radius^2 + s^2 - g_radius^2) / (2*mi.radius*s));
            len = 2 * (integral(gfunprime,mu,pi)+mi.m*mi.alpha/2 ...
                      +mu*mi.alpha/(2*pi)* atan(mi.k*(mi.radius)-tan(mi.m)));
            synapse = density...
                      * (q*pi*len-spineVolume*deadSpace*(len/mi.totlength));

        elseif s < (mi.radius - g_radius)
            len = (integral(gfunprime,0,2*pi)+mi.m*mi.alpha);
            synapse = density...
                      * (q*pi*len-spineVolume*deadSpace*(len/mi.totlength)); 

        else
            synapse = 0;
        end    

    % for MC radius less than the GC radius
    else    
        if s < (g_radius + mi.radius) && s >= sqrt(g_radius^2 + mi.radius^2)
            mu = acos((mi.radius^2 + s^2 - g_radius^2) / (2*mi.radius*s)); 
            len = 2 *(mu *mi.alpha/(2*pi)* atan(mi.k*(mi.radius)-tan(mi.m)) - integral(gfun, 0, mu));
            synapse = density...
                     *(q*pi*len -spineVolume*deadSpace*(len/mi.totlength));


        elseif s >= g_radius && s < sqrt(g_radius^2 + mi.radius^2)
            mu = acos((mi.radius^2 + s^2 - g_radius^2) / (2*mi.radius*s));
            gamma = asin(g_radius/s);  
            len = 2 * (mu *mi.alpha/(2*pi)* atan(mi.k*(mi.radius)-tan(mi.m))...
                      - integral(gfun, 0, gamma) + integral(gfunprime, mu, gamma));
            synapse = density...
                      *(q*pi*len-spineVolume*deadSpace*(len/mi.totlength));

        elseif s >= (g_radius-mi.radius) && s < g_radius
            mu = acos((mi.radius^2 + s^2 - g_radius^2) / (2*mi.radius*s));
            len = 2 * (integral(gfunprime,mu,pi)+ mi.m*mi.alpha/2 ...
                      +mu*mi.alpha/(2*pi)* atan(mi.k*(mi.radius)-tan(mi.m)));
            synapse = density...
                      * (q*pi*len-spineVolume*deadSpace*(len/mi.totlength));

        elseif s < (g_radius-mi.radius)
            len = mi.totlength;
            % since when integrating over the whole MC, len = mi.totlength
            synapse = density * (q*pi*len-spineVolume*deadSpace);
        else
            synapse = 0;
        end  
    end

    % if average number of synapses is less than 0, it is equal to 0
    if synapse < 0
        synapse = 0;
    end
end
