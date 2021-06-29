classdef granule %defines a granule cell
    properties
        x % x-location of bottom vertex
        y % y-location of bottom vertex
        z0 % z-location of bottom vertex
        top_x % x-location of top face
        top_y % y-location of top face
        zmax % z-location of top face
        rmax % radius of top face
        numSpines % total number of spines
        availableSpines % number of spines available for synapse
        top_shift % center point of top face
        theta % angle relative to 0 of line formed by bottom and top vertex centers
    end
    methods
        % Assign properties for a newly generated granule cell
        function g = assignProperties(g, space_rmax, neurogen_mode, mitrals)
            % Layers of the OB space
            EPLthickness = 131;
            MCLthickness = 36;
            IPLthickness = 27; 
            
            th = EPLthickness + MCLthickness + IPLthickness;
            
            if neurogen_mode == "None"
                % Set the location of the bottom vertex and constrain it to
                % fall within bounds (randomly)
                d = space_rmax+1;
                while d > space_rmax
                    g.x = -space_rmax + 2*rand*space_rmax;
                    g.y = -space_rmax + 2*rand*space_rmax;
                    d = norm([g.x,g.y]);
                end 
                
                % Set the z-position of the bottom vertex
                g.z0 = (MCLthickness + IPLthickness) * rand;
            
                % Set the z-position of the top face
                g.zmax = MCLthickness + IPLthickness + 1/2*EPLthickness + (1/2*EPLthickness*rand);   
            else
                % if doing neurogenesis, then find a target mitral cell
                % to place granule cell near to and place bottom vertex
                % as below
                load('centralities.mat', 'Eigenvector', 'Closeness', 'Betweenness');
                centrality_type = neurogen_mode;
                if centrality_type == "Eigenvector"
                    c = Eigenvector;
                elseif centrality_type == "Closeness"
                    c = Closeness;
                elseif centrality_type == "Betweenness"
                    c = Betweenness;
                end 

                % we set mu such that 80% of the samples fall in [0, 10], or in other words
                % 80% of the time we choose the top 11 mitral cells. I chose 80 because of
                % the observed '80-20' rule. upon solving, we get mu = ~6.21

                mu = 6.21;
                r = floor(exprnd(mu)) + 1;
                m = c(r) + 1;
                mc = mitrals(:, m);
                d = space_rmax+1;
                while d > space_rmax
                    g.x = -space_rmax/10 + mc.x + 2*rand*space_rmax/10;
                    g.y = -space_rmax/10 + mc.y + 2*rand*space_rmax/10;
                    d = norm([g.x,g.y]);
                end
                
                g.z0 = (MCLthickness + IPLthickness) * rand;
                g.zmax = MCLthickness + IPLthickness + 1/2*EPLthickness + (1/2*EPLthickness*rand); 
                % height ev is 0.668 th so if mc.z is near to 0.332th then
                % select a g.z0 near to it. Otherwise, choose one further 
                % so that gc can have 'full' height. 
                % g.z0 = max(mc.z - 2*rand*abs((0.332 * th) - mc.z), MCLthickness+IPLthickness);
                
                % EV of height is 1/2 (mcl + ipl) + 3/4 (epl)
                % so we do rand of th + 0.25 epl which has the same EV
                % cap it at the total height to keep within bounds
                % g.zmax = min(g.z0 + rand * th + 0.25 * EPLthickness, th);
                
                % this way, we have a decent probability of a new gc
                % attaching to a 'desirable' mc
            end
            
            % Set the location and radius of the top face and 
            % constrain the whole face to fall within bounds
            acceptable_top = false;            
            while acceptable_top == false

                % Set the xy distance of the center of the top face from
                % the bottom vertex
                g.top_shift = rand*50;
                g.theta = rand*2*pi;
                g.top_x = g.x + g.top_shift * cos(g.theta);
                g.top_y = g.y + g.top_shift * sin(g.theta);
                
                % Set the radius of the top face
                g.rmax = 0;
                while g.rmax < 25 || g.rmax > 170
                    g.rmax = normrnd(86,29);
                end
                
                % Check to see if the top face falls within bounds of the
                % system
                phi = atan2(g.top_y,g.top_x); 
                if norm([g.top_x + g.rmax*cos(phi), g.top_y + g.rmax*sin(phi)]) < space_rmax
                    acceptable_top = true;
                end
            end


            % determine the upper and lower bounds of the spine
            % distribution based on the volume
            volume = pi/3 * g.rmax^2 * (g.zmax -g.z0);
            upper_bound = 357.7 * atan(2.653/1000000 * volume);
            lower_bound = 39.31 * atan(1.043/100000*volume);
            g.numSpines = round(lower_bound + rand*(upper_bound-lower_bound));
            
            
            % determine the number of available spines for synapse
            zmin = MCLthickness+IPLthickness;
            densityFunction = @(z) (z-g.z0).*(z-g.zmax);
            g.availableSpines = round(-6*g.numSpines/((g.zmax-g.z0)^3)*integral(densityFunction, zmin, g.zmax));
        end
        
        % Calculate the radius of the GC cone 
        % at the height of the relevant MC
        function r = calculateRadius(g, miZ)
            r = g.rmax/(g.zmax - g.z0) * (miZ - g.z0);
        end
        
        % Calculate the spine density of the GC cone 
        % at the height of the relevant MC
        function d = calculateDensity(g, miZ)
            d = -6*g.numSpines/(pi*g.rmax^2*(g.zmax-g.z0)) * (miZ-g.zmax)/(miZ-g.z0);
        end
        
        % Calculate the location of the center of the GC spine distribution 
        % at the height of the relevant MC
        function [x,y] = calculateLocation(g, miZ)
            r = g.top_shift/(g.zmax-g.z0) * (miZ-g.z0);
            x = g.x + r*cos(g.theta);
            y = g.y + r*sin(g.theta);
        end
        
    end                                                     
end
