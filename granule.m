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
        function g = assignProperties(g, space_rmax)
            % Layers of the OB space
            EPLthickness = 131;
            MCLthickness = 36;
            IPLthickness = 27; 
            
            
            % Set the location of the bottom vertex and constrain it to
            % fall within bounds
            d = space_rmax+1;
            while d > space_rmax
                g.x = -space_rmax + 2*rand*space_rmax;
                g.y = -space_rmax + 2*rand*space_rmax;
                d = norm([g.x,g.y]);
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

            % Set the z-position of the bottom vertex
            g.z0 = (MCLthickness + IPLthickness) * rand;
            
            % Set the z-position of the top face
            g.zmax = MCLthickness + IPLthickness + 1/2*EPLthickness + (1/2*EPLthickness*rand);    

            
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
