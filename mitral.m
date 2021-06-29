classdef mitral %defines a mitral cell
    properties
        x % x-location of center of disk
        y % y-location of center of disk
        z % z-location of disk
        radius % radius of lateral dendrite field
        glomerulus % identity of glomerulus to which MC belongs
        totlength % total length of dendrite available to the MC
        type % MC cell type
        % Parameters for lateral dendrite density equation
        alpha
        k
        m
        number % MC "name tag"
    end
    methods
        function m = assignProperties(m, gl, glomX, glomY, rmax, num)
            % Thicknesses of the OB space layers 
            EPLthickness = 131;
            MCLthickness = 36;
            IPLthickness = 27; 
            
            % assign glomerulus
            m.glomerulus = gl;
            
            % randomly assign MC as type I or type II
            if rand < (2/3)
                m.type = 1;
            else
                m.type = 2;
            end
            
            % set MC radius
            m.radius = 75 + rand*(800-75);
            
            % set MC height based on type
            if m.type == 1
                m.z = IPLthickness + MCLthickness + (EPLthickness * 0.5 * rand);         
            elseif m.type == 2
                m.z = IPLthickness + MCLthickness + EPLthickness * 0.4 + (EPLthickness * 0.4 * rand);
            end
            
            % set total length of dendrite 
            factor = 2*0.003827/3 + rand*0.003827*(2/3);
            m.totlength = factor * pi * m.radius^2;

            % set variables gamma and xi to determine alpha, m, and k
            gamma = 0.2 + rand*0.1;
            xi = 1/3 + (4/5-1/3)*rand;
            m.k = 1/(gamma*m.radius)*sqrt(1/xi-1);
            m.m = atan(sqrt(1/xi-1));
            m.alpha = m.totlength/(atan(m.k*m.radius-tan(m.m))+m.m);
            
            % create distribution of distance from glomerulus
            glom_distance = makedist('logistic','mu', 78.3618, 'sigma', 23.0732);
            glom_distance = truncate(glom_distance, 0, 300);

            % determine MC location based on glomerulus location
            location = rmax+1;
            while location > rmax
                r = random(glom_distance,1,1);
                theta = rand * 2 * pi;
                m.x = glomX + r*cos(theta);
                m.y = glomY + r*sin(theta);
                location = sqrt(m.x^2 + m.y^2);
            end 
            
            % assign label
            m.number = num;
       end
    end
end

