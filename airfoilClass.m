classdef airfoilClass < handle
    %airfoil - This class stores mid-points, end-points, cP graph, panel
    %angles, etc
    %   Detailed explanation goes here
    %%%%%%%%%%%%%%%%%%%%%%%%%% ANGLE PHI MAYBE NOT GETTING CALCULATED
    %%%%%%%%%%%%%%%%%%%%%%%%%% CORRECTLY
    properties
        xz = struct('Coords',{},'Collocation',{},'Endpoints',{},'NormPanels',{});
        panelCoords  = struct('Coords',{},'Collocation',{},...
            'Endpoints',{},'NormPanels',{});
        NumPanels
        FoilName
        FigHandle
        Betas
        RefLength
        Phi
        PanelSize   
    end
    properties (Dependent)
        Phis
        BetaDeg
        PanelAngles
        
    end
    
    methods
        function a = get.Phis(obj)
            a = obj.Betas - pi/2;
        end
        function deg = get.BetaDeg(obj)
            deg = rad2deg(obj.Betas);
        end
        function theta = get.PanelAngles(obj)
            x_endpoints = obj.xz.Endpoints.x;
            z_endpoints = obj.xz.Endpoints.z;
            for j = 1:length(x_endpoints)-1
                dx = x_endpoints(j+1)-x_endpoints(j);
                dz = z_endpoints(j+1)-z_endpoints(j);
                theta(j) = atan2(dz,dx);
            end
        end    
        % Rotate panel by alpha radians
        function coords = panel2xz(obj,alpha)
            xz_coords = obj.xz.Coords;
            rotMat = [cos(alpha) sin(alpha);
                -sin(alpha) cos(alpha)];
            coords = rotMat*xz_coords'; 
            obj.panelCoords(1).Coords = coords';
            obj.panelCoords.Collocation = rotMat*[obj.xz.Collocation.x;...
                obj.xz.Collocation.z];
            obj.panelCoords.Endpoints = rotMat*[obj.xz.Endpoints.x;...
                obj.xz.Endpoints.z];
            obj.panelCoords.NormPanels = rotMat*[obj.xz.NormPanels.x;...
                obj.xz.NormPanels.z]';
        end
        % Not really working correctly
        function coords = xz2panel(obj,alpha)
            xz_coords = obj.xz.Coords;
            rotMat = [cos(alpha) -sin(alpha);
                sin(alpha) cos(alpha)];
            
            for i = 1:length(xz_coords)-1
                xVec(i) = xz_coords(i+1,1)- xz_coords(i,1);
                zVec(i) = xz_coords(i+1,2)- xz_coords(i,2);
            end
            
            coords = rotMat*[xVec;zVec]; 
            coords = coords';
%             obj.panelCoords(1).Coords = coords';
%             obj.panelCoords.Collocation = rotMat*[obj.xz.Collocation.x;...
%                 obj.xz.Collocation.z];
%             obj.panelCoords.Endpoints = rotMat*[obj.xz.Endpoints.x{:};...
%                 obj.xz.Endpoints.z{:}];
%             obj.panelCoords.NormPanels = rotMat*[obj.xz.NormPanels.x;...
%                 obj.xz.NormPanels.z]';
        end
        
        % Initialization good
        function obj = airfoilClass(foilName,foilPoints,reflength)
            %airfoilClass Construct an instance of this class
            %   Detailed explanation goes here
            obj.FoilName = foilName;
            obj.RefLength = reflength;
            foilpointsClockwise = zeros(length(foilPoints),1);
            foilpointsClockwise(:,1) = flipud(foilPoints(:,1));
            foilpointsClockwise(:,2) = flipud(foilPoints(:,2));

            obj.xz(1).Coords = foilpointsClockwise;      
        end
 
        % Discretization good
        function obj = discretize(obj,numPanels)
            obj.NumPanels = numPanels;
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if length(obj.xz.Coords(:,1)) < obj.NumPanels
                fprintf(['Number of panels requested (',num2str(obj.NumPanels),...
                    ') is greater than\nthe number of airfoil nodes available (',...
                    num2str(length(obj.xz.Coords)),').\n\n']);
                fprintf(['Mesh will be initialized using the max number of panels (',...
                    num2str(length(obj.xz.Coords)),').\n\n']);
            end
            xVec = obj.xz.Coords(:,1);
            zVec = obj.xz.Coords(:,2);
            panelSize = ceil(length(xVec)/numPanels);
            count = 1;
            Endpoints.x = [];
            Endpoints.z = [];
            
            for i = 1:panelSize:length(xVec)-1
                Endpoints.x(count) = xVec(i);
                Endpoints.z(count) = zVec(i);
                count = count+1;
            end
            Endpoints.x(end+1) = xVec(end);
            Endpoints.z(end+1) = zVec(end);
            obj.xz.Endpoints.x = Endpoints.x;
            obj.xz.Endpoints.z = Endpoints.z;
            
            for j = 1:length(Endpoints.x)-1
                obj.xz.Collocation.x(j) = mean([obj.xz.Endpoints.x(j) obj.xz.Endpoints.x(j+1)]);
                obj.xz.Collocation.z(j) = mean([obj.xz.Endpoints.z(j) obj.xz.Endpoints.z(j+1)]);
                            
                dx = (obj.xz.Endpoints.x(j+1)-obj.xz.Endpoints.x(j));
                dz = (obj.xz.Endpoints.z(j+1)-obj.xz.Endpoints.z(j));
                obj.PanelSize(j) = (dx^2+dz^2)^.5;
                x = obj.xz.Collocation.x(j);
                normal = (dx/dz)*x+obj.xz.Collocation.z(j);
                normVecx{j} = [obj.xz.Collocation.x(j) obj.xz.Collocation.x(j)-dz/3];
                normVecz{j} = [obj.xz.Collocation.z(j) obj.xz.Collocation.z(j)+dx/3];
                
                obj.xz.NormPanels.x(j,:) = normVecx{j};
                obj.xz.NormPanels.z(j,:) = normVecz{j};
                
                phi(j) = atan2(dz,dx);
                if (phi(j) < 0)                                                        % Make all panel angles positive [deg]
                    phi(j) = phi(j) + 2*pi;
                end
                obj.Phi(j) = phi(j);
            end
            

                
            status = 'Discretization successful';
            
        end
        function beta = calcBetas(obj,velocityDir,plotVel)
            normVec = obj.NormPanels;
            CO = obj.xz.Collocation;

            if plotVel == 1
                fig = obj.FigHandle;
                figure(fig)
            end
            velocityDir = velocityDir/norm(velocityDir);
            for i = 1:length(CO.x)
                origin = [normVec.x(i,1) normVec.z(i,1) 0];
                vec1 = [normVec.x(i,2) normVec.z(i,2) 0] - origin;
                vec1 = vec1/norm(vec1);
                vec2 = [velocityDir 0];
%                 vec2 = origin - [velocityDir 0];
                if plotVel == 1
                    plot([origin(1) origin(1)-velocityDir(1)/30],...
                        [origin(2) origin(2)+velocityDir(2)/30],'b',...
                        'linewidth',1)
                end
                angle(i) = atan2(norm(cross(vec1,vec2)),dot(vec1,vec2));
                if angle(i)>pi 
                    angle(i)=angle(i)-2*pi; 
                end
                if angle(i)<0
                    angle(i)=angle(i)+2+pi;
                end
%                 
%                 beta(i) = acos(dot(velocityDir,normVec));
%                 beta(i) = atan2d(x1*y2-y1*x2,x1*x2+y1*y2);
            end
            obj.Betas = angle;
        end
        
        % Plot works
        function fig = plotPanels(obj,plotKey,varargin)
           
%             fig = figure();
            ref_length = obj.RefLength;
            if strcmp(plotKey,'xy') == 1
                plot(obj.xz.Coords(:,1)./ref_length,obj.xz.Coords(:,2)./ref_length,'r-','linewidth',3),hold on
                plot(obj.xz.Endpoints.x./ref_length,obj.xz.Endpoints.z./ref_length,'g-','linewidth',2)
                plot(obj.xz.Collocation.x./ref_length,obj.xz.Collocation.z./ref_length,'ko',...
                    'markerfacecolor','k','markersize',2)

                title(obj.FoilName)
                legend('Airfoil Points','Discretized Panels','Collocation','location',...
                    'best');
                grid on
                if size(varargin,1) >= 1
                    if strcmp(varargin{1},'equal') == 1
                        axis equal
                    end
                end
                fill(obj.xz.Endpoints.x./ref_length,obj.xz.Endpoints.z./ref_length,'g')
%                 obj.FigHandle = fig;
            end
            
            if strcmp(plotKey,'panelCoords') == 1
                plot(obj.panelCoords.Coords(:,1),obj.panelCoords.Coords(:,2),'r-','linewidth',3),hold on
                plot(obj.panelCoords.Endpoints(1,:),obj.panelCoords.Endpoints(2,:),'g-','linewidth',2)
                plot(obj.panelCoords.Collocation(1,:),obj.panelCoords.Collocation(2,:),'ko',...
                    'markerfacecolor','k','markersize',2)
                title(obj.FoilName)
                legend('Airfoil Points','Discretized Panels','Collocation','location',...
                    'best');
                grid on
                if size(varargin,1) >= 1
                    if strcmp(varargin{1},'equal') == 1
                        axis equal
                    end
                end
%                 axis equal
                fill(obj.panelCoords.Endpoints(1,:),obj.panelCoords.Endpoints(2,:),'g')
%                 obj.FigHandle = fig;
            end
        end
        
        % Plot normals works
        function plotNormals(obj)
%             fig = obj.FigHandle;
%             figure(fig)
            for i = 1:length(obj.xz.NormPanels.x)
            plot(obj.xz.NormPanels.x(i,:)./obj.RefLength,obj.xz.NormPanels.z(i,:)./obj.RefLength,'k')
            legend off
            end
        end
        

    end
end

