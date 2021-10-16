classdef airfoilClass < handle
    %airfoil - This class stores mid-points, end-points, cP graph, panel
    %angles, etc
    %   Detailed explanation goes here
    %%%%%%%%%%%%%%%%%%%%%%%%%% ANGLE PHI MAYBE NOT GETTING CALCULATED
    %%%%%%%%%%%%%%%%%%%%%%%%%% CORRECTLY
    properties
        Datapoints
        Collocation ={};
        Endpoints = {};
        NormPanels = {};
        NumPanels
        FoilName
        FigHandle
        Betas
   
        
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
            x_endpoints = obj.Endpoints.x{1};
            z_endpoints = obj.Endpoints.z{1};
            for j = 1:length(x_endpoints)-1
                dx = x_endpoints(j+1)-x_endpoints(j);
                dz = z_endpoints(j+1)-z_endpoints(j);
                theta(j) = atan2(dz,dx);
            end
        end
            
        function obj = airfoilClass(foilName,foilPoints)
            %airfoilClass Construct an instance of this class
            %   Detailed explanation goes here
            obj.FoilName = foilName;
            obj.Datapoints = foilPoints;
        end
        function obj = discretize(obj,numPanels)
            obj.NumPanels = numPanels;
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            xVec = obj.Datapoints(:,1);
            zVec = obj.Datapoints(:,2);
            panelSize = ceil(length(xVec)/numPanels);
            count = 1;
            endPoints.x = {};
            endPoints.z = {};
            
            for i = 1:panelSize:length(xVec)
                endPoints.x{1}(count) = xVec(i);
                endPoints.z{1}(count) = zVec(i);
                count = count+1;
            end
            endPoints.x{1}(end+1) = xVec(end);
            endPoints.z{1}(end+1) = zVec(end);
            obj.Endpoints.x = endPoints.x;
            obj.Endpoints.z = endPoints.z;
            
            for j = 1:length(endPoints.x{1})-1
                obj.Collocation.x(j) = mean([obj.Endpoints.x{1}(j) obj.Endpoints.x{1}(j+1)]);
                obj.Collocation.z(j) = mean([obj.Endpoints.z{1}(j) obj.Endpoints.z{1}(j+1)]);
                            
                dx = (obj.Endpoints.x{1}(j+1)-obj.Endpoints.x{1}(j));
                dz = (obj.Endpoints.z{1}(j+1)-obj.Endpoints.z{1}(j));
                x = obj.Collocation.x(j);
                normal = (dx/dz)*x+obj.Collocation.z(j);
                normVecx{j} = [obj.Collocation.x(j) obj.Collocation.x(j)+dz/3];
                normVecz{j} = [obj.Collocation.z(j) obj.Collocation.z(j)-dx/3];
                
                obj.NormPanels.x(j,:) = normVecx{j};
                obj.NormPanels.z(j,:) = normVecz{j};
                
                phi = rad2deg(atan2((obj.Endpoints.z{1}(j+1)-obj.Endpoints.z{1}(j)),...
                    (obj.Endpoints.x{1}(j+1)-obj.Endpoints.z{1}(j))));

            end
            

                
            status = 'Discretization successful';
            
        end
        function beta = calcBetas(obj,velocityDir,plotVel)
            normVec = obj.NormPanels;
            CO = obj.Collocation;

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
        
        function fig = plotPanels(obj)
            fig = figure();
            plot(obj.Datapoints(:,1),obj.Datapoints(:,2),'r-','linewidth',3),hold on
            plot(obj.Endpoints.x{1},obj.Endpoints.z{1},'g-','linewidth',2)
            plot(obj.Collocation.x,obj.Collocation.z,'ko',...
                'markerfacecolor','k','markersize',2)

            title(obj.FoilName)
            legend('Airfoil Points','Discretized Panels','Collocation','location',...
                'best');
            grid on
            axis equal
            fill(obj.Endpoints.x{1},obj.Endpoints.z{1},'g')
            obj.FigHandle = fig;
        end
        
        function plotNormals(obj)
            fig = obj.FigHandle;
            figure(fig)
            for i = 1:length(obj.NormPanels.x)
            plot(obj.NormPanels.x(i,:),obj.NormPanels.z(i,:),'k')
            legend off
            end
        end
        
            

    end
end

