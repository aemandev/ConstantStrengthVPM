%{
Write a program to calculate the lift coefficient and pitching moment coefficient, 
as well as plot the pressure coefficient (cp) distribution for an airfoil using 
the constant strength vortex panel method. Compare the results from your program 
to the results from XFLR5's INVISCID analysis of the same airfoil. Discuss your results.

The program should be able to read in an input file which describes 
the airfoil geometry in XFLR5 format.

You should be able to change the angle of attack and your program should 
work for reasonable angles of attack (45 degrees is NOT reasonable. 
15 degrees probably is, but it depends on your airfoil shape).

Your program should work for any reasonable airfoil shape. 
Some good airfoils to test are the NACA 2412 (naca2412.dat), 
and the S1223 airfoil (s1223.dat).

Note that for the constant strength vortex method, the problem is overconstained 
so you will have to ignore the zero normal flow condition for one panel. 
You should try selecting different panels to ignore. 
Compare and discuss your results.

Once your program works for the above requirements, you can try enhancing 
your program to do any (or all) of the following for extra credit:

Modify your program to use the linear strength vortex panel method.

Calculate drag. This should be close to zero for inviscid flow, but 
calculating this will give you a better idea of how much error may be 
caused by numerical approximation.

Calculate lift and pitching moment forces instead of dimensionless 
coefficients.

Alternative methods: If you can't get the constant strength vortex panel 
method to work, you can use a different method instead with a slight grade 
penalty. Other methods include point vortices (one vortex per panel instead 
of a vortex sheet), Smith-Hess method (source panels plus vortices), or 
constant strength source panels (note that this method won't produce lift 
but you can get a cp distribution).
%}
close all; clc
clear CP A B DL TH VEL b EP CO

if exist('xyAirfoilData','var') ~= 1
    
%% Get Airfoil .dat files
    foilPath = './Data Files';
% fileList = dir([foilPath '/*.dat'])
    file = uigetfile([foilPath '/*.dat']);

%% Load File

    opts = detectImportOptions([foilPath '/' file]);
    xyAirfoilData = readmatrix([foilPath '/' file],'NumHeaderLines',1);
end

%% Create Panels and Plot
ref_length = 1;
NACA = airfoilClass(file,xyAirfoilData,ref_length);
numPanels = 200;
NACA.discretize(numPanels);

% Plot initial airfoil
figure()
NACA.plotPanels('xy','equal');
NACA.plotNormals;

% Angle of attack
alpha = 10;

% Calculate CP, CL, M (pitching moment), and velocity
[CP,CL,M,VEL] = VP(NACA,alpha);

% Plot CP
figure()
subplot(2,1,1)

% Plot the rotated airfoil by *alpha* radians
NACA.panel2xz(deg2rad(alpha));
NACA.plotPanels('panelCoords');



% NACA.plotNormals
subplot(2,1,2)
NACA.plotPanels('panelCoords');

% Plot CP values 
plot(NACA.xz.Collocation.x,-CP,'bo'),grid on, grid minor
% set(gca, 'YDir','reverse')
% axis equal


function [CP,CL,M0,VEL] = VP(airfoilObj,alpha)
AL = deg2rad(alpha);
M = length(airfoilObj.xz.Collocation.x);
for i = 1:M
    for j = 1:M

        X2T = airfoilObj.xz.Endpoints.x(j+1)-airfoilObj.xz.Endpoints.x(j);
        Z2T = airfoilObj.xz.Endpoints.z(j+1)-airfoilObj.xz.Endpoints.z(j);
        XT = airfoilObj.xz.Collocation.x(i)-airfoilObj.xz.Endpoints.x(j);
        ZT = airfoilObj.xz.Collocation.z(i)-airfoilObj.xz.Endpoints.z(j);
 
        X2 = X2T*cos(airfoilObj.Phi(j))+Z2T*sin(airfoilObj.Phi(j));
        Z2 = 0;
        X = XT*cos(airfoilObj.Phi(j))+ZT*sin(airfoilObj.Phi(j));
        Z = -XT*sin(airfoilObj.Phi(j))+ZT*cos(airfoilObj.Phi(j));

        if i == 1
            DL(j) = X2;
        end
        
        if i == j
            up = 0.5;
            wp = 0;
        else
            [up,wp] = VOR2DC(X,Z,X2,Z2);
        end

        U = up*cos(-airfoilObj.Phi(j))+wp*sin(-airfoilObj.Phi(j));
        W = -up*sin(-airfoilObj.Phi(j))+wp*cos(-airfoilObj.Phi(j));
        
        A(i,j) = -U*sin(airfoilObj.Phi(i))+W*cos(airfoilObj.Phi(i));
        B(i,j) = U*cos(airfoilObj.Phi(i))+W*sin(airfoilObj.Phi(i));
    end

%     b(i)= sin(airfoilObj.Phi(i)-AL);
    b(i) = cos(AL)*sin(airfoilObj.Phi(i))-sin(AL)*cos(airfoilObj.Phi(i));
    
end

% Apply Kutta B.C.
b = b';
A(M-1,:) = 0;
A(M-1,M-1) = 1;
A(M-1,1) = 1;
b(M-1) = 0;

% Solve for Gamma
G = A\b;

% Calculate V and Cp
CL = 0;
M0 = 0;
Vinf = 1;
for i = 2:M
    TEMP = 0;
    for j = 1:M
        TEMP=TEMP+B(i,j)*G(j);
    end
    VEL(i)=TEMP+cos(AL)*cos(airfoilObj.Phi(i))+sin(AL)*sin(airfoilObj.Phi(i));
    CL=CL+VEL(i)*airfoilObj.PanelSize(i);
    CP(i) = 1-((VEL(i)+VEL(i-1))/2/Vinf)^2;
    
% %     M0 = M0 + VEL(i)*DL(i)*airfoilObj.xz.Collocation.x(i)*cos(airfoilObj.Phi(i));
    
end
M0 = sum(CP.*airfoilObj.xz.Collocation.x.*airfoilObj.PanelSize.*cos(airfoilObj.Phi));
fprintf(['Results with ',num2str(length(airfoilObj.Phi)),' panels at ',num2str(alpha),...
    ' degrees angle of ','attack:\nCL = ',num2str(CL),'\n'])
end

function [up,wp] = VOR2DC(X,Z,X2,Z2)
% 
% R1 = X^2+Z^2;
% R2 = (X-X2)^2+Z^2;
% TH1 = atan2(Z,X);
% TH2 = atan2(Z-Z2,X-X2);

R1 = sqrt(X^2+Z^2);
R2 = sqrt((X-X2)^2+Z^2);
TH1 = atan2(Z,X);
TH2 = atan2(Z,X-X2);


up = (1/(2*pi))*(TH2-TH1);
wp = (1/(2*pi))*log(R2/R1);


end



