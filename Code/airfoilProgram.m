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

if exist('xyAirfoilData','var') ~= 1
    
%% Get Airfoil .dat files
    foilPath = '../Data Files';
% fileList = dir([foilPath '/*.dat'])
    file = uigetfile([foilPath '/*.dat']);

%% Load File

    opts = detectImportOptions([foilPath '/' file]);
    xyAirfoilData = readmatrix([foilPath '/' file],'NumHeaderLines',1);
end


%% Create Panels and Plot
numPanels = 60;
velDir = [1 0];
NACA = airfoilClass(file,xyAirfoilData);
NACA.discretize(numPanels);
plot = NACA.plotPanels;
NACA.plotNormals;


%% Solve Vortex Panel Method
gamma = 1;
[u,w] = VOR2DC(gamma,NACA);

function [u,w] = VOR2DC(gamma,airfoilObj)
coll_x = airfoilObj.Collocation.x;
coll_y = airfoilObj.Collocation.x;

theta1 = atan2((z-z2),(x-x2));
theta2 = atan2((z-z1),(x-x1));
R2 = (x-x2)^2 + (z-z2)^2;
R1 = (x-x1)^2+(z-z1)^2;

for i = 1:length(coll_x)
    for j = 1:length(coll_x)
        up = (gamma/(2*pi))*(theta2-theta1);
        wp = (-gamma/(4*pi))*log(R1/R2);
        
    end
end

end


