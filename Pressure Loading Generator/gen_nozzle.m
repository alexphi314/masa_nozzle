%Generate Nozzle Geometry based on piecewise sections
%Make use of premade geometry functions:
%gen_line(Dy,slope,res)-->          Generates a linear section
%gen_circ(R,theta_i,theta_f,res)--> Generates a cicular section
%gen_bell(Dy, a, b, res)-->         Generates a parabolic section

%Script fits each sectiomn together automatically. Each section should be
%described as height (distance from nozzle's center axis to the inner wall)
%as a function of y (the distance from the begining of the converging
%section)

%Example for a converging section with an angle of -70deg from y axis
% N = 1; <---Section Number
% Dy = 0.00635; <---Length in [m]
% slope = -tand(70); <---dH/dy
% res = 10; %<--- Number of steps in a section
% [sect_ls(N).L, sect_ls(N).H] = gen_line(Dy, slope, res); <--- Propper Assignment fromat for prebuilt functions

clear;
close all;
addpath('./Supporting_Functions/');
res = 25;

%% Mojave NOZZLE v1 %%
% Naming convention from RPA output %
H1 = 0.078200; %The initial height of the nozzle in [m]
len_conv = 30.5294e-3; %m %length of section 1
R1 = 0.049935; %fillet radius leading to throat
rt = 0.033045; %m
theta_n = 23.72; %deg
Rn = 0.012720; %m
% Rao Bell contour coefficients
a = 0.750332;
b = -0.128206;
c = -0.098494;
d = 0.0256073;
x41 = 0.00531717;
x42 = 0.193768;

%SECTION 1
N = 1;
slope = -1;
[sect_ls(N).L, sect_ls(N).H] = gen_line(len_conv,slope,res);

%SECTION #2
N = 2;
theta_f = 270;
theta_i = 225;
[sect_ls(N).L, sect_ls(N).H] = gen_circ(R1,theta_i,theta_f,res); %Length vector

%SECTION #3
N = 3;
theta_i = 270;
theta_f = theta_i+theta_n;
[sect_ls(N).L, sect_ls(N).H] = gen_circ(Rn, theta_i, theta_f, res);

%SECTION #4
N = 4;
[sect_ls(N).L, sect_ls(N).H] = gen_bell(x41, x42, a, b, c, d, res);

%% Formats Output %% DO NOT CHANGE
for i = [1:size(sect_ls,2)-1]
    disp(min(sect_ls(i).H)*39.3701);
    sect_ls(i+1).L = sect_ls(i+1).L + sect_ls(i).L(end);
    sect_ls(i).L = sect_ls(i).L(1:end-1);
    sect_ls(i+1).H = sect_ls(i+1).H + sect_ls(i).H(end);
    sect_ls(i).H = sect_ls(i).H(1:end-1);
end
y = [];
H = [];
for i = [1:size(sect_ls,2)]
    y = cat(2,y,sect_ls(i).L);
    H = cat(2,H,sect_ls(i).H);
end
H = H + H1;



plot(y*39.3701,H*39.3701);
xlabel('in');
ylabel('in');
axis equal
A = pi.*H.^2;
geometry = [y.',A.'];
save('nozzle_geometry.txt','geometry','-tabs','-ascii');