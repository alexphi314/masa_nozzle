function NozzleGenerator
clc
clear all
close all

%% Initial Conditions
prompt = {'Throat Radius m','Exit Radius m','Chamber Radius m'};
title = 'Nozzle Paramaters';
definput = {'.033045','.097945','.0782'};
Radii = inputdlg(prompt,title,[1 40],definput);

rt = str2double(Radii{1});
re = str2double(Radii{2});
rc = str2double(Radii{3});

prompt = {'/theta inflection deg','/theta exit deg'};
title = 'Nozzle Angles';
Angles = inputdlg(prompt,title,[1 40]);

thetaIn = str2double(Angles{1});
thetaEx = str2double(Angles{2});
% At = 0.672;   %Area of nozzle throat (m^2)
% AR = 16;      %Nozzle expansion ratio
% Ac = 4*At;    %Combustion chamber area (m^2)
% thetaIn = 27; %Nozzle Inflection Angle (deg)
% thetaEx = 10; %Nozzle Exit Angle (deg)

%% PreCalculations
% 
% Ae = At*AR;
% 
% rc = sqrt(Ac/pi());
% rt = sqrt(At/pi());
% re = sqrt(Ae/pi());
% 
% fprintf("Combustion Chamber Radius = %g m\n",rc);
% fprintf("Nozzle Throat Radius = %g m\n",rt);
% fprintf("Nozzle Exit Radius = %g m\n\n",re);

dr = rc - rt;
thetaC = acosd((1.5*rt - dr)/(1.5*rt))
L100 = (re-rt)/tand(15);
L80 = 0.8*L100;

%% Establish Domains
%Region 1, Pre-Throat Rao Circle 1.5rt, NOTE: x=0 marks throat location
domain1 = linspace(-1.5*rt*sind(thetaC),0,100);
field1 = [];

%Region 2, Post-Throat Rao Circle 0.4rt
domain2 = linspace(0,0.4*rt*sind(thetaIn),100);
field2 = [];

%Region 3, Post-Throat Rao Parabola with specified inflection and exit
domain3 = linspace(0.4*rt*sind(thetaIn),L80,100);
field3 = [];

%% System of Equations

c1 = sym('c1');
c2 = sym('c2');
c3 = sym('c3');
c4 = sym('c4');

xi = 0.4*rt*sind(thetaIn);
yi = f2(xi);
dyi = sind(thetaIn)/cosd(thetaIn);

xe = L80;
ye = re;
dye = sind(thetaEx)/cosd(thetaEx);

eqn1 = c1*xi + c2 + (c3*xi + c4)^(1/2) == yi;
eqn2 = c1 + (1/2)*c3*((c3*xi + c4)^(-1/2)) == dyi;
eqn3 = c1*xe + c2 + (c3*xe + c4)^(1/2) == ye;
eqn4 = c1 + (1/2)*c3*((c3*xe + c4)^(-1/2)) == dye;

sol = solve([eqn1, eqn2, eqn3,eqn4], [c1, c2, c3, c4]);


c1sol = sol.c1;
c2sol = sol.c2;
c3sol = sol.c3;
c4sol = sol.c4;


c1sol = vpa(sol.c1,8);
c2sol = vpa(sol.c2,8);
c3sol = vpa(sol.c3,8);
c4sol = vpa(sol.c4,8);

%m2in = 39.3701; %1 m = 39 in
fprintf("Rao Parabola Fit Coefficients:\n");
fprintf("c1 = %g\nc2 = %g\nc3 = %g\nc4 = %g\n\n",c1sol,c2sol,c3sol,c4sol);

fprintf("Rao Parabola Contour (x in m):\n");
fprintf("f(x) = %g*x + %g + (%g*x + %g)^(1/2)\n",c1sol,c2sol,c3sol,c4sol);
fprintf("Start: x = %g m\n",domain3(1));
fprintf("End:   x = %g m\n",domain3(end));

%% Calculate Solution

for i=1:length(domain1)
    x = domain1(i);
   local =  [x,f1(x)];
   field1 = [field1;local];
end

for i=1:length(domain2)
    x = domain2(i);
   local =  [x,f2(x)];
   field2 = [field2;local];
end

for i=1:length(domain3)
    x = domain3(i);
   local =  [x,f3(x)];
   field3 = [field3;local];
end


%% Plot

hold on
plot(field1(:,1), field1(:,2),'r','LineWidth',2)
plot(field2(:,1), field2(:,2),'r','LineWidth',2)
plot(field3(:,1), field3(:,2),'r','LineWidth',2)

plot(field1(:,1), -field1(:,2),'r','LineWidth',2)
plot(field2(:,1), -field2(:,2),'r','LineWidth',2)
plot(field3(:,1), -field3(:,2),'r','LineWidth',2)

xlim([-2 5])
ylim([-3 3])
grid on
daspect([1 1 1])
ylabel('(m)');
xlabel('(m)');
title('F1 Engine with Raos Method');

%% Region Functions

    function y = f1(x)
        y = 2.5*rt - 1.5*rt*cos(asin((-x)/(1.5*rt)));
    end

    function y = f2(x)
        y = 1.4*rt - 0.4*rt*cos(asin(x/(0.4*rt)));
    end

    function y = f3(x)
        y = c1sol*x + c2sol + (c3sol*x + c4sol)^(1/2);
    end


end