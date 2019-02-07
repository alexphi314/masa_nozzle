close all;
clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Nozzle Exit Area Optimization Code
%
% Project: Big Rocket
%
% Purpose: to determine maximum impulse as a function of exit radius
%
% Changelog
% 2018-09-06 Katie Lerond and Alex Philpott: Initial Creation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize Input Values %%
m_pi = 593 ; %mass of initial propellant (in lbs)
dry_m= 442; %dry mass (in lbs)
SA = pi.* ((13.6).^2)./4; %projected surface area of rocket (in^2)
agl = 1188.72; %elevation of las cruces (in m)
cham_p = 450; %chamber pressure (PSI)
cham_T = 3300; %chamber temperature (K)
burn_t = 40; %burn duration (seconds)
m_dot_ox = 3.61; %Ox flow rate (kg/s)
m_dot_f = 2.89; %Fuel flow rate (kg/s)
g = 9.81; % m/s
gam = 1.2; %% from RPA
air_gam = 1.4;
air_R = 287;
R = 404; %J/(kg-K) From RPA
throat_diam = 66e-3; %m
throat_r = throat_diam/2;
A_star = pi*throat_diam.^2./4;

%% Convert to SI %%
m_pi = m_pi * 0.453592; %kg
dry_m = dry_m * 0.453592; %kg
SA = SA * 0.00064516; %m^2
cham_p = cham_p * 6894.76; %Pa

%% Read in CSV to correlate velocity with Cd
load('rev4e.mat');
elem_num = find(Rev4e1.Timesec > 40);
mach_ref = Rev4e1.Mach(1:elem_num);
cd_ref = Rev4e1.CD(1:elem_num);

%Get only unique entries in mach_ref
[mach_ref_unique, mi] = unique(mach_ref);
cd_ref_unique = cd_ref(mi);

%% Pre-loop Calcs %%
m_dot = m_dot_ox + m_dot_f; %total mass flow rate
n = 40; %Number of cases to run
re_span = linspace(3.25,4.25,n).*0.0254; %m values to span in exit radius
Js = zeros(1,n); %Impulse result array
res = zeros(1,n); %Validtion
Ues = zeros(1,n);
Tes = zeros(1,n);
Tgs = zeros(1,n);
Pes = zeros(1,n);
Pas = zeros(1,n);
M_rs = [];
cds = [];

%Copy mach ref and cd ref to parallel pools
mrC = parallel.pool.Constant(mach_ref_unique);
crC = parallel.pool.Constant(cd_ref_unique);
m_check_value = mach_ref_unique(100);
cd_check_value = cd_ref_unique(100);
fprintf("Now running Alex Philpott's exit area optimizer version 1.0\n");
parfor k = 1:length(re_span)
    fprintf('Now running case %i\n',k);
    re = re_span(k);
    Ae = pi*re.^2;
    A_rat = Ae./A_star;
    
    %% Initialize Trackers %%
    % Assumption: rocket travels vertically, straight up
    elev = agl;
    vel = 0;
    dt = 1;
    t = 0;
    Ts = [];
    
    %% Check that the constant array passed in has correct values
%     if (m_check_value ~= mrC.Value(100))
%         fprintf('Failed mach array check on loop %i! %d != %d',k,m_check_value,mrC.Value(100));
%     elseif (cd_check_value ~= crC.Value(100))
%         fprintf('Failed cd array check %i! %d != %d\n',k,cd_check_value,crC.Value(100));
%     else
%         fprintf('Check passed on loop %i\n',k);
%     end

    while t <= burn_t
        %% Calculate Forces %%

        %Standard Atmosphere Model
        [pa, ta, rho] = EvalStdAtm(elev);

        % Update mass
        m_p = m_pi - m_dot.*t; 
        m = dry_m + m_p; %total mass: dry mass plus propellant
        
        aa = sqrt(air_R*air_gam*ta);
        M_r = vel./aa;
        cd = interp1(mrC.Value(:,:), crC.Value(:,:), M_r);
        if M_r > max(mach_ref_unique)
            cd = cd_ref_unique(end);
        elseif M_r < min(mach_ref_unique)
            cd = cd_ref_unique(1);
        end
        
        %Store for plot
        if k == 1
            cds = [cds cd];
            M_rs = [M_rs M_r];
        end
        D = 0.5.*cd.*rho.*(vel.^2).*SA;
        W = m.*g;
        
        %% Actual Thrust Calculation
        Me = solve_Me(A_rat, gam);
        
        To = cham_T;
        Te = To./(1 + (gam-1)./2.*Me.^2);
        pe = cham_p.*(Te./To).^(gam./(gam-1));
        %fprintf('Exit pressure: %d Pa\n',pe);
        a = sqrt(gam*R*Te);
        Ue = Me.*a;
        Ueq = Ue + (pe - pa).*Ae./m_dot;
        T = m_dot.*Ueq;
        %fprintf('Thrust: %.3d lbf\n',T.*0.224809);
        Ts = [Ts T];
        
        if t == 0
            Ues(1,k) = Ue;
            Tes(1,k) = Te;
            Tgs(1,k) = T;
            Pes(1,k) = pe;
            Pas(1,k) = pa;
        end
        
        F = T - D - W;

        %% Update Acceleration, Velocity, Elev
        a = F./m;
        
        if a < 0
            break;
        end

        vel = vel + a.*dt;

        elev = elev + vel.*dt;

        t = t + dt;
    end

    %% Total Impulse Calculation
    J = trapz(Ts)*0.224809*dt; %lbf-s
    res(1,k) = re.*39.3701; %in
    Js(1,k) = J;
    
    %fprintf('Final elev for case %i: %f m\n',k,elev);
end  

%% Plot
hold on;
plot(res.*2, Js, 'b', 'DisplayName','Impulse (lbf-s)');
title('Total Impulse vs. Exit Diameter');
xlabel('Nozzle Exit Diameter (in)');
ylabel('Impulse (lb-s)'); % measured in Ns - make sure to have units right!
legend('location','Northwest');
print('exitAreaParallel.png','-dpng');

figure; hold on;
plot(res.*2, Tes, 'r', 'DisplayName', 'Te (K)');
plot(res.*2, Ues, 'k', 'DisplayName', 'Ue (m/s)');
plot(res.*2, Tgs, 'g', 'DisplayName', 'T (N)');
plot(res.*2, Pes, 'c', 'DisplayName', 'pe (Pa)');
xlabel('Nozzle Exit Diameter (in)');
legend('location','Northwest');
print('miscValues.png','-dpng');

figure; hold on;
plot(res.*2, Pes, 'c', 'DisplayName', 'Nozzle Exit Pressure (Pa)');
plot(res.*2, Pas, 'm', 'DisplayName', 'Ambient Pressure (Pa)');
plot(res.*2, Tgs, 'g', 'DisplayName', 'Thrust (N)');
xlabel('Nozzle Exit Diameter (in)');
ylabel('Selected Values at 1000 m above MSL (elevation of Las Cruces)');
legend('location','Northwest');
print('PressureThrust.png','-dpng');

figure; hold on;
plot(mach_ref, cd_ref, 'c', 'DisplayName', 'Reference Cd vs. Mach Number Curve');
plot(M_rs, cds, 'm', 'DisplayName', 'Calculated Cd vs. Mach Number Curve');
xlabel('Mach Number');
ylabel('Cd');
legend('location','Southwest');
print('cd_mach.png','-dpng');

%% Test Setup
[mval mindx] = max(Js);
fout = fopen('output.txt','wt');
fprintf(fout,'Max impulse (%d lbf-s) occurs at exit area %d (exit diameter %3.3f in) in2\n',mval,res(mindx)^2*pi, res(mindx).*2);
fprintf(fout,'Max impulse expansion ratio is %3.3f\n',res(mindx)^2/(throat_r*39.3701)^2);
max_re = 4.724*0.0254/2; %m %AA test nozzle exit radius
%max_re = res(mindx)*0.0254;
max_ae = pi*max_re^2;
max_A_rat = max_ae./A_star;
[aa_pa, ta, rho] = EvalStdAtm(0); %MSL elev

max_Me = solve_Me(max_A_rat, gam);

aa_T = calc_T(max_Me,aa_pa,gam,cham_T,cham_p,max_ae,m_dot,R);
fprintf(fout,'Thrust at MSL elevation with AA test config is %d lbf\n',aa_T*0.224809);

%optim_alt = bisection(@EvalStdAtm,max_pe,0,100e3);
%fprintf(fout,'Test Nozzle optimized for %3.3f ft\n',optim_alt.*3.28084);

figure;
plot(res.*2, (Js - mval)/mval*100, 'b', 'DisplayName','Impulse (lbf-s)');
title('Percent Drop from Max Impulse');
xlabel('Nozzle Exit Diameter (in)');
ylabel('Percent Change from Max Impulse (%)'); % measured in Ns - make sure to have units right!
legend('location','Northwest');
print('exitAreaPercent.png','-dpng');

fclose(fout);
fprintf('View output.txt for nozzle sizing parameters.\n');

%% Functions
function T = calc_T(Me,pa,gam,To,po,ae,m_dot,R)
% T = calc_T(Me,pa,gam,To,po,ae,m_dot,R)
Te = To./(1 + (gam-1)./2.*Me.^2);
pe = po.*(Te./To).^(gam./(gam-1));
a = sqrt(gam*R*Te);
Ue = Me.*a;
Ueq = Ue + (pe - pa).*ae./m_dot;
T = m_dot.*Ueq;
end

function [p, T, rho] = EvalStdAtm(h)
% Alex Philpott
% Given altitude, return pressure, temperature, and density
% h: altitude (m)
% 
% p: pressure (Pa)
% T: temperature (K)
% rho: density (kg/m3)

hl = [0, 11000, 20000, 32000, 47000, 51000, 71000, 85000]; % altitudes at layer endpoints [m] 
Tl = [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.65]; % temperature at start of each layer [K] 
Ll = [-.0065, 0, 0.0010, 0.0028, 0, -0.0028, -0.002, 0]; % lapse rate for each layer [K/m]
pl = [101325,22632.0554587517,5474.88465973091,868.017647755643,110.906115784329,66.9387500974116,3.95641034818857,0.363411441797265]; %pressure Pa

g = 9.80665; % gravitational acceleration [m/s?2]  
R = 287.053; % gas constant for air [J/kg.K] 

for l=1:7 % loop over layers 
    lower_h = hl(l);
    upper_h = hl(l+1);
    if h > lower_h && h <= upper_h
        i = l;
        break;
    end
end

if h > hl(end)
    i = 8;
elseif h <= 0
    i = 1;
end

if Ll(i) == 0
    T = Tl(i);
    p = pl(i)*exp(-g*(h-hl(i))./(R*T));
else
    T = Tl(i) + Ll(i)*(h - hl(i));
    p = pl(i)*(T/Tl(i)).^(-g/(R*Ll(i)));
end

rho = p./(R*T);
end

function Me = solve_Me(A_rat, gam)
%Alex Philpott
%Given expansion ratio and gamma, calculate exit mach number
syms Me;
eqn1 = A_rat == (((2./(gam+1)).*(1 + (gam-1).*Me.^2/2)).^((gam+1)./2./(gam-1)))./Me;
soln = solve(eqn1, Me);
Me = 0;
for i = 1:length(soln)
    me = double(soln(i,1));
    if isreal(me) && me > 1
        Me = me;
        break;
    end
end
end

function h = bisection(func, pe, low_h, max_h)
%Given function handle, value, and initial height, implement binary search
%to find that value

h_guess = (max_h-low_h)/2;
pe_out = func(h_guess);

while abs(pe - pe_out) > 1e-7
    if pe_out < pe
        max_h = max_h - (max_h-low_h)/2;
    else
        low_h = low_h + (max_h-low_h)/2;
    end
    h_guess = low_h+(max_h-low_h)/2;
    pe_out = func(h_guess);
end
h = h_guess;
end
