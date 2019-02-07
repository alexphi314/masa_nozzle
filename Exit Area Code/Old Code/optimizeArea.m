close all;
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
cham_T = 600; %chamber temperature (K)
burn_t = 40; %burn duration (seconds)
cd = 0.2;
m_dot_ox = 3.61; %Ox flow rate (kg/s)
m_dot_f = 2.89; %Fuel flow rate (kg/s)
g = 9.81; % m/s
gam = 1.4; %% check
R = 287.053;
throat_diam = 69.29e-3; %m
A_star = pi*throat_diam.^2./4;

%% Convert to SI %%
m_pi = m_pi * 0.453592; %kg
dry_m = dry_m * 0.453592; %kg
SA = SA * 0.00064516; %m^2
cham_p = cham_p * 6894.76; %Pa

%% Pre-loop Calcs %%
m_dot = m_dot_ox + m_dot_f;
n = 15;
re_span = linspace(0.127,0.254,n);
Js = zeros(1,n);
res = zeros(1,n);
for k = 1:length(re_span)
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

    while t <= burn_t
        %% Calculate Forces %%

        %Standard Atmosphere Model
        [pe, temp, rho] = EvalStdAtm(elev);

        % Update mass
        m_p = m_pi - m_dot.*t; 
        m = dry_m + m_p; %total mass: dry mass plus propellant

        D = 0.5.*cd.*rho.*(vel.^2).*SA;
        W = m.*g;
        
        %% Actual Thrust Calculation
        Me = solve_Me(A_rat, gam);
        
        To = cham_T;
        Te = To./(1 + (gam-1)./2.*Me.^2);
        a = sqrt(gam*R*Te);
        Ue = Me.*a;
        T = m_dot.*Ue;
        Ts = [Ts T];
        
        F = T - D - W;

        %% Update Acceleration, Velocity, Elev
        a = F./m;

        vel = vel + a.*dt;

        elev = elev + vel.*dt;

        t = t + dt;
    end

    %% Total Impulse Calculation
    J = trapz(Ts);
    res(1,k) = re;
    Js(1,k) = J;
end  

%% Plot
plot(res, Js, 'b', 'DisplayName','Realistic Impulse');
 
title('Total Impulse as a Function of Exit Radius');
xlabel('Nozzle Exit Radius (m)');
ylabel('Impulse (Ns)'); % measured in Ns - make sure to have units right!
 
print('exitArea.png','-dpng');

%% Functions
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
elseif h == 0
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
