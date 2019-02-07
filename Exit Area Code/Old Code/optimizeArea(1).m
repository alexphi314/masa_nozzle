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
burn_t = 40; %burn duration (seconds)
cd = 0.2;
m_dot_ox = 3.61; %Ox flow rate (kg/s)
m_dot_f = 2.89; %Fuel flow rate (kg/s)
g = 9.81; % m/s
qr_fuel = 26.7e6; %J/kg_fuel
fp_rat = 4/9;
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
Qr = qr_fuel * fp_rat;

%% Initialize Trackers %%
% Assumption: rocket travels vertically, straight up
elev = 0;
vel = 0;
dt = 0.1;
t = 0;
Ts = [];
rs = [];
 
while t <= burn_t
    %% Calculate Forces %%
    
    %Standard Atmosphere Model
    [pe, temp, rho] = EvalStdAtm(elev);
    
    % Update mass
    m_p = m_pi - m_dot.*t; 
    m = dry_m + m_p; %total mass: dry mass plus propellant

    D = 0.5.*cd.*rho.*(vel.^2).*SA;
    g_coef = (gam-1)/gam;
    Ue = sqrt(2*Qr*(1-(pe./cham_p).^g_coef));
    T = m_dot.*Ue;
    Ts = [Ts T];
    W = m.*g;

    F = T - D - W;
    
    %% Calculate nozzle size
    a = sqrt(gam*R*temp);
    Me = Ue./a;
    A_rat = ((2./(gam+1)).*(1 + (gam-1).*Me.^2/2)).^((gam+1)./2./(gam-1))./Me;
    Ae = A_rat .* A_star;
    re = sqrt(Ae./pi);
    rs = [rs re];
    
    %% Update Acceleration, Velocity, Elev
    a = F./m;

    vel = vel + a.*dt;

    elev = elev + vel.*dt;

    t = t + dt;
end
    
    
%% Plot
plot(rs, Ts);
 
title('Maximum Thrust as a Function of Exit Radius');
xlabel('Nozzle Exit Radius (m)');
ylabel('Thrust (N)'); % measured in N - make sure to have units right!
 
print('exitArea.png','-dpng');
