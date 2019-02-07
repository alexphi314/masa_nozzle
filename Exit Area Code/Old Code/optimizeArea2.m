% Purpose: to determine maximum impulse as a function of exit radius

% inputs chamber pressure, burn time, coefficient of drag, propellant mass flow (mdot), respectively


function [ ] = optimizeArea(cp, burn_t, cd, m_dot)

m_pi = 593 ; %mass of initial propellant (in lbs)
dry_m= 442; %dry mass (in lbs)
SA = pi.* ((13.6).^2)./4; %projected surface area of rocket (in^2)
agl = 1188.72; %elevation of las cruces (in m)
 
g = 9.8; % in m/s
 
 
%% MAIN %%
 
A_e = new_exit_area; %*this should go at end* updating the exit area
 
%binary search here%
 
ti = 0;
 
%Ground level case%
%using elevation of las cruces - find air pressure in surrounding area and
%then air pressure needed for at nozzle exit, which gives exit area
 
    while (tf <= burn_t){
    %get thrust,drag, weight to find total force
    %find total impulse (force per unit time)
    tf = ti + 0.1; %length of time between each iteration?
    r = density at time tf;
    v = velocity at time tf; %at the exit of nozzle
    a = altitude at time tf;
    
    
        %%Standard Atmosphere Model%% to determine ambient pressure
        %temp measured in degrees Celsius
        if(a < 11000){
            T = 15.04 - .00649.*a;
            P_a = 101.29.*(((T + 273.1)./ 288.08).^5.256);
	}
        elseif (11000 <= a < 25000){
            T = -56.46;
            P_a = 22.65.* e.^(1.73 - .000157.*a);
	}
        elseif (a >= 25000){
            T = -131.21 + .0029.*a;
            P_a = 2.488.*(((T + 273.1)./216.6).^-11.388);
	}
        end
        
    
    m_p = m_pi - m_dot.*tf; 
%mass of initial propellant - mass of propellant lost (this works if m_dot is constant)

    m = dry_m + m_p; %total mass: dry mass plus propellant
    
    drag = cd.*(0.5).*r.*(v.^2).*SA; %where r = density (not radius)
    thrust = %m_dot.*v + A_e.*(P_e-P_a); 
%last term needed? (exit pressure - ambient pressure) should = 0

    total_force = thrust - drag - m.*g; %if moving directly vertical
    j = total_force.*(tf-ti); %total impulse
    
    a = total_force./m; %acceleration
    
    v = a.* (tf-ti); %velocity
    
    s = v.* (tf-ti); %position
    
    
    
    ti = tf; %changes so that what was the initial time is now the final, and we can go through another iteration
    } %end of while loop here
    
    %ending conditions??%
    % use exit pressure to get exit radius <- this is a work in progress
    
    
%lastly make graph that shows max impulse as a function of exit radius
 
plot(exit_rad, j);
 
title('Maximum Impulse as a Function of Exit Radius');
xlabel('Exit Radius (m)');
ylabel('Impulse (J)'); % measured in kg*m/s - make sure to have units right!
 
print('exitArea.png','-dpng');
 
 
    end
