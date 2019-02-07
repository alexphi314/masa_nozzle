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