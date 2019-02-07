function [L, H] = gen_circ(R, theta_i, theta_f, res)
%Generate a circular section of the nozzle
%R: Radius of the circle
%theta_i: Initial angle in degrees
%theta_f: Final angle in degrees
%res: number of points
theta = linspace(theta_i,theta_f,res);
L = R.*cosd(theta);
H = R.*sind(theta);

L = L - L(1);
H = H - H(1);
end