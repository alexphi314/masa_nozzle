function [L, H] = gen_bell(x1,x2, a, b, c, d, res)
%Generates a parabolic section of the nozzle.
%y1: starting y value
%y2: starting y value
%Dy: The total changed in y for the section (i.e. the length)
%a: The quardatic coefficient in h(y) = ay^2 + by
%b: The linear coefficient in h(y) = ay^2 + by
%res: the number of points in the section

L = linspace(x1,x2,res);
%H = (-b + sqrt(b^2 + 4*a.*(L-c)))/(2*a);
H = a.*L + b + (c.*L + d).^(1/2);
L = L - L(1);
H = H - H(1);
end