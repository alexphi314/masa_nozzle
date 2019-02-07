function [L, H] = gen_line(Dy, slope, res)
%Generate linear nozzle section
%Dy: Total change in y of section (length)
%slope: The slope of the line (i.e. dh/dy)
%res: The number of points in the section

L = linspace(0,Dy,res);
H = L.*slope;

end