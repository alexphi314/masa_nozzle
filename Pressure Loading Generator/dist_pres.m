function [] = dist_pres(geometry,FOS)
%Generates a tabulated pressure distrobution which is used to do FEA
%(imported to ANSYS) geometry generated by gen_nozzle.m
close all;

geometry = load(geometry);
gamma = 1.1902;
p0 = 3.103e+6; %Pa
T0 = 3300; %K
A_t = min(geometry(:,2));
converge = geometry(1:find(geometry(:,2)==A_t),2);
diverge = geometry(find(geometry(:,2)==A_t)+1:end,2);
for i = [1:size(converge,1)]
    M_vect(i) = fsolve(@(M) A_t/M*((2/(gamma+1))*(1+(gamma-1)/2*M^2))^((gamma+1)/(2*(gamma-1))) - converge(i),.1);
    if round(M_vect(i),4) == 0
        M_vect(i) = fsolve(@(M) A_t/M*((2/(gamma+1))*(1+(gamma-1)/2*M^2))^((gamma+1)/(2*(gamma-1))) - converge(i),1);
    end
end

for i = [1:size(diverge,1)]
    M_vect(i+size(converge,1)) = fsolve(@(M) A_t/M*((2/(gamma+1))*(1+(gamma-1)/2*M^2))^((gamma+1)/(2*(gamma-1))) - diverge(i),3);
    if round(M_vect(i),4) == 0
        M_vect(i) = fsolve(@(M) A_t/M*((2/(gamma+1))*(1+(gamma-1)/2*M^2))^((gamma+1)/(2*(gamma-1))) - diverge(i),1.2);
    end
end

T = T0./(1 + ((gamma - 1)/2).*M_vect.^2); %K
p = p0.*(T./T0).^(gamma/(gamma-1));
T = T - 273.15; %C
figure;
plot(geometry(:,1).',p.*FOS);
figure;
plot(geometry(:,1).',T.*FOS);
p_curve = [geometry(:,1),FOS.*p.'];
t_curve = [geometry(:,1),FOS.*T.'];
save('p_curve.txt','p_curve','-tabs','-ascii');
save('t_curve.txt','t_curve','-tabs','-ascii');
