clc
clear

%fix() integer part of division.
%rem() r
text = "Perfil NACA a analizar:  ";
naca = input(text);

c = 1;
f = fix(naca/1000)/100;                     %Max camber.
p = fix(rem(naca, 1000)/100)/10;            %Position max camber. 
t = rem(fix(rem(naca, 1000)), 100)/100;     %Max thickness.
alpha = 4;                                      %Angle of attack.
u_inf = 1;                                      %Freestream.
dens = 1;                                       %Density.

% pan = 150;        %N panels require N + 1 points
% point = pan+1;

lim_pan = 200;
pan=1;
n_pan = 1:lim_pan;
cl_dist = zeros(1, lim_pan);
cm_dist = zeros(1, lim_pan);

cl_real = 0.66673;
cm_le_real = -0.21981;

while pan<=200
    point=pan+1;
    [cl_dist(1, pan), cm_dist(1, pan)] = naca_calculus_autom(point, pan, f, p, t, alpha, 0, 0);
    pan = pan+1;
end

error_cl = (cl_real-cl_dist)./cl_dist*100;
error_cm = (cm_le_real-cm_dist)./cm_dist*100;

%Cl plot
figure
colororder({'black','black'});
yyaxis left
plot(n_pan, cl_dist, 'b');
ylabel('Cl', 'Fontsize', 14);
xlabel('Número paneles', 'Fontsize', 14);

hold on
plot(1:1:lim_pan, ones(1, 200)*cl_real, 'black');
yyaxis right

plot(n_pan, error_cl(1,:), 'r');
ylabel('% Error', 'Fontsize', 14);
legend('Cl', 'Cl = ' + string(cl_real), '% Error');
title('Cálculo de Cl DVM vs TAT');
hold off

%Cm plot
figure
colororder({'black','black'});
yyaxis left
plot(n_pan, cm_dist, 'b');
ylabel('Cm_{LE}', 'Fontsize', 14);
xlabel('Número paneles', 'Fontsize', 14);

hold on
plot(1:1:lim_pan, ones(1, 200)*cm_le_real, 'black');
yyaxis right

plot(n_pan, error_cm(1,:), 'r');
ylabel('% Error', 'Fontsize', 14);
legend('Cm_{LE}', 'Cm = ' + string(cm_le_real), '% Error');
title('Cálculo de Cm en borde ataque DVM vs TAT');
hold off

