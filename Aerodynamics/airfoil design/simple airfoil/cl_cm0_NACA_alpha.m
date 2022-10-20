clc;
clear; close all;


text = "Perfil NACA a analizar:  ";
naca = input(text);

c = 1;
f = fix(naca/1000)/100;                     %Max camber.
p = fix(rem(naca, 1000)/100)/10;            %Position max camber. 
t = rem(fix(rem(naca, 1000)), 100)/100;     %Max thickness.
incr_alpha = 0.01;
n_alpha = -5:incr_alpha:10;
alpha =-5;                                      %Angle of attack.
u_inf = 1;                                      %Freestream.
dens = 1;                                       %Density.

leng = (10-alpha)/incr_alpha;
cl_dist = zeros(1, leng);
cm_dist = zeros(1, leng);
cm_0 = zeros(1, leng);

%% Flap parameters and definition of number of panels 
pan = 200;
point = pan+1;
i=1;

%% Airfoil analysis
while alpha<=10
    point=pan+1;
    [cl_dist(1, i), cm_dist(1, i)] = cl_cm0_NACA_alpha_function(point, pan, f, p, t, alpha, 0, 0);
    cm_0(1, i) = cl_dist(1, i)*0.25 + cm_dist(1, i);
    if int32(alpha) == 5
       m = (cl_dist(1,i)-cl_dist(1,(i-1)))/(alpha-(alpha-incr_alpha));
       al_0 = -cl_dist(1,(i-1))/m + (alpha-incr_alpha);   
    end
    alpha = alpha+incr_alpha;
    i=i+1;
end


%% Postprocess
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure
hold on
title("\textbf{C_l vs \alpha");
plot(n_alpha, cl_dist(1,:), 'b', 'LineWidth', 1);
xlabel('$\alpha$ $\left[\mathrm{^\circ}\right]$')
ylabel('$C_l$ $\left[\mathrm{-}\right]$')
title('Cl  vs \alpha')
grid on
hold off

figure
hold on
title("\textbf{C_l vs \alpha");
plot(n_alpha, cm_0(1,:), 'b');
plot(-5:1:10, zeros(16,1), 'black');
xlabel('$\alpha$ $\left[\mathrm{^\circ}\right]$')
ylabel('$C_{m_0}$ $\left[\mathrm{-}\right]$')
grid on
hold off