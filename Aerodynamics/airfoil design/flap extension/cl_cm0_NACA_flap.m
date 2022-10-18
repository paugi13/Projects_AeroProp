clc;
clear; close all;
format long;

%% CODE FOR FLAP ANALYSIS 

text = "Perfil NACA a analizar:  ";
naca = input(text);

c = 1;
f = fix(naca/1000)/100;                     %Max camber.
p = fix(rem(naca, 1000)/100)/10;            %Position max camber. 
t = rem(fix(rem(naca, 1000)), 100)/100;     %Max thickness.
incr_alpha = 0.01;
n_alpha = -5:incr_alpha:10;
alpha =-5.0;                                      %Angle of attack.
u_inf = 1;                                      %Freestream.
dens = 1; %Density.

leng = (10-alpha)/incr_alpha;
cl_dist = zeros(1, leng);
cm_dist = zeros(1, leng);
cm_0 = zeros(1, leng);


%% Flap parameters and definition of number of panels 
xh = 0.8;   % hinge position out of 1. 
def = 20;   % flap deflection

% Number of panels. Depending on this number the results are more or less
% precise.
        %N panels require N + 1 points
pan = 150;
point= pan+1;
i=1;

%% Airfoil with flap building and analysis. 
while alpha<=10
    [cl_dist(1, i), cm_dist(1, i), pos] = cl_cm0_2408_flap_function(point, pan, f, p, t, alpha, xh, def);
    cm_0(1, i) = cl_dist(1, i)*0.25 + cm_dist(1, i);
    if alpha == 5.0
       m = (cl_dist(1,i)-cl_dist(1,(i-1)))/(alpha-(alpha-incr_alpha));
       al_0 = -cl_dist(1,(i-1))/m + (alpha-incr_alpha);   
    end
    alpha = alpha+incr_alpha;
    i=i+1;
end

i=1;
n = [0,0,1];
vec_tan = zeros(pan, 2);
vec_norm = zeros(pan, 2);
while i<=pan
    
    vec_tan(i, 1) = pos(i+1,1) - pos(i,1);
    vec_tan(i, 2) = pos(i+1,2) - pos(i,2);
    vec_incr = [vec_tan(i, 1), vec_tan(i, 2), 0];
    v_normal = cross(n, vec_incr);
    vec_norm(i, 1) = v_normal(1)/norm(v_normal);
    vec_norm(i, 2) = v_normal(2)/norm(v_normal);
    
     i = i+1;
end

i = 1;
cp = zeros(pan, 2);
vortex_point = zeros(pan, 2);
while i<=pan
    cp(i, 1) = pos(i, 1) + 0.75*vec_tan(i, 1);
    cp(i, 2) = pos(i, 2) + 0.75*vec_tan(i, 2);
    vortex_point(i, 1) = pos(i, 1) + 0.25*vec_tan(i, 1);
    vortex_point(i, 2) = pos(i, 2) + 0.25*vec_tan(i, 2);
    i = i+1;
end

%% Postprocess
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure
quiver(cp(:,1), cp(:,2), vec_norm(:,1), vec_norm(:,2));
hold on
plot(pos(:,1),pos(:,2))
axis equal
xlabel('x/c');
ylabel('z');
title("\textbf{Plot $y$ vs. $x$ -- Example}");
hold off

% alpha - Cl
figure
hold on
title("\textbf{C_l vs \alpha");
plot(n_alpha, cl_dist(1,:), 'b', 'LineWidth', 1);
xlabel('$\alpha$ $\left[\mathrm{^\circ}\right]$')
ylabel('$C_l$ $\left[\mathrm{-}\right]$')
grid on
hold off

%alpha - Cm0
figure
hold on
title("\textbf{C_{m_0} vs \alpha");
plot(n_alpha, cm_0(1,:), 'b', 'LineWidth', 1);
xlabel('$\alpha$ $\left[\mathrm{^\circ}\right]$')
ylabel('$C_{m_0}$ $\left[\mathrm{-}\right]$')
grid on
axis equal
hold off