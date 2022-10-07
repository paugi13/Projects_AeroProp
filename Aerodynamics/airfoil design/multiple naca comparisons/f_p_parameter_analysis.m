clc
clear

%% CODE THAT ENABLES COMPARISONS BETWEEN NACA AIRFOILS GIVEN YY FORM XXYY NOM.

%fix() integer part of division.
%rem() r
%p is not given 

text = "Ultimos dos digitos del perfil NACA:  ";
t = input(text);
c = 1;
f = [0 0.02 0.04 0.06];
p = 0.05;
t = t/100;  
alpha = [-1 2];                                 %Two aleatory angles.
u_inf = 1;                                      %Freestream.
dens = 1;                                       %Density.

% pan = 150;        %N panels require N + 1 points
% point = pan+1;
incr_p = 0.01;
pan = 200;
point=pan+1;
cm_dist = zeros(1, int32((0.6-p)/incr_p));
p_vec = zeros(1, int32((0.6-p)/incr_p));
al_0 = zeros(4, int32((0.6-p)/incr_p));
cm_0 = zeros(4, int32((0.6-p)/incr_p));
xh = 0;
def = 0;

i=1;
j=1;

while j<=4
    while p<=0.6

%It shouldnt matter which cm we pick as cm_0 should not be related with
%alpha
[cl_1, cm_dist(1, i)] = naca_calculus_autom(point, pan, f(j), p, t, alpha(1), xh, def);
[cl_2, aux] = naca_calculus_autom(point, pan, f(j), p, t, alpha(2), xh, def);
%aux ha no real use.
m = (cl_2-cl_1)/(alpha(2)-alpha(1));
al_0(j, i) = -cl_1/m + alpha(1);          %From a line's eq.
cm_0(j, i) = cl_1*0.25 + cm_dist(1, i);   %cm_0 is always the same. It should
                                          %matter the angle of attack from which
p_vec(1, i) = p;                          %it is calculated
p = p+incr_p;
i=i+1;                                          
    end
    j = j+1;
    i = 1;
    p = 0.05;
end


%% PLOTTING
figure
plot(p_vec(1,:), al_0(1,:), 'r');
hold on
plot(p_vec(1,:), al_0(2,:), 'b');
plot(p_vec(1,:), al_0(3,:), 'color', [0.6350 0.0780 0.1840]);
plot(p_vec(1,:), al_0(4,:), 'color', [0.4660 0.6740 0.1880] );
% plot(0:0.1:0.6, zeros(7,1), 'black');
xlabel('p/c',  'Fontsize', 14);
ylabel('\alpha_{l0}', 'Fontsize', 14);
legend('f = 0', 'f = 0.02', 'f = 0.04', 'f = 0.06' , 'location', 'southwest');
title('\alpha_{l0} en función de p y f');
hold off
% 

figure
plot(p_vec(1,:), cm_0(1,:), 'r');
hold on
plot(p_vec(1,:), cm_0(2,:),'b');
plot(p_vec(1,:), cm_0(3,:), 'color', [0.6350 0.0780 0.1840]);
plot(p_vec(1,:), cm_0(4,:), 'color', [0.4660 0.6740 0.1880]);
% plot(0:0.1:0.6, zeros(7,1), 'black');
xlabel('p/c', 'Fontsize', 14);
ylabel('C_{m0}', 'Fontsize', 14);
legend('f = 0', 'f = 0.02', 'f = 0.04', 'f = 0.06', 'location', 'southwest');
title('Momento libre en función de p y f');
hold off