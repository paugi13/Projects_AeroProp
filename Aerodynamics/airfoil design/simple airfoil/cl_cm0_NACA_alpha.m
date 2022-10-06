clc

%fix() integer part of division.
%rem() r
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

% pan = 150;        %N panels require N + 1 points
% 
leng = (10-alpha)/incr_alpha;
cl_dist = zeros(1, leng);
cm_dist = zeros(1, leng);
cm_0 = zeros(1, leng);

pan = 200;
point = pan+1;
i=1;

while alpha<=10
    point=pan+1;
    [cl_dist(1, i), cm_dist(1, i)] = cl_cm0_2408_alpha_function(point, pan, f, p, t, alpha, 0, 0);
    cm_0(1, i) = cl_dist(1, i)*0.25 + cm_dist(1, i);
    if int32(alpha) == 5
       m = (cl_dist(1,i)-cl_dist(1,(i-1)))/(alpha-(alpha-incr_alpha));
       al_0 = -cl_dist(1,(i-1))/m + (alpha-incr_alpha);   
    end
    alpha = alpha+incr_alpha;
    i=i+1;
end

figure
plot(n_alpha, cl_dist(1,:), 'b');
hold on
xlabel('\alpha [º]', 'Fontsize', 14)
ylabel('Cl', 'Fontsize', 14);
plot(-5:1:10, zeros(16,1), 'black');
title('Cl  vs \alpha')
hold off

figure

plot(n_alpha, cm_0(1,:), 'b');
hold on
plot(-5:1:10, zeros(16,1), 'black');
xlabel('\alpha [º]', 'Fontsize', 14)
ylabel('C_{m0}', 'Fontsize', 14);
title('C_{m0}  vs \alpha');
hold off