clc
clear

%fix() integer part of division.
%rem() r
text = "Perfil NACA a analizar:  ";
naca = input(text);

f = fix(naca/1000)/100;                     %Max camber.
p = fix(rem(naca, 1000)/100)/10;            %Position max camber. 
t = rem(fix(rem(naca, 1000)), 100)/100;     %Max thickness.
alpha = 0;                                   %Angle of attack.
u_inf = 1;                                      %Freestream.
dens = 1; 
% xh = 0.8;
def = 0;
lim_def = 20;

incr_def = 0.2;
leng = lim_def/incr_def;

pan = 150;
point= pan+1;
i=1;
E = zeros(31,1);

cl_hinge = zeros(leng, 31);
cm_hinge = zeros(leng, 31);
force_hing = zeros(leng, 31);
moment_hinge = zeros(leng, 31);
def_vec = zeros(leng, 1);
j=0;
h = 1;

%Expected behaviour: linear cl due to flaps deflection angles. Same
%situation as an airfoil with increasing angle of attack. 
while j <= 0.3
    while def <=lim_def
        [cl_hinge(i, h), cm_hinge(i, h), force_hing(i, h), moment_hinge(i, h)] = ... 
            forces_hinge_flap(point, pan, f, p, t, alpha, 1-j, def);
        def_vec(i, 1) = def;
        def = def + incr_def;
        i = i+1;
    end
    def = 0;
    i=1;
    E(h)= j;
    h=h+1;
    j=j+0.01;
end

%% Afegit
vec_pend = zeros(31,1);
vec_cm_pend = zeros(31,1);

for i = 1:30
vec_pend(i)= (cl_hinge(6,i)-cl_hinge(5,i))/(def_vec(6)-def_vec(5));
vec_cm_pend(i)= (cm_hinge(6,i)-cm_hinge(5,i))/(def_vec(6)-def_vec(5));
end

figure
plot(E, vec_pend);
hold on
plot(E, vec_cm_pend);
axis equal




