clc

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
E = [0.15 0.2 0.25 0.3];

cl_hinge = zeros(leng, 4);
cm_hinge = zeros(leng, 4);
force_hing = zeros(leng, 4);
moment_hinge = zeros(leng, 4);
def_vec = zeros(leng, 1);

%Expected behaviour: linear cl due to flaps deflection angles. Same
%situation as an airfoil with increasing angle of attack. 
for j = 1:4
    while def <=lim_def
        [cl_hinge(i, j), cm_hinge(i, j), force_hing(i, j), moment_hinge(i, j)]...
            = forces_hinge_flap(point, pan, f, p, t, alpha, 1-E(j), def);
        def_vec(i, 1) = def;
        def = def + incr_def;
        i = i+1;
    end
    def = 0;
    i=1;
end

%Coefficients
figure
colororder({'black','black'});
plot(def_vec(:,1), cl_hinge(:,1),'b','LineStyle','-' );
hold on
plot(def_vec(:,1), cl_hinge(:,2), 'b', 'LineStyle','--');
plot(def_vec(:,1), cl_hinge(:,3),'b','LineStyle',':' );
plot(def_vec(:,1), cl_hinge(:,4),'b','LineStyle','-.');
xlabel('Ángulo de deflexión [º]', 'Fontsize', 14);
yyaxis left
ylabel('Cl', 'Fontsize', 14);

yyaxis right
plot(def_vec(:,1), cm_hinge(:,1),  'r');
plot(def_vec(:,1), cm_hinge(:,2), 'r');
plot(def_vec(:,1), cm_hinge(:,3), 'r');
plot(def_vec(:,1), cm_hinge(:,4), 'r');
ylabel('Cm', 'Fontsize', 14);
title('Representación de coeficientes');
legend('E = 0.15', 'E = 0.2', 'E = 0.25', 'E = 0.3');
hold off

%Forces and moments
figure
colororder({'black','black'});
plot(def_vec(:,1),  force_hing(:,1), 'b','LineStyle','-' );
hold on
plot(def_vec(:,1), force_hing(:,2), 'b','LineStyle','--' );
plot(def_vec(:,1), force_hing(:,3), 'b','LineStyle',':' );
plot(def_vec(:,1), force_hing(:,4), 'b','LineStyle','-.' );
xlabel('Ángulo de deflexión [º]', 'Fontsize', 14);
yyaxis left
ylabel('Lift [N/m]', 'Fontsize', 14);

yyaxis right
plot(def_vec(:,1), moment_hinge(:,1), 'r');
plot(def_vec(:,1), moment_hinge(:,2), 'r');
plot(def_vec(:,1), moment_hinge(:,3), 'r');
plot(def_vec(:,1), moment_hinge(:,4), 'r');
ylabel('Momento [N]', 'Fontsize', 14);
title('Representación de fuerzas y momentos');
legend('E = 0.15', 'E = 0.2', 'E = 0.25', 'E = 0.3');
hold off