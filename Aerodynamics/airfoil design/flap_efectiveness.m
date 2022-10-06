clc

% -----------------------------------------------------------------
%       cl_cm0_2408_alpha must be run before to obtain al_0
% -----------------------------------------------------------------


%fix() integer part of division.
%rem() r
text = "Perfil NACA a analizar:  ";
naca = input(text);

c = 1;
f = fix(naca/1000)/100;                     %Max camber.
p = fix(rem(naca, 1000)/100)/10;            %Position max camber. 
t = rem(fix(rem(naca, 1000)), 100)/100;     %Max thickness
alpha =[5 6];                                      %Angle of attack.
u_inf = 1;                                      %Freestream.
dens = 1;                                       %Density.
def = 20;
E = [0 0.05 0.1 0.15 0.2 0.25 0.3];
effec_exp = [0 0.125 0.25 0.35 0.45 0.525 0.6];


        %N panels require N + 1 points
pan = 150;
point= pan+1;
i=1;
effectiveness = zeros(1,7);
while i<=7
[cl_1, cm_dist(1, i), pos] = cl_cm0_2408_flap_function(point, pan, f, p, t, alpha(1), 1-E(i), def);
[cl_2, cm_dist(1, i), pos] = cl_cm0_2408_flap_function(point, pan, f, p, t, alpha(2), 1-E(i), def);
m = (cl_2-cl_1)/(alpha(2)-(alpha(1)));
effectiveness(1, i) = (al_0-(-cl_1/m + alpha(1)))/def;  

i=i+1;
end

lim_e = 0.300;
e = 0;
incr_e = 0.05;

i = 1;
effect_theory = zeros(1,int32(lim_e/incr_e));
e_vec = zeros(1,int32(lim_e/incr_e));

while e<=lim_e
[cl_1, cm_dist(1, i), pos] = naca_calculus_autom(point, pan, f, p, t, alpha(1), 1-e, def);
[cl_2, cm_dist(1, i), pos] = cl_cm0_2408_flap_function(point, pan, f, p, t, alpha(2), 1-e, def);
m = (cl_2-cl_1)/(alpha(2)-(alpha(1)));
effect_theory(1, i) = (al_0-(-cl_1/m + alpha(1)))/def;  
e_vec(1,i) = e;

e = e + incr_e;
i=i+1;
end

fact_corr = effec_exp(1,:)./effectiveness(1,:);

figure
colororder({'black','black'});
yyaxis left
plot(e_vec(1,:), effect_theory(1,:), 'b');
ylabel('\Delta\alpha_0 / \Delta\eta');
hold on
plot(E(1,:), effec_exp(1,:), 'r', 'LineStyle', '-');
yyaxis right
scatter(E(1,:), fact_corr(1,:),[], [0.8500 0.3250 0.0980],  'filled');
ylabel('Factor de correción');
title('Eficiencia de flap para \eta = 20º');
xlabel('E');
legend('Eficiencia teórica', 'Datos experimentales', 'Factor de corrección', 'Location','northwest');
hold off


