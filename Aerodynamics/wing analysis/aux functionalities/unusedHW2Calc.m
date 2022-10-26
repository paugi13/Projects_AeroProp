%% NOT USED ORIGINAL LLWING CODE

%momento de cabeceo 
pi=3.14;
k=0.75;
bf=20;
alpha=[-10*pi/180,20*pi/180];
CMf=258.2996*alpha/(3.02*0.5*(100/3.6)^2*1.225);

figure
hold on
plot(alpha,CMf);
title 'Coeficiente de momento de cabeceo'
xlabel '\alpha[rad]'
ylabel 'M_F [Nm]'
hold off

%constantes
a=0.3;
b=0.6;
lh=2;
SB=pi*(b/2)^2;
% Sb=0.25*SB;
Ss=3.06;
Re=3.91e6;
Cf=0.074/(Re^(1/5));
q=0.5*1.225*(100/3.6)^2;
Vf=3.02;

%Area de b en función de al angulo de ataque.
%S_balpha=pi.*0.3*sqrt((1./(1+(tand(ALPHA).^2)/0.09)).*(1+tand(ALPHA).^2));
%figure
%plot(ALPHA,S_balpha);

S_bal = pi*0.3.*sqrt((1+(tand(ALPHA)).^2).*(1./(1/0.09+(tand(ALPHA)).^2)));
figure
hold on
plot(ALPHA,S_bal);
hold off
%drag parasito

%dfb=Cf*(1+60/((lh/b)^3)+0.0025*(lh/b))*(Ss/SB);
%db=sqrt(Sb/0.7854);
%Cdb=0.029*(db/b)^3/sqrt(Cdfb);

% Fuselage drag
% Cd0=Cf*(1+60/((lh/b)^3)+0.0025*(lh/b))*(Ss/SB);


%% Part 2. 
%1. OK
%drag parasito en función de alpha
Cd01=Cf*(1+60/((lh/b)^3)+0.0025*(lh/b))*(Ss./S_bal);

%momento de cabeceo en función de alfa
Mf=2*k*q*Vf*ALPHA*pi/180;
C_mf=Mf/(q*Vf);             %Munk Method -> Da momentos demasiado altos
C_mf2 = (pi*k*q*2*0.09*ALPHA*pi/180*4/3)/(q*Vf);    %Por integración del pdf 3.3. 
figure
plot(ALPHA,C_mf)
title('Momento de cabeceo')
xlabel('\alpha [º]')
ylabel('Cm_f')
grid on

cm_0 = 0.1197;  %Por interpolación
x_ca_c = -1.498;    %Pendiente de la recta
X_ac = -1*x_ca_c*mac*bf;    %Cálculo a partir de la mean aerodynamic chord.

%5. Cdo - k model. (CON FUSELAJE)
v = linspace(-1.5, 1, 150);
pol = polyval([0.021021 -0.0064986 0.1103], v);
CBW= Cd01+force_coeff(11,:);
figure
plot(force_coeff(7,:),CBW)
hold on
plot(v, pol, 'r');
grid on
xlabel('C_L')
ylabel('CD_{BW}')
hold off

%6. CMcg-CL. OK 

Xcg=1.3;
%Distancias adimensionalizadas.
Cm_cg=cm_0-force_coeff(7,:)*(X_ac - Xcg)/(mac*bf)+C_mf2;  
figure
plot(force_coeff(7,:),Cm_cg);
title 'Coeficiente de momento respecto el centro de gravedad'
xlabel 'C_L'
ylabel 'CM_{cg}'
axis equal
grid on 

