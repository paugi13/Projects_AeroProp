%% CODE FOR HPT ANALYSIS

v1 = Mi*sqrt(gamma_c*Ti*287)*sind(alpha_1);
Cp = 287*gamma_c/(gamma_c-1);

syms x
eq = Cp*AT_real == x^2 - x*2*v1;
U = solve(eq,x);
U = double(U(2));


%% 3

M1_5 = 0.4*sqrt(Ti/(Ti+4*AT_real));
T1_5 = Ti + 4*AT_real;
omega = U/rm;
Uh = omega*rh;
Ut = omega*rt;
aux = 2*M1_5*sqrt(gamma_c*287*T1_5);
V1_5=M1_5*sqrt(gamma_c*287*T1_5);
k1=V1_5*sind(30)*rm;
v1_h=k1/rh;
v1_t=k1/rt;
C_p = U*(U-aux*sind(30));
alpha_h=asind(v1_h/V1_5);
alpha_t=asind(v1_t/V1_5);

%% 4%% 5
R=0.5;
gamma_t=1.3;
w=omega;
tau_t=0.8036*0.8958;
Tt_4=1777.8;
AT_turb=(Tt_4)*(1-tau_t)/3;
C_pt=gamma_t/(gamma_t-1)*287;
U_t=sqrt(C_pt*AT_turb);
rm_t=U_t/w;

%alpha_2=45;
v2_t=U_t;
X=1:0.1:1.5;



alpha_2=45;
V2_t=v2_t/sind(alpha_2);
T2=Tt_4-1/(2*C_pt)*V2_t^2;
Mt=U_t/sqrt(gamma_t*287*T2);
M2=V2_t/sqrt(gamma_t*287*T2);
%f=1/((gamma_t-1)*Mt)*(1-tau_t)==(M2*sind(alpha_2))/(1+(gamma_t-1)/2*M2^2);
mfp_M2=sqrt(gamma_t)*(M2)/(1+(gamma_t-1)/2*M2^2)^((gamma_t+1)/(2*(gamma_t-1)));
mfp_M1_2=sqrt(-(sind(alpha_2)^2-1).*X.^2)*mfp_M2;
mfp_M1=sqrt((2*gamma_t)/(gamma_t+1))*((gamma_t+1)/(3*gamma_t-1))^((gamma_t+1)/(2*(gamma_t-1))).*X;
syms  M1
f=sqrt(gamma_t)*(M1)/((1+(gamma_t-1)/2*M1^2)^((gamma_t+1)/(2*(gamma_t-1))))-0.45==0;
S=vpasolve(f, M1);