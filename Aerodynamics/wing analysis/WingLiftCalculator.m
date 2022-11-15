%% Wing lift calculus
% Main goal: to calculate stall speed

opt = input('Add flap to wing? Y/N (1/0): ');
wantedAoA = input('Wing incidence angle?: ');

numSt = buildStringAD(wantedAoA);
if opt == 0
    direct = join(['wing analysis/workspaces/wingLiftdist', numSt]);
    load(direct);
else
    direct = join(['wing analysis/workspaces/wingLiftdistFlap', numSt]);
    load(direct);
end

rho = 1.225;
MTOW = 39850*9.81;
S = 65.258;

syms x
eq = MTOW == 0.5*rho*S*wingCL*x^2;

sol = double(solve(eq, x));
stallSpeed = sol(2)*3.6;


