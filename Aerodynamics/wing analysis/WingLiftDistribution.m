%% Lift computation

%% Load workspace for desired incidence angle
% Only one analysis at the time

wantedAoA = input('Wing incidence angle?: ');
direct = join(['wing analysis/workspaces/wingLiftdist', num2str(wantedAoA)]);
disp(join(['wing analysis/workspaces/wingLiftdist', num2str(wantedAoA)]));
load(direct);

%% Input data
% cruise regime: MTOW?
MTOW = 39850*9.81;
cruiseWeight = MTOW * 0.85;
cruiseSpeed = 857/3.6;
rho = 0.363918;
q = 0.5*rho*cruiseSpeed^2;


%% Calculations
liftDist = 2*q*wingChord*cl_local(N/2+1:end, 12).*wingChord(N/2+1:end);

d = 5;
dist = @(x) 0;
for i = 1:length(polinomialFit)
    dist = dist + polinomialFit(i)*x^d;
    d = d - 1;
end

TotalLift = 2*q*integral(dist, spanCoords(1), spanCoords(end));


%% Plot lift distribution
fig1 = figure(1);
hold on
title("\textbf{Local $C_l$ vs. Spanwise station }");
plot(spanCoords, Cl_Values, 'b', 'LineWidth', 1);
xlabel("$x/b$ $\left[\mathrm{-}\right]$");
ylabel("$C_l$ $\left[\mathrm{-}\right]$");
grid on;
grid minor;
box on;
legend('$\alpha = 5^{\circ}$', 'location', 'south');
hold off