%% Lift computation

%% Load workspaces
% Only one analysis at the time
load('wing analysis/workspaces/wingLiftdist5');


%% Input data
% cruise regime: MTOW?
MTOW = 39850*9.81;
cruiseSpeed = 857/3.6;
rho = 0.363918;
q = 0.5*rho*cruiseSpeed^2;

liftDist = 2*q*wingChord*cl_local(N/2+1:end, 12).*wingChord(N/2+1:end);

%% Plot lift distribution
fig1 = figure(1);
hold on
title("\textbf{Local $C_l$ vs. Spanwise station }");
plot(spanCoords, cl_local(N/2+1:end, 12).*wingChord(N/2+1:end), 'b', 'LineWidth', 1);
xlabel("$x/b$ $\left[\mathrm{-}\right]$");
ylabel("$C_l$ $\left[\mathrm{-}\right]$");
grid on;
grid minor;
box on;
legend('$\alpha = 5^{\circ}$', 'location', 'south');
hold off