%% Engine positioning
% Script to calculate potencial engine locations

heightTren = 1.65;       % [m]
desiredSpanPos = 2;     % [m]
diedral = 6;            % [m]
diamNacelle = 1.5;      % [m]
reqDist = 0.1778;       % [m]
zWing = 0.05;

minHeight = zWing + heightTren + desiredSpanPos*tand(diedral) - diamNacelle;

if minHeight < reqDist
    error('Nacelle position does not meet the requirements');
else
    disp(['Distance = ', num2str(minHeight), ' m']);
    margin = minHeight - reqDist;
    disp(['Margin = ', num2str(margin), ' m']);
end



