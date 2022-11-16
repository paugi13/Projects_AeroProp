%% Engine positioning
% Script to calculate potencial engine locations

heightTren = 1.4;       % [m]
desiredSpanPos = 3;     % [m]
diedral = 6;            % [m]
diamNacelle = 1.5;      % [m]
reqDist = 0.1778;       % [m]

minHeight = heightTren + desiredSpanPos*tand(diedral) - diamNacelle;

if minHeight < reqDist
    error('Nacelle position does not meet the requirements');
else
    disp(['Distance = ', num2str(minHeight), ' m']);
    margin = minHeight - reqDist;
    disp(['Margin = ', num2str(margin), ' m']);
end



