%% Engine positioning
% Script to calculate potencial engine locations

heightTren = 1.5;
desiredSpanPos = 3;
diedral = 6;
diamNacelle = 1.5;

minHeight = heightTren + desiredSpanPos*tand(diedral) - diamNacelle;
disp(minHeight);



