ALPHA = [ -10. -8.0 -4.0 0. 4.0 8.0 10. ];

sb1 = pi*0.3.*sqrt((1+(tand(ALPHA)).^2).*(1./(1/0.09+(tand(ALPHA)).^2)));
figure
plot(ALPHA,sb1);