clc;clear;close all;
data = readmatrix('stress_deformation.csv');
stress = data(2:end,8)/320e6;
plastic_strain = data(2:end,5);
F = data(2:end,2);
d = data(2:end,4);
psie = data(2:end,6);
psip = data(2:end,7);
coal = data(2:end,3);
plot(F,stress);
title("stress strain")
figure
plot(F,d);
title("Damage")
figure
plot(psie);
title("psie");
figure
plot(psip);
title("psip");
figure
plot(coal);
title("Coal")
figure
plot(F+plastic_strain,plastic_strain)
title("Strain vs plastic_strain")