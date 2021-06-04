clc;clear;close all;
data = readmatrix('stress_deformation.csv');
stress = data(:,5)/320e6;
plastic_strain = data(:,4);
F = data(:,2)-1;
d = data(:,3);
plot(F(2:end)+plastic_strain(2:end),stress(2:end));
figure
plot(d);