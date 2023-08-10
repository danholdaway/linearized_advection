close all
clear
clc
mydir = pwd;

jac

x = 1:64;


figure
contour(J3rd)
colorbar

figure
contour(J3rdul)
colorbar

figure
contour(J3rdrl)
colorbar

figure
contour(Jppmunl)
colorbar

figure
contour(Jlin)
colorbar

figure
contour(Jcwl)
colorbar

figure
contour(Jcsl)
colorbar

figure
contour(Jslice)
colorbar

figure
contour(Jslicebs)
colorbar








cd(mydir)