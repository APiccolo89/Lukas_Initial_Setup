%%%%%% Simple sigmoid function parametrisation
close all;
clf;
clear all; 
C1 =-100;
C2 = 100;



x = -1000:0.1:1000;
gr = 0.5;
v  = 1.0; 

%Boundary = C1+(C2-C1)./(1+exp(-gr*(x-0))).^(1/v);
%Boundary2 = C1+(C2-C1)./(1+exp(-gr*(x-0))).^(1/v);
Boundary3 = C1+(C2-C1)./(1+exp(-0.018*(x-0))).^(1/0.01);


figure(1)
%plot(x,Boundary,'Color','r');
hold on
%plot(x,Boundary2,"Color",'k');
plot(x,Boundary3,'Color','b','LineStyle','-',LineWidth=1.3);
xline(-350)
xline(350)
yline(-800)
yline(-300,LineWidth=1.4)
xline(0)





