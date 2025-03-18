%Script to practice using latex interpreter commands when plotting to make pretty plots
% Eric Pardyjak 1/30/2025
clear all; close all;
%Make up some data to plot
cp = 1004; %J/kg-K specific heat const pressure air
m_dot = 1; %kg/s mass flow rate
T_hot = [230 250 289 320]; % degrees C
T_cold = [200 220 235 319]; % degrees C

DeltaT = T_hot-T_cold; %compute the temperature difference

%Compute the heat flux 
Q_dot = m_dot*cp*(T_hot-T_cold)/1000; %Heat transfer rate in kWatts or J/s
Q_err = 0.25*Q_dot.*ones(size(DeltaT)); %Calculate fake error

%Make a new figure with lots of cool math fonts in it using latex
%interpreter
figure;
%plot(DeltaT,Q_dot,'ro','MarkerFaceColor','r');
errorbar(DeltaT,Q_dot,Q_err,'ro','MarkerFaceColor','r'); %scatter plot with errorbars
set(gca,'XMinorTick','on','YMinorTick','on'); %add minor ticks on x and y axis
set(gca,'FontSize',16,'FontName','Times'); %set x and y axis font sizes of numbers on ticks
xlabel('Temperature difference $\Delta T$ ($^{\circ}$C)','Interpreter','latex','FontSize',16);
ylabel('$\dot{Q}$ (kW)','Interpreter','latex','FontSize',16);
legend('efficiency, $\eta$','interpreter','latex','Location','northwest','FontSize',16);
legend('boxoff');
title('Heat flux as a function on temperature difference','Interpreter','latex','FontSize',20);
text(20,20,'testing adding text to the plot, $\Gamma$ or $\gamma$','Interpreter','latex','FontSize',16);
text(20,15,'$\dot{m}_{in}$ is mass flow rate in','Interpreter','latex','FontSize',16); 
text(20,5,'$\frac{\partial T}{\partial x}$ partials','Interpreter','latex','FontSize',16); 
