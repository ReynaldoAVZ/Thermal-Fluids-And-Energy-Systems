%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reynaldo Villarreal Zambrano
% January 22, 2025
% TFES Lab 1 - Cooling Tower Lab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc, close all

%% Declare constants

% declare system constants
T_amb = 20.3; % ambient temperature in lab (oC)
P_atm = 868 / 10; % barometric pressure in lab (measured in mbar -> kPa)
Q_dot_in = 1.6; % input power to water heaters (kW)
D_makeup = 7; % inside diameter of makeup water tank (measured in cm -> m)
z = [0, 24.8, 48.3, 71.8, 100] / 100; % height of segments of cooling tower (cm -> m)

% parse data into separate physical quantities
m_win = [20, 30, 40]; % inlet water flow speed (kg/s)
T1 = [21.8, 21.8, 21.9]; % T1, air inlet temperature, dry bulb (oC)
T2 = [10.8, 10.8, 10.9]; % T2, air inlet temperature, wet bulb (oC)
T3 = [23.1, 24.3, 23]; % T3, air outlet temperature, dry bulb (oC)
T4 = [21.4, 22.3, 21.5]; % T4, air outlet temperature, wet bulb (oC)
T5 = [35.9, 32.8, 28.2]; % T5, water inlet temperature (oC)
T6 = [15.8, 17.1, 17.7]; % T6, water outlet temperature (oC)
t1 = [20.4, 21.9, 21.7]; % t1, air temperature at H, wet bulb (oC)
t2 = [19.3, 21.6, 21.4]; % t2, air temperature at H, dry bulb (oC)
t3 = [25.4, 26.7, 24.9]; % t3, water temperature at H (oC)
t4 = [17.0, 18.7, 19.2]; % t4, air temperature at G, wet bulb (oC)
t5 = [18.1, 19.4, 19.9]; % t5, air temperature at G, dry bulb (oC)
t6 = [20.3, 22.7, 22.1]; % t6, water temperature at G (oC)
t7 = [13.9, 15.4, 16.1]; % t7, air temperature at F, wet bulb (oC)
t8 = [15.6, 17, 17.7]; % t8, air temperature at F, dry bulb (oC)
t9 = [16.1, 18.6, 19]; % t9, water temperature at F(oC)
% pressure drop at air inlet (mm H20)
P_b = [10, 10, 10]; % pressure drop at air outlet (mm H20)

% declare collected/given data
T_water_30 = [T6(2), t9(2), t6(2), t3(2), T5(2)];
T_wetbulb_30 = [T2(2), t7(2), t4(2), t1(2), T4(2)];
T_drybulb_30 = [T1(2), t8(2), t5(2), t2(2), T3(2)];

T_water_20 = [T6(1), t9(1), t6(1), t3(1), T5(1)];
T_wetbulb_20 = [T2(1), t7(1), t4(1), t1(1), T4(1)];
T_drybulb_20 = [T1(1), t8(1), t5(1), t2(1), T3(1)];

T_water_40 = [T6(3), t9(3), t6(3), t3(3), T5(3)];
T_wetbulb_40 = [T2(3), t7(3), t4(3), t1(3), T4(3)];
T_drybulb_40 = [T1(3), t8(3), t5(3), t2(3), T3(3)];

T_water_all = [T_water_20; T_water_30; T_water_40];
T_wetbulb_all = [T_wetbulb_20; T_wetbulb_30; T_wetbulb_40];
T_drybulb_all = [T_drybulb_20; T_drybulb_30; T_drybulb_40];

% experiment info
L_initial = [6; 6.25; 6.2];
L_final = [4.6; 5; 5];
t_experiment = [260, 142, 129];
%% 1a.
% On a single figure, plot water temperature (Tw) and wet bulb temperature
% (Twb), for the case of m_dot_w_in = 30 g/s, as a function of cooling 
% tower height (z) at the 5 positions: A=0 cm, F=24.8 cm, G=48.3 cm, 
% H=71.8 cm, B=100 cm. Plot height on the x-axis in units of m, and 
% temperature on the y-axis in units of C. Use the following marker styles:
% (red) Tw; (blue) Twb. Do NOT connect the markers with a line. Include a
% legend. On the same gure, draw the Range (R) and Approach (A) using 
% vertical lines with double arrows; and, label the two lines as R and A 
% appropriately. Draw these lines in Matlab or using other software.

% plot water temperature (Tw) with respect to height
figure;
plot(z, T_water_30, 'rs', 'MarkerFaceColor', 'r');
hold on;

% plot bump temperature (Twb) with respect to height
plot(z, T_wetbulb_30, 'bo', 'MarkerFaceColor', 'b');
hold on;

% declare plot formatting info
xlabel('Height (m)', 'Interpreter', 'latex');
ylabel('Temperature (\(^\circ\mathrm{C}\))', 'Interpreter', 'latex');
xlim([min(z) max(z)]);
ylim([0 max(T_water_30)+10]);
title('Temperature (\(^\circ\mathrm{C}\)) vs Height (m)', 'Interpreter', 'latex');

% find range (R) and approach (A)
R = max(T_water_30) - min(T_water_30); % Range of temperatures for water temperature (inlet & outlet temps)
A = min(T_water_30) - min(T_wetbulb_30); % Approach (how close our wetbulb gets to the water temp)

% plot range (R) and approach (A) on same plot w/ text
plot([0 0],[min(T_water_30) max(T_water_30)],'--r','Linewidth',1.4)
plot([0 0], [min(T_water_30) min(T_wetbulb_30)], '--b','Linewidth',1.4)

% add legend
legend('Water Temperature ($T_w$)', 'Wet Bulb Temperature ($T_{wb}$)', 'Range ($R$)', 'Approach ($A$)', 'Interpreter', 'latex');
hold off;

%% 1b. 
% Plot cooling tower effciency in terms of a percentage on the y-axis as a 
% function of water inlet ow rate ( m_dot_in) in units of g/s on the x-axis.

% declare empty arrays to contain values
Range_all = [];
Approach_all = [];
efficiency_all = [];

% calculate range and approach for all 3 mass flow rates
for i = 1:length(m_win)
    Range_all(i) = max(T_water_all(i, :)) - min(T_water_all(i, :)); % Range of temperatures for water temperature (inlet & outlet temps)
    Approach_all(i) = min(T_water_all(i, :)) - min(T_wetbulb_all(i, :)); % Approach of all mass flow rates
end

% calculate 3 efficiencies
efficiency_all = (Range_all ./ (Range_all + Approach_all)) * 100; % efficiencies in percentage form

% plot mass flow rate vs efficiencies
colors = lines(length(m_win)); % Generate a set of unique colors
figure;
for i = 1:length(m_win)
    plot(m_win(i), efficiency_all(i), 'o', 'MarkerFaceColor', colors(i, :), 'MarkerEdgeColor', colors(i, :));
    hold on;
end

% plot formatting
xlabel('$\dot{m}_{in}$ (g/s)', 'Interpreter', 'latex');
ylabel('Efficiency ($\eta$) (\%)', 'Interpreter', 'latex');
xlim([min(m_win) max(m_win)]);
ylim([0 100]);
title('Cooling Tower Efficiency ($\eta$) (\%) vs Mass Flow Rate $\dot{m}_{in}$ (g/s)', 'Interpreter', 'latex');
legend('$\dot{m}_{1in}$ (g/s) efficiency ($\eta$) (\%)','$\dot{m}_{2in}$ (g/s) efficiency ($\eta$) (\%)','$\dot{m}_{3in}$ (g/s) efficiency ($\eta$) (\%)', 'Interpreter', 'latex');
hold off;

%% 1c.
% plot specific humidity (omega) as a function of cooling tower height (z)
% where height is on the x-axis in units of m and humidity on the y-axis
% in units of kg water vapor/kg of dry air

% declare variables to contain data
specific_humidity_all = []; % contains all specific humidities for 3 mass flow rates at different heights
specific_volume_all = [];
enthalpy_all = [];
% calculate the specific humidity
for i = 1:length(m_win) % for every mass flow rate (3 in our case)
    for j = 1:length(T_water_all(1, :)) % for all temperature recordings (5 per experiment)
        [Tdb, w, phi, h, Tdp, v, Twb] = Psychrometrics('tdb', T_drybulb_all(i, j), 'twb', T_wetbulb_all(i, j), 'p', P_atm); % calculate properties using Psychormetric function
        specific_humidity_all(i, j) = w; % place our specific humidity value such that we can plot later
        specific_volume_all(i, j) = v; % place our specific volume value that will be used later (1d)
        enthalpy_all(i, j) = h; % place our enthalpy of dry air that will be used later (1f)
    end
end

% plot our specific humidities vs height
figure;
plot(z, specific_humidity_all(1, :), "Color", 'green', 'Marker', 'diamond', 'LineStyle','none');
hold on;
plot(z, specific_humidity_all(2, :), "Color", 'blue', 'Marker', 'o', 'LineStyle','none');
plot(z, specific_humidity_all(3, :), "Color", 'red', 'Marker', 'square', 'LineStyle','none');
hold off;

% plot formatting
xlabel('Tower Height (m)', 'Interpreter', 'latex');
ylabel('Specific Humidity ($\omega$)', 'Interpreter', 'latex');
xlim([0 1]);
title('Specific Humidity ($\omega$) vs Tower Height (m)', 'Interpreter', 'latex');
legend('$\dot{m}_{1in}$ (g/s) $\approx$  20','$\dot{m}_{2in}$ (g/s) $\approx$  30','$\dot{m}_{3in}$ (g/s) $\approx$  40', 'Interpreter', 'latex', 'Location','northwest');

%% 1d.
% plot the dry bulb air temperature (Tdb) as a function of cooling tower
% height (z) at the 5 positions. Plot height on the x-axis in units of m 
% and dry bulb air temperature on y-axis in units of oC. Use same style 
% and legend as in plot 1c.

% plot our dry bulb temperatures vs height
figure;
plot(z, T_drybulb_all(1, :), "Color", 'green', 'Marker', 'diamond', 'LineStyle','none');
hold on;
plot(z, T_drybulb_all(2, :), "Color", 'blue', 'Marker', 'o', 'LineStyle','none');
plot(z, T_drybulb_all(3, :), "Color", 'red', 'Marker', 'square', 'LineStyle','none');
hold off;

% plot formatting
xlabel('Tower Height (m)', 'Interpreter', 'latex');
ylabel('Dry Bulb Temperature ($T_{db}$) (\(^\circ\mathrm{C}\))', 'Interpreter', 'latex');
xlim([0 1]);
title('Dry Bulb Temperature ($T_{db}$) (\(^\circ\mathrm{C}\)) vs Tower Height (m)', 'Interpreter', 'latex');
legend('$\dot{m}_{1in}$ (g/s) $\approx$  20','$\dot{m}_{2in}$ (g/s) $\approx$  30','$\dot{m}_{3in}$ (g/s) $\approx$  40', 'Interpreter', 'latex', 'Location','northwest');

%% 1e.
% plot the ratio of water outlet mass flow rate to water inlet mass flow
% rate (m_wout / m_win) on the y axis as a function of inlet water
% temperature (T_win) in units of oC on x-axis

% declare constants
C = 0.0137; % system constant (comes from manufacturer of machine)
m_a = []; % array to contain calculated air flow rates

% calculate m_a
for i = 1:length(m_win) % for every mass flow rate (3 in our case)
    m_a(i) = C * sqrt(P_b(i) / (((1 + specific_humidity_all(i, 5)) * specific_volume_all(i, 5))));
end

% calculate m_v_in and m_v_out
m_vin = [];
m_vout = [];
for i = 1:length(m_win)
    m_vin(i) = specific_humidity_all(i, 1) * m_a(i);
    m_vout(i) = specific_humidity_all(i, 5) * m_a(i);
end

% calculate mass flow rate out
m_wout = m_win + m_vin - m_vout;

% calculate ratio
m_ratio = m_wout ./ m_win;

% find makeup water loss
area = pi * (D_makeup/2)^2;
delta_h = (L_initial - L_final) .* 2.54;
volume_water = delta_h .* area; % cm^3
mass_water = volume_water;
rate_of_water_loss = mass_water ./ t_experiment';
avg_water_loss = mean(rate_of_water_loss);

% plot
figure
plot(T_water_all(:, 1), m_ratio, 'bo', 'MarkerFaceColor', 'b'); % blue circles

% add labels
xlabel('Inlet Water Temperature, $T_{w,\mathrm{in}}$ ($^\circ$C)', 'Interpreter', 'latex');
ylabel('Mass Flow Rate Ratio, $\frac{\dot{m}_{w,\mathrm{out}}}{\dot{m}_{w,\mathrm{in}}}$', 'Interpreter', 'latex');

% add title
title('Mass Flow Rate Ratio vs Inlet Water Temperature', 'Interpreter', 'latex');

% grid and formatting
grid on;

%% 1f.
% On a single figure, plot the heat transfer rates ( Qa and Qamb) in units of
% kW on the y-axis as a function of inlet water temperature (Twin) in units
% of C on the x-axis. Use a different marker style for Qa and Qamb. Include
% a legend

% calculate Q_a
Q_a = m_a' .* (enthalpy_all(:, 5) - enthalpy_all(:, 1)) ./ 1000;

% calculate Q_amb
Q_amb = (Q_dot_in + m_a' .* (enthalpy_all(:, 1) - enthalpy_all(:, 5))) ./ 1000;

% plot
figure;
plot(T_water_all(:, 1), Q_amb, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', '$Q_\mathrm{amb}$'); % red circles with line
hold on;
plot(T_water_all(:, 1), Q_a, 'bs', 'MarkerFaceColor', 'b', 'DisplayName', '$Q_\mathrm{a}$'); % blue squares with line
hold off;
% add labels
xlabel('Inlet Water Temperature, $T_\mathrm{win}$ ($^\circ$C)', 'Interpreter', 'latex');
ylabel('Heat Transfer Rate (kW)', 'Interpreter', 'latex');

% add legend
legend('Interpreter', 'latex', 'Location', 'best');

% add title
title('Heat Transfer Rates vs Inlet Water Temperature', 'Interpreter', 'latex');

% grid for better visualization
grid on;