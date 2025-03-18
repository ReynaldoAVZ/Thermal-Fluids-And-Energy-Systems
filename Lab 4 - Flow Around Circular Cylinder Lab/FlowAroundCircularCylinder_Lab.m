%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reynaldo Villarreal Zambrano
% February 19, 2025
% TFES Lab 4 - Flow Around Circular Cylinder Lab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc, close all

%% Declare constants

% Diameter
D = .75; % in
D = D * .0254; % in -> m

% Ideal gas law constants
R = 287.05; % J / Kg * K

T = 20.3; % degrees C
T = T + 273.15; % C -> K

P_atm = 868; % mbar 
P_atm = P_atm * 100; % mbar -> Pa

% calculate air density
rho = P_atm / (R * T); % kg/m^3

P_convert = 133.3224; % conversion from mmHg to Pa

v = 1.83*10^-5;
u = sqrt(2*(.7828)*133.32/rho);
Re_calc = u*D/v;
%% Calculate U_inf
yD_up = [0, 1, 2, 3, 4]; 
Uinfm = zeros(size(yD_up)); 

% Loop through all data files for U_inf
for i = 1:length(yD_up)
    % Construct filename
    FileName = ['Uinf_yD', num2str(yD_up(i))];

    % Read and convert pressure data
    Delta_p = readmatrix(FileName) * P_convert; 

    % Compute U_inf values
    Uinf = sqrt(2 .* Delta_p ./ rho);
    Uinfm(i) = mean(Uinf);
end

% Compute average U_inf
Uinf_avg = mean(Uinfm);

%% Calculate Uwake Data
yD = ["00", "01", "02", "03", "04", "05", "06", "08", "10", "12", ...
      "14", "16", "18", "20", "22", "24", "26", "30", "35", "40"];

% Preallocate arrays for performance
umean = zeros(size(yD));
uStd = zeros(size(yD));
uErr = zeros(size(yD));

for i = 1:length(yD)
    % Construct filename
    FileName = ['Uwake_yD' + yD(i)];

    % Read and convert pressure data
    Delta_p = readmatrix(FileName) * P_convert; 

    % Compute Uwake values
    Uwake = sqrt(2 .* Delta_p ./ rho);
    
    % Store computed values
    umean(i) = mean(Uwake);
    uStd(i) = std(Uwake);
    
    % Compute 95% confidence interval error bars
    uErr(i) = 2 * uStd(i) / sqrt(length(Delta_p)); 
end

% Convert yD to numeric values
yD_numeric = str2double(yD) / 10;

%% Plot Figure 1a: Normalized Mean Horizontal Velocity
figure;
errorbar(yD_numeric, umean / Uinf_avg, uErr / Uinf_avg, 'bo');
hold on;
yline(1);
hold off;
xlabel('Vertical Distance in Wake ($y/D$)', 'Interpreter', 'latex');
ylabel('Normalized Mean Horizontal Velocity ($\bar{U} / U_{\infty}$)', 'Interpreter', 'latex');
title('Normalized Mean Horizontal Velocity vs. Vertical Distance', 'Interpreter', 'latex');
legend('95\% Confidence Interval', '$\bar{U} / U_{\infty}$ = 1', 'Interpreter', 'latex');
grid on;

%% Plot Figure 1b: Normalized Turbulence Intensity
figure;
plot(yD_numeric, uStd ./ umean, 'bo');

xlabel('Vertical Distance in Wake ($y/D$)', 'Interpreter', 'latex');
ylabel('Normalized Turbulence Intensity ($TI = u^{\prime} / \bar{U}$)', 'Interpreter', 'latex');
title('Normalized Turbulence Intensity vs. Vertical Distance', 'Interpreter', 'latex');
legend('Turbulence Intensity', 'Interpreter', 'latex');
grid on;

%% calculate Pinf data

yD_up = [0, 1, 2, 3, 4];

% Loop through all data files for U_inf
for i = 1:length(yD_up)
    % Construct filename
    FileName = ['Pinf_yD', num2str(yD_up(i))];

    % Read and convert pressure data
    Delta_p = readmatrix(FileName) * P_convert;

    % Compute U_inf values
    Pinf(i) = mean((P_atm - Delta_p));
end
Pinf_avg = mean(Pinf);

%% Calculate Pwake
yD = ["00", "01", "02", "03", "04", "05", "06", "08", "10", "12", ...
      "14", "16", "18", "20", "22", "24", "26", "30", "35", "40"];

% Preallocate arrays for performance
Pmean = zeros(size(yD));
PStd = zeros(size(yD));
PErr = zeros(size(yD));

for i = 1:length(yD)
    % Construct filename
    FileName = ['Pwake_yD' + yD(i)];

    % Read and convert pressure data
    Delta_p = readmatrix(FileName) * P_convert; 

    % Compute pressure wake values
    P = (P_atm - Delta_p);
    Pmean(i) = mean(P);
    PStd(i) = std(P);

    % Compute 95% confidence interval error bars
    PErr(i) = 2 * PStd(i) / sqrt(length(Delta_p)); 
end

% Convert yD to numeric values
yD_numeric = str2double(yD) / 10;

%% Perform integration
yUp = yD_up * D;
T1 = 2 * trapz(yD_up, Pinf);
T2 = 2 * trapz(yD_numeric, Pmean);
T3 = 2 * rho * trapz(yD_up, Uinfm.^2);
T4 = 2 * rho * trapz(yD_numeric, umean.^2);

Fd = T1 - T2 + T3 - T4;

Cd_calculated_P = Fd / ( (1/2) * rho * Uinf_avg.^2);

%% Calculate Pcylinder

% define array of polar angles examined in the experiment
theta = [0:5:90,100:20:180];

% total number of data files collected
N = length(theta);

% create arrays to store the mean and standard deviation
pmean = zeros(size(theta));
pstd = zeros(size(theta));

% loop through all data files
for (i=1:N)
    % filename of ith file
    FileName=['Pcyc_deg',num2str(theta(i))];

    % read in data from the file
    Delta_p = readmatrix(FileName) * P_convert;
    
    P_cyl = (P_atm - Delta_p);

    % calculate mean and std 
    P_cylm(i) = mean(P_cyl);
    p_cylstd(i) = std(P_cyl);
    P_cylErr(i) = 2 * p_cylstd(i) / sqrt(length(Delta_p)); 

end

Cp = (P_cylm - Pinf_avg) / ( (1/2) * rho * Uinf_avg.^2);
Cd_calc_momentum = trapz((theta .* pi / 180), Cp.*cos(theta .* pi ./ 180));
%% Plot Figure 1c: Coefficient of Pressure vs. Angular Position

figure;
errorbar(theta, Cp, P_cylErr, 'bo');
hold on;
yline(0);
hold off;
% Format the plot
xlabel('Angular Position ($\theta$) (degrees)', 'Interpreter', 'latex');
ylabel('Coefficient of Pressure ($C_p$)', 'Interpreter', 'latex');
title('Coefficient of Pressure vs. Angular Position', 'Interpreter', 'latex');
legend('95\% Confidence Interval','$C_p$ = 0', 'Interpreter', 'latex');

%% Plot Figure 1d: Drag Coefficient vs Reynolds Number
data = readmatrix('CD_RE_Textbook.dat'); % Read the published results

Cd = data(:,2); % Drag coefficient, C_d (vertical axis)
Re = data(:,1); % Reynolds number, Re_D (horizontal axis)

figure;
loglog(Re, Cd, 'k-', 'LineWidth', 1.5); % Published results as a solid black line
hold on;

% Plot Cd vs ReD based on static pressure analysis (solid red square)
loglog(Re_calc, Cd_calculated_P, 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 5);

% Plot Cd vs ReD based on momentum analysis (solid blue circle)
loglog(Re_calc, Cd_calc_momentum, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5);

% Format the plot
xlabel('Reynolds Number ($Re_D$)', 'Interpreter', 'latex');
ylabel('Drag Coefficient ($C_D$)', 'Interpreter', 'latex');
title('Drag Coefficient vs. Reynolds Number', 'Interpreter', 'latex');
grid on;

legend('Published Results', 'Static Pressure Analysis', 'Momentum Analysis', 'Location', 'Best');
hold off;

%% 2b

[maxUbar, idxUbar] = max((2 * abs(uStd))/umean);
maxUbar = maxUbar * 100;
maxUbarPos = yD_numeric(idxUbar);

disp('Max uncertainty in u_bar: ');
disp(maxUbar);
disp('Max uncertainty position in u_bar: ');
disp(maxUbarPos);

[maxPcyl, idxPcyl] = max((2 * abs(p_cylstd)) / (P_cylm - Pinf_avg));
maxPcyl = maxPcyl * 100;
maxPcylPos = theta(idxPcyl);

disp('Max uncertainty in Cp: ');
disp(maxPcyl);
disp('Max uncertainty position in Cp: ');
disp(maxPcylPos);

%% 2d.
interpolatedCd = interp1(Re(35:37), Cd(35:37), Re_calc);

errorCdMassMomentum = (abs(interpolatedCd - Cd_calc_momentum) / interpolatedCd) * 100;
errorCdStaticPressure = (abs(interpolatedCd - Cd_calculated_P) / interpolatedCd) * 100;

disp('The percent error between published vs momentum calc: ');
disp(errorCdMassMomentum);

disp('The percent error between published vs static pressure calc: ');
disp(errorCdStaticPressure);
