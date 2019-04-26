% How to find "Adjust" parameter for aligning two signals?
% 1) Choose a corresponding point with respect to "Amplitude" in both signals
% 2) Find their indices in the vectors of signals
% 3) Find the data for those indices in the time vectors
% 4) Calculate the difference between the two time elements

close all; clear all; clc;

load('C:\Users\shaya\Desktop\Third_Work\ADC Data\RAMP_ALL.mat');

E1_Time = 1e6 * E1_Time;

E2_AOUT_Defense_Time = 1e6 * E2_AOUT_Defense_Time;

E2_AOUT_Orig_Time = 1e6 * E2_AOUT_Orig_Time;

E2_AOUT_Trojan_Time = 1e6 * E2_AOUT_Trojan_Time;

E2_Input_Time = 1e6 * E2_Input_Time;

% Plot 1

Adjust_1 = E1_Time(74126) - E1_Time(72178);

figure();
plot(E1_Time + Adjust_1,E1_Input,'--k','LineWidth',1.5);
hold on;
plot(E1_Time,E1_AOUT_Orig,'g','LineWidth',1.5);
hold off;

ylabel('Signal Amplitude (V)', 'FontSize', 14);
xlabel(['Time (','\mu','s)'], 'FontSize', 14);
set(gca, 'fontsize', 12);
AX=legend(sprintf('Signal Generator'),sprintf('Healthy ADC'),'Location','SouthEast');
LEG = findobj(AX,'type','text');
set(LEG,'FontSize',12);
xlim([0 95]);
ylim([-0.31,0.31]);

% Plot 2 (Red)

Adjust_2 = E1_Time(74126) - E1_Time(72178);

figure();
plot(E1_Time + Adjust_2,E1_Input,'--k','LineWidth',1.5);
hold on;
plot(E1_Time,E1_AOUT_Trojan,'r','LineWidth',1.5);
hold off;

ylabel('Signal Amplitude (V)', 'FontSize', 14);
xlabel(['Time (','\mu','s)'], 'FontSize', 14);
set(gca, 'fontsize', 12);
AX=legend(sprintf('Signal Generator'),sprintf('ADC + Attack 1'),'Location','SouthEast');
LEG = findobj(AX,'type','text');
set(LEG,'FontSize',12);
xlim([0 95]);
ylim([-0.31,0.31]);

% Plot 3 (Blue)

Adjust_3 = E1_Time(77454) - E1_Time(72178);

figure();
plot(E1_Time + Adjust_3,E1_Input,'--k','LineWidth',1.5);
hold on;
plot(E1_Time,E1_AOUT_Defense,'b','LineWidth',1.5);
hold off;

ylabel('Signal Amplitude (V)', 'FontSize', 14);
xlabel(['Time (','\mu','s)'], 'FontSize', 14);
set(gca, 'fontsize', 12);
AX=legend(sprintf('Signal Generator'),sprintf('ADC + Attack 1 + Defense 1'),'Location','SouthEast');
LEG = findobj(AX,'type','text');
set(LEG,'FontSize',12);
xlim([0 95]);
ylim([-0.31,0.31]);

% Plot 4 (Red)

Adjust_4 = E2_AOUT_Trojan_Time(93558) - E2_Input_Time(62804);

figure();
plot(E2_Input_Time + Adjust_4,E2_Input_Data,'--k','LineWidth',1.5);
hold on;
plot(E2_AOUT_Trojan_Time,E2_AOUT_Trojan_Data,'r','LineWidth',1.5);
hold off;

ylabel('Signal Amplitude (V)', 'FontSize', 14);
xlabel(['Time (','\mu','s)'], 'FontSize', 14);
set(gca, 'fontsize', 12);
AX=legend(sprintf('Signal Generator'),sprintf('ADC + Attack 2'),'Location','SouthEast');
LEG = findobj(AX,'type','text');
set(LEG,'FontSize',12);
xlim([0 95]);
ylim([-0.31,0.31]);

% Plot 5 (Blue)

Adjust_5 = E2_AOUT_Defense_Time(120612) - E2_Input_Time(62804);

figure();
plot(E2_Input_Time + Adjust_5,E2_Input_Data,'--k','LineWidth',1.5);
hold on;
plot(E2_AOUT_Defense_Time,E2_AOUT_Defense_Data,'b','LineWidth',1.5);
hold off;

ylabel('Signal Amplitude (V)', 'FontSize', 14);
xlabel(['Time (','\mu','s)'], 'FontSize', 14);
set(gca, 'fontsize', 12);
AX=legend(sprintf('Signal Generator'),sprintf('ADC + Attack 2 + Defense 2'),'Location','SouthEast');
LEG = findobj(AX,'type','text');
set(LEG,'FontSize',12);
xlim([0 95]);
ylim([-0.31,0.31]);

% All Plots

figure();
plot(E1_Time + Adjust_2,E1_Input,'--k','LineWidth',1.5);
hold on;
plot(E1_Time,E1_AOUT_Trojan,'r','LineWidth',1.5);
hold off;
ylabel('Signal Amplitude (V)', 'FontSize', 14);
xlabel(['Time (','\mu','s)'], 'FontSize', 14);
set(gca, 'fontsize', 12);
title('ADC + Attack 1','FontSize', 14);
xlim([20 40]);
ylim([-0.22 -0.02]);

figure();
plot(E1_Time + Adjust_3,E1_Input,'--k','LineWidth',1.5);
hold on;
plot(E1_Time,E1_AOUT_Defense,'b','LineWidth',1.5);
hold off;
ylabel('Signal Amplitude (V)', 'FontSize', 14);
xlabel(['Time (','\mu','s)'], 'FontSize', 14);
set(gca, 'fontsize', 12);
title('ADC + Attack 1 + Defense 1','FontSize', 14);
xlim([20 40]);
ylim([-0.22 -0.02]);

figure();
plot(E2_Input_Time + Adjust_4,E2_Input_Data,'--k','LineWidth',1.5);
hold on;
plot(E2_AOUT_Trojan_Time,E2_AOUT_Trojan_Data,'r','LineWidth',1.5);
hold off;
ylabel('Signal Amplitude (V)', 'FontSize', 14);
xlabel(['Time (','\mu','s)'], 'FontSize', 14);
set(gca, 'fontsize', 12);
title('ADC + Attack 2','FontSize', 14);
xlim([20 40]);
ylim([-0.22 -0.02]);

figure();
plot(E2_Input_Time + Adjust_5,E2_Input_Data,'--k','LineWidth',1.5);
hold on;
plot(E2_AOUT_Defense_Time,E2_AOUT_Defense_Data,'b','LineWidth',1.5);
hold off;
ylabel('Signal Amplitude (V)', 'FontSize', 14);
xlabel(['Time (','\mu','s)'], 'FontSize', 14);
set(gca, 'fontsize', 12);
title('ADC + Attack 2 + Defense 2','FontSize', 14);
xlim([20 40]);
ylim([-0.22 -0.02]);
