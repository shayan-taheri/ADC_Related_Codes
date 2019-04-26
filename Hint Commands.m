uiopen('C:\Users\shaya\Desktop\EXP_1_RAMP.csv',1) % Loading Data

plot(VTAOUTX,VTAOUTY); % A simple testing plot

sum(VTAOUTX - VTAOUT_DX); % Check to see if there is any difference in DATA

sum(AOUT_D_Time - AOUT_T_Time); % Check to see if there is any difference in TIME

% Second to Microsecond Conversion
Time_Unit = 1e6 * Time;

clc; % Clear window

% A series of testing differences in signals
A = sum(VTAOUT_DX - VTAOUT_TX); B = sum(VTAOUTX - vVIPresulttranvVINresulttranX); C = sum(VTAOUTX - VTAOUT_TX);

save('C:\Users\shaya\Desktop\ADC Data\EXP_1_RAMP.mat'); % Saving Data
clear all; % Clear all the variables
load('C:\Users\shaya\Desktop\ADC Data\EXP_1_RAMP.mat'); % Loading Data

% Find "Not A Number (NAN)" in signal
indices = find(isnan(AOUTX) == 1);
[I,J] = ind2sub(size(AOUTX),indices);

% Truncating the data vector (signal) for eliminating the NAN elements
AOUT_Orig_Data = AOUT_Orig_Data(1:372964);
AOUT_Orig_Time = AOUT_Orig_Time(1:372964);

% Processing Functions for Signal/Data
std(X)
mean(X)
gt(mean(Diff1),mean(Diff2))
crosscorr(X1,X2)
xcorr(X1,X2)

% How to find "Adjust" parameter for aligning two signals?
% 1) Choose a corresponding point with respect to "Amplitude" in both signals
% 2) Find their indices in the vectors of signals
% 3) Find the data for those indices in the time vectors
% 4) Calculate the difference between the two time elements

Adjust_1 = E1_Time(74126) - E1_Time(72178);
figure(); plot(E1_Time + Adjust_1,E1_Input,'k'); hold on; plot(E1_Time,E1_AOUT_Orig);

Adjust_5 = E2_AOUT_Defense_Time(120612) - E2_Input_Time(62804);
figure(); plot(E2_Input_Time + Adjust_5,E2_Input_Data,'k'); hold on; plot(E2_AOUT_Defense_Time,E2_AOUT_Defense_Data);