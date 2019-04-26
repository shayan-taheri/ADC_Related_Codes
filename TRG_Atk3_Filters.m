% Filtering for Defense

close all; clear all; clc;

load('C:\Users\shaya\Desktop\Delta Sigma ADC\Code\TRG\Ref.mat');

load('C:\Users\shaya\Desktop\Delta Sigma ADC\Code\TRG\Atk3.mat');

Ref_Time = 1e6 * Ref_Time;

Atk3_Time = 1e6 * Atk3_Time;

% "Interpolation" for making the size of data vectors equal
% Ref_Sized = spline(Ref_Time,Ref_Vout,Atk3_Time); --> Not good for sizing 
Ref_Sized = interp1(Ref_Time,Ref_Vout,Atk3_Time); % Doing Interpolation

% Modified Residual: Reference Data - Filtered Data
% -------------------------------------------------------------------------

% General Filter Model: Atk3_Filtered = FUNCTION(Atk3_Vout);

% A) Simple Linear Regression
% Inputs: "Desired Data" (x) and "Test/Noisy Data" (y)
x = Ref_Sized; % Ref. Data in All Cases (except Case C)
y = Atk3_Vout; % Test/Noisy Data in All Cases
deg = 1; % Degree of Regression (N-Degree Polynomial)
p = polyfit(x,y,deg);
yfit = polyval(p,x);
yfit = yfit - mean(yfit); % Offset Elimination

% scale_ft = mean(Ref_Sized(1.6e5:1.8e5))/mean(yfit(1.6e5:1.8e5)); % OLD
% Scaling all data
scale_ft = 0.3/mean(abs(yfit(1.42e5:1.45e5))); % For Attack 3 Model (Clean Part)
yfit = yfit * scale_ft;

scale_ft = 0.3/mean(abs(x(1.6e5:1.8e5)));
x = x * scale_ft;
% ---- R-Square Calculation ----
yresid = x - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(x)-1) * var(x);
rsqA = 1 - (SSresid/SStotal);
yfitA = yfit;
% rsq_adj = 1 - (SSresid/SStotal)*(length(y)-1)/(length(y)-length(p));

% ------------------------------------------------------------------

% B) Generalized Linear Regression for Different Distribution
% Inputs: "Desired Data" (x) and "Test/Noisy Data" (y)
% Different Options for Distribution Argument
% and Their Corresponding Link Functions
% 1) 'gamma' --> 'reciprocal'
% 2) 'inverse gaussian' --> p (a number). Default: p = -2
% 3) 'normal' --> 'identity'
% 4) 'poisson' --> 'log'
% *) 'binomial' --> 'logit' , Note: It is used for binary data.
x = Ref_Sized; % Ref. Data in All Cases (except Case C)
y = Atk3_Vout; % Test/Noisy Data in All Cases
distribution = 'gamma'; % Type of Distribution
link_function = 'reciprocal'; % Corresponding Link Function
% ---- Preparing Data for Process ----
x = x';
y = y';
y = y + 1; % Prepare data for "Gamma" distribution
% ---- Fitting Data and Getting Values of Coefficients ----
coef_vals = glmfit(x,y,distribution);
% ---- Extracting Estimated Response Data ----
yfit = glmval(coef_vals,x,link_function);
yfit = yfit - 1; % Inverse operation due to distribution
yfit = yfit';
yfit = yfit - mean(yfit); % Offset Elimination

% Scaling all data
scale_ft = 0.3/mean(abs(yfit(1.42e5:1.45e5))); % For Attack 3 Model (Atk3)
yfit = yfit * scale_ft;

scale_ft = 0.3/mean(abs(x(1.6e5:1.8e5)));
x = x * scale_ft;
% ---- R-Square Calculation ----
yresid = x - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(x)-1) * var(x);
rsqB = 1 - (SSresid/SStotal);
yfitB = yfit;

% ------------------------------------------------------------------

% C) Moving Average Filter
% Inputs: "Desired Data" (z), "Test/Noisy Data" (y), and Its Time (x)
y = Atk3_Vout; % Test/Noisy Data in All Cases
x = Atk3_Time; % Time Data in Case C
z = Ref_Sized; % Ref. Data in Case C
wind_size = 15000; % Window Size for Moving Average
Initial_Time = x;
x = x'; % Forward data preparation
y = y'; % Forward data preparation
y = tsmovavg(y,'s',wind_size);
nan_ind = isnan(y);
y = y(~nan_ind);
x = x(~nan_ind);
x = x'; % Backward data preparation
x(1) = 0;
y = y'; % Backward data preparation
yfit = interp1(x,y,Initial_Time); % Doing Interpolation

for i = 1:length(yfit)
    if(yfit(i) > 0)
        yfit(i) = 0.3;
    else if (yfit(i) < 0)
        yfit(i) = -0.3;
        end
    end
end

% Scaling all data
scale_ft = 0.3/mean(abs(z(1.6e5:1.8e5)));
z = z * scale_ft;
% ---- R-Square Calculation ----
yresid = z - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(z)-1) * var(z);
rsqC = 1 - (SSresid/SStotal);
yfitC = yfit;

% ------------------------------------------------------------------

% D) Robust Regression Filter (LAR)
% Inputs: "Desired Data" (x) and "Test/Noisy Data" (y)
x = Ref_Sized; % Ref. Data in All Cases (except Case C)
y = Atk3_Vout; % Test/Noisy Data in All Cases
[sdf1,yinf1] = fit(x,y,'poly1','Robust','LAR');
yfit = feval(sdf1,x);

% Scaling all data
scale_ft = 0.3/mean(abs(yfit(1.42e5:1.45e5))); % For Attack 3 Model (Atk3)
yfit = yfit * scale_ft;

scale_ft = 0.3/mean(abs(x(1.6e5:1.8e5)));
x = x * scale_ft;
% ---- R-Square Calculation ----
yresid = x - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(x)-1) * var(x);
rsqD = 1 - (SSresid/SStotal);
yfitD = yfit;

% ------------------------------------------------------------------

% E) Robust Regression Filter (Bisquare)
% Inputs: "Desired Data" (x) and "Test/Noisy Data" (y)
x = Ref_Sized; % Ref. Data in All Cases (except Case C)
y = Atk3_Vout; % Test/Noisy Data in All Cases
[sdf2,yinf2] = fit(x,y,'poly1','Robust','Bisquare');
yfit = feval(sdf2,x);

yfit = yfit - mean(yfit); % Offset Elimination

% Scaling all data
scale_ft = 0.3/mean(abs(yfit(1.42e5:1.45e5))); % For Attack 3 Model (Atk3)
yfit = yfit * scale_ft;

scale_ft = 0.3/mean(abs(x(1.6e5:1.8e5)));
x = x * scale_ft;
% ---- R-Square Calculation ----
yresid = x - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(x)-1) * var(x);
rsqE = 1 - (SSresid/SStotal);
yfitE = yfit;

% ------------------------------------------------------------------

% F) Median Filtering
% Inputs: "Desired Data" (x) and "Test/Noisy Data" (y)
x = Ref_Sized; % Ref. Data in All Cases (except Case C)
y = Atk3_Vout; % Test/Noisy Data in All Cases
yfit = medfilt1(y);

% Scaling all data
scale_ft = 0.3/mean(abs(yfit(1.42e5:1.45e5))); % For Attack 3 Model (Atk3)
yfit = yfit * scale_ft;

scale_ft = 0.3/mean(abs(x(1.6e5:1.8e5)));
x = x * scale_ft;
% ---- R-Square Calculation ----
yresid = x - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(x)-1) * var(x);
rsqF = 1 - (SSresid/SStotal);
yfitF = yfit;

% ------------------------------------------------------------------

% G) Savitzky-Golay Filtering
% Inputs: "Desired Data" (x) and "Test/Noisy Data" (y)
x = Ref_Sized; % Ref. Data in All Cases (except Case C)
y = Atk3_Vout; % Test/Noisy Data in All Cases
% Inputs: "Desired Data" (x) and "Test/Noisy Data" (y)
framelen = 10001;
order = 1;
yfit = sgolayfilt(y,order,framelen);

for i = 1:length(yfit)
    if(yfit(i) > 0)
        yfit(i) = 0.3;
    else if (yfit(i) < 0)
        yfit(i) = -0.3;
        end
    end
end

% Scaling all data
scale_ft = 0.3/mean(abs(yfit(1.42e5:1.45e5))); % For Attack 3 Model (Atk3)
yfit = yfit * scale_ft;

scale_ft = 0.3/mean(abs(x(1.6e5:1.8e5)));
x = x * scale_ft;
% ---- R-Square Calculation ----
yresid = x - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(x)-1) * var(x);
rsqG = 1 - (SSresid/SStotal);
yfitG = yfit;

% ------------------------------------------------------------------

% H) ARMAX Filtering
% Inputs: "Desired Data" (x) and "Test/Noisy Data" (y)
x = Ref_Sized; % Ref. Data in All Cases (except Case C)
y = Atk3_Vout; % Test/Noisy Data in All Cases
% Inputs: "Desired Data" (x) and "Test/Noisy Data" (y)
% ---- Data Object Creation ----
data = iddata(y,x,1);
% ---- Best Value for Goodness of Fit ----
best_AR = 10;
best_X = 10;
best_MA = 10;
best_Delay = 0; % There is no delay between input and output (Assumption).
% ---- Getting System ----
system = armax(data,[best_AR,best_X,best_MA,best_Delay]);
% coef_vect = getpvec(system);
[yfit_obj,~,~] = compare(data,system);
yfit = yfit_obj.OutputData;

% Scaling all data
scale_ft = 0.3/mean(abs(yfit(1.42e5:1.45e5))); % For Attack 3 Model (Atk3)
yfit = yfit * scale_ft;

scale_ft = 0.3/mean(abs(x(1.6e5:1.8e5)));
x = x * scale_ft;

for i = 1:length(yfit)
    if(yfit(i) > 0)
        yfit(i) = 0.3;
    else if (yfit(i) < 0)
        yfit(i) = -0.3;
        end
    end
end

% ---- R-Square Calculation ----
yresid = x - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(x)-1) * var(x);
rsqH = 1 - (SSresid/SStotal);
yfitH = yfit;

% ------------------------------------------------------------------

% I) Ridge Regression
% Inputs: "Desired Data" (x) and "Test/Noisy Data" (y)
x = Ref_Sized; % Ref. Data in All Cases (except Case C)
y = Atk3_Vout; % Test/Noisy Data in All Cases
k = 1e-1; % 1e-10:1e-1:1e1
coeff = ridge(y,x,k,0);
yfit = coeff(2) * x + coeff(1);
yfit = yfit - mean(yfit); % Offset Elimination

% Scaling all data
scale_ft = 0.3/mean(abs(yfit(1.42e5:1.45e5))); % For Attack 3 Model (Atk3)
yfit = yfit * scale_ft;

scale_ft = 0.3/mean(abs(x(1.6e5:1.8e5)));
x = x * scale_ft;

yresid = x - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(x)-1) * var(x);
rsqI = 1 - (SSresid/SStotal);
yfitI = yfit;

% ------------------------------------------------------------------

% J) Critical Region Filtering
% Inputs: "Input Time" (Ref_Time) and "Input Signal" (Ref_Vs)
% Inputs: "Output Time" (Atk3_Time) and "Output Signal" (Atk3_Vout)
x = Ref_Sized; % Ref. Data in All Cases (except Case C)
y = Ref_Vs; % with "Ref_Time"
w = Atk3_Vout; % with "Atk3_Time"

y = interp1(Ref_Time,Ref_Vs,Atk3_Time); % Doing Interpolation

% Scaling all data before filtering
scale_ft = 0.3/mean(abs(w(1.42e5:1.45e5))); % For Attack 3 Model (Atk3)
w = w * scale_ft;
scale_ft = 0.3/mean(abs(x(1.6e5:1.8e5)));
x = x * scale_ft;

% Applying Filtering Process
for ind = 1:length(Atk3_Time)
    if (y(ind) >= 0.2)
        w(ind) = abs(min(y));
    else if (y(ind) <= -0.25)
        w(ind) = min(y);
        end
    end
    disp(['Completed Percentage = ',num2str(100*ind/length(Atk3_Time))]);
end

yfit = w;
yresid = x - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(x)-1) * var(x);
rsqJ = 1 - (SSresid/SStotal);
yfitJ = yfit;

% ------------------------------------------------------------------
% ------------------------------------------------------------------

% Plots for Checking Results
figure; plot(Ref_Sized);
figure; plot(yfitA);
figure; plot(yfitB);
figure; plot(yfitC);
figure; plot(yfitD);
figure; plot(yfitE);
figure; plot(yfitF);
figure; plot(yfitG);
figure; plot(yfitH);
figure; plot(yfitI);
figure; plot(yfitJ);

% ---- Desired Plotting Results ----
figure();
plot(x,y,'*b');
hold on;
plot(x,yfit,'-r');
xlabel('X', 'FontSize', 14);
ylabel('Y', 'FontSize', 14);
Rsquare = num2str(rsqXX);
title(['XX Filter Analysis - R^2 = ',Rsquare], 'FontSize', 14);
set(gca, 'fontsize', 12);
fig_leg = legend ('Actual Data','Estimated/Filtered Data');
set(fig_leg,'FontSize',10);
hold off;
