%This script visualizes the timing simulations and fits a square root function
%timePoints were collected by running NetworkTest.m 

%load data
timePoints = load('timePoints.mat');
timePoints = timePoints.timePoints;

%make histogram plot
figure
histogram(timePoints(:,end),15)

%make time variability plot
devs = std(timePoints,1);
means = mean(timePoints,1);
fit_coeff = sqrt(means)'\devs';
figure;
box off;
plot(means,devs,'*');
hold on
plot(0:0.01:630,fit_coeff*sqrt(0:0.01:630))
xlabel('Mean time [ms]', 'FontSize',25)
ylabel('Standard deviation [ms]', 'FontSize',25)
legend('simulation','fit')

%the rms error of the fit
RMSE = sqrt(mean((devs-fit_coeff*sqrt(means)).^2));

%try a linear fit to see how much larger the error is (constraint: has to
%start from the origin)
devs = std(timePoints,1);
means = mean(timePoints,1);
fit_coeff_lin = means'\devs';
figure;
box off;
plot(means,devs,'*');
hold on
plot(0:0.01:630,fit_coeff_lin*(0:0.01:630))
xlabel('Mean time [ms]', 'FontSize',25)
ylabel('Standard deviation [ms]', 'FontSize',25)
legend('data','fit')

%the rms error of the fit
RMSE_lin = sqrt(mean((devs-fit_coeff_lin*means).^2));
