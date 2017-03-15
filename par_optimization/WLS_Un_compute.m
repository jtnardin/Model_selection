clear all; clc

format long

welllet = 'F';
 well = 2;

%load in data and best-fits
load('cell_data_1d_struct_mod.mat')
data = cell_data_1d_mod{5,well-1}';
n = numel(data); 


load(['MLE_EST_' welllet num2str(well) '_DV_final.mat'])
JD = zeros(6,1);

for i = 1:6
    JD(i) = sum((res_all{i}./weight_matrix_all{i}).^2);
end


load(['MLE_EST_' welllet num2str(well) '_D0_final.mat'])
JD0 = zeros(6,1);

for i = 1:6
    JD0(i) = sum((res_all{i}./weight_matrix_all{i}).^2);
end

%calculate U_n
Un_WLS = n*(JD0 - JD)./JD

chi2cdf(Un_WLS(1:2:end),2)
chi2cdf(Un_WLS(2:2:end),6)

%do the same for model 1 vs. model 2
Un_WLS_1v2 = n*(JD0(1:3) - JD0(4:6))./JD0(4:6)

chi2cdf(Un_WLS_1v2,2)

