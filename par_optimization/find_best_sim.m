%find_best_sim.m written 3-13-17 by JTN to find
%from multiple best-fit simulations, the one that best-fit 
%the data

clear all; clc

welllet = 'F';
well = 4;
Dtype = '0';

%load in data
load('cell_data_1d_struct_mod.mat')

%choose data
data = cell_data_1d_mod{5,well-1}';
%data grids
xdata = linspace(0,1,540);
tdata = 0:1/3:1/3*(size(data,1)-1);

J_final = zeros(6,1);
q_final = cell(6,1);
WLS_weights_final = cell(6,1);
weight_matrix_final = cell(6,1);

for i = 1:6 %do for all models and grid sizes
    
    load(['MLE_EST_' welllet num2str(well) '_WLS_D' Dtype '_par_' num2str(i) '.mat'])
    
    best_fit = J_all == min(J_all);
    
    J_final(i) = J_all(best_fit);
    q_final{i} = q_all{best_fit};
    WLS_weights_final{i} = WLS_weight_all{best_fit};
    weight_matrix_final{i} = weight_matrix_all{best_fit};
       
end

save(['MLE_EST_' welllet num2str(well) '_D' Dtype '_final.mat'],...
    'J_final','q_final','WLS_weights_final','weight_matrix_final')