%main_file_compare_w_data_MLE.m written 2-28-17 by JTN to select data set to use, model, and
%parameter and send to cost function

clear all; clc

welllet = 'F';
 well = 2;

Dtype = '0';

%load in data and best-fits
load('cell_data_1d_struct_mod.mat')
load(['MLE_EST_' welllet num2str(well) '_WLS_D' Dtype '.mat'])

%load in data
data = cell_data_1d_mod{5,well-1}';
xdata = linspace(0,1,540);
tdata = 0:1/3:1/3*(size(data,1)-1);


xnsize = [25 50 100];
model_sim = 1:2;

count = 1;

