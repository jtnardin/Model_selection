%written 3-14-17 by JTN

clear all; clc

welllet = 'F';
well = 2;

Dtype = 'V';

%load in data and best-fits
load('cell_data_1d_struct_mod.mat')
load(['MLE_EST_' welllet num2str(well) '_D' Dtype '_final.mat'])


data = cell_data_1d_mod{5,well-1}';


N = numel(data);


q_length = zeros(length(q_final),1);

for i = 1:length(q_final)
    q_length(i) = length(q_final{i});
end



%Write parameter table
clear input

input.data = zeros(2,q_length(1) + q_length(2));

count = 1;

for i = 1:3
        input.data(i,1:q_length(1)) = q_final{count};
        input.data(i,q_length(1)+1:q_length(1) + q_length(4)) = q_final{count+3};
        count = count+1;
end

input.tableColLabels = {'D','v','a','D','$v_1$','$v_2$','$v_3$','$v_4$','$v_5$','a'};
input.tableRowLabels = {'$x_n=25$','$x_n=50$','$x_n=100$'};


% we want a complete LaTex document
input.makeCompleteLatexDocument = 1;

% generate LaTex code
latex = latexTable(input);

% save LaTex code as file
fid=fopen([welllet num2str(well) '_params_WLS.tex'],'w');
[nrows,ncols] = size(latex);
for row = 1:nrows
    fprintf(fid,'%s\n',latex{row,:});
end
fclose(fid);
% fprintf('\n... your LaTex code has been saved as ''MyLatex.tex'' in your working directory\n');
