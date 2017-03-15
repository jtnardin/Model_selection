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


AIC = AIC_compute_WLS_f(N,J_final,weight_matrix_all,res_all,q_length);


%%%%% write AIC, J table
input.data = zeros(3,6);

count = 1;
for i = 1:2
    for j = 1:3
        input.data(j,3*(i-1)+1) = J_final(count);
        input.data(j,3*(i-1)+2) = WLS_SV_all(count);
        input.data(j,3*i) = AIC(count);
        count = count+1;
    end
end

input.tableColLabels = {'J','\sigma','AIC','J','\sigma','AIC'};
input.tableRowLabels = {'$x_n=25$','$x_n=50$','$x_n=50$'};


% we want a complete LaTex document
input.makeCompleteLatexDocument = 1;

% generate LaTex code
latex = latexTable(input);

% save LaTex code as file
fid=fopen([welllet num2str(well) '_JsigmaWLS.tex'],'w');
[nrows,ncols] = size(latex);
for row = 1:nrows
    fprintf(fid,'%s\n',latex{row,:});
end
fclose(fid);
% fprintf('\n... your LaTex code has been saved as ''MyLatex.tex'' in your working directory\n');
