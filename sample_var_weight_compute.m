%compute the sample variance and weight matrix for a given simulations

clear all; clc

welllet = 'F';
well = 5;

Dtype = '0';

%load in data and best-fits
load('cell_data_1d_struct_mod.mat')
load(['MLE_EST_' welllet num2str(well) '_WLS_D' Dtype '.mat'])
% load(['WLS_weight_estimate' welllet num2str(well) '.mat'])


%load in data
data = cell_data_1d_mod{5,well-1}';
xdata = linspace(0,1,540);
tdata = 0:1/3:1/3*(size(data,1)-1);

xnsize = [25 50 100];
model_sim = 1:2;

count = 1;

model_sims = cell(6,1);

WLS_SV_all = zeros(6,1);
WLS_weight_all = cell(6,1);
weight_matrix_all = cell(6,1);
res_all = cell(6,1);

%model simulations


for i = 1:3
    for j = 1:2


        q = q_all{count};

        xn = xnsize(i);
        dt = 1e-3;

        [x,t] = grid_generate(xn,xdata(1),xdata(end),dt,tdata(1),tdata(end));
        tn = length(t);
        dx = x(2)-x(1);
        %interior, boundary points
        [x_int,xbd_0,xbd_1] = int_bd_def(xn);


        %initial condition
        %IC = interp1(xdata,smooth(data(1,:)),x,'spline');

        LE_loc = leading_edge_calc(smooth(data(1,:)),xdata,0.8,0);

        %boundary conditions
        BC_x_0 = @(t) 1;
        BC_x_1 = @(t) 0;


        %sparse matrix as a function for computation

        if strcmp(Dtype,'V')

            A_pos = @(se,sw,D,v,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[(-D+-v+v.*sw/2); ...
                (2*D+v-v.*se/2-v.*sw/2); (-D+v.*se/2)],xn,xn);

            A_neg = @(se,sw,D,v,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[(-D-v.*sw/2); ...
                (2*D+v.*se/2+v.*sw/2-v); (-D+v-v.*se/2)],xn,xn);
    
        elseif strcmp(Dtype,'0')

            A_pos = @(se,sw,v,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[(-v+v.*sw/2); ...
                (v-v.*se/2-v.*sw/2); (v.*se/2)],xn,xn);

            A_neg = @(se,sw,v,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[(-v.*sw/2); ...
                (v.*se/2+v.*sw/2-v); (v-v.*se/2)],xn,xn);
            
        else
            error('incorrect Dtype variable')
            
        end
    
            
       tic

       
       if strcmp(Dtype,'0')
            [J,WLS_SV_all(count),WLS_weight_all{count},weight_matrix_all{count},res_all{count}] ...
                = MLE_cost_D0(data,q,x,dx,xn,x_int,xbd_0,xbd_1,...
                t,dt,tn,tdata,xdata,LE_loc,BC_x_0,BC_x_1,A_pos,A_neg,j);
       elseif strcmp(Dtype,'V')
           [J,WLS_SV_all(count),WLS_weight_all{count},weight_matrix_all{count}] ...
                = MLE_cost_DV(data,q,x,dx,xn,x_int,xbd_0,xbd_1,...
                t,dt,tn,tdata,xdata,LE_loc,BC_x_0,BC_x_1,A_pos,A_neg,j);
       end
       
       
       
       toc     

       count = count + 1;
   
    end
end

save(['MLE_EST_' welllet num2str(well) '_WLS_D' Dtype '.mat'],'WLS_SV_all','WLS_weight_all','res_all','weight_matrix_all','-append')
