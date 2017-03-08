%main_file_sim_janus_D0 2-7-17 by JTN to perform OLD
%optimization with D = 0

% clear all; clc


welllet = 'F';
well = 2;


%load data
load('cell_data_1d_struct_mod.mat')


%load in data
data = cell_data_1d_mod{5,well-1}';

xdata = linspace(0,1,540);
tdata = 0:1/3:1/3*(size(data,1)-1);

   
xnsize = [25 50 100];
model_sim = 1:2;

simnum = length(xnsize)*length(model_sim);

q_all = cell(simnum,1);
J_all = zeros(simnum,1);
WLS_weight_all = cell(simnum,1);
weight_matrix_all = cell(simnum,1);

count = 1;

for i = 1:length(xnsize)
    for j = 1:length(model_sim)

        %generate grids for computation
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

        
        
        
        if j == 1 %u_t = Du_xx - v_x
            q0 = [.01,.01,1];
        elseif j == 2
            %parameter value
            q0 = [.01,.01*ones(1,5),1];
        end

        %sparse matrix as a function for computation

        
        A_pos = @(se,sw,D,v,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[(-D+-v+v.*sw/2); ...
            (2*D+v-v.*se/2-v.*sw/2); (-D+v.*se/2)],xn,xn);

        A_neg = @(se,sw,D,v,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[(-D-v.*sw/2); ...
            (2*D+v.*se/2+v.*sw/2-v); (-D+v-v.*se/2)],xn,xn);


        options = optimset('display','iter','maxiter',100);

        tic

        [q,J,WLS_weight,weight_matrix] = fmincon(@(q) MLE_cost_DV(data,q,x,dx,xn,x_int,xbd_0,xbd_1,t,dt,tn,tdata,xdata,LE_loc,...
                     BC_x_0,BC_x_1,A_pos,A_neg,j),q0,[],[],[],[],[zeros(1,length(q0)-1),.7],[inf*ones(1,length(q0)-1),1.2],[],options);

        toc
        
        q_all{count} = q;
        J_all(count) = J;
        WLS_weight_all{count} = WLS_weight;
        weight_matrix_all{count} = weight_matrix;
        save(['/lustre/janus_scratch/jona8898/spline_vel/MLE_EST_' welllet num2str(well) '_OLS_DV.mat' ],...
            'q_all','J_all','WLS_weight_all','weight_matrix_all')
    
        count = count+1;
        
    end
end
