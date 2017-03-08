function [J] = MLE_cost_DV_set_weights(data,q,x,dx,xn,x_int,xbd_0,xbd_1,t,dt,tn,tdata,xdata,LE_loc,...
        BC_x_0,BC_x_1,A_pos,A_neg,modelno)

    N = numel(data);
    weight_f = [1 5*ones(1,9) 1];
    
    
    model = diff_conv_f_vel_spline_OLS_DV(q,x,dx,xn,x_int,xbd_0,xbd_1,t,dt,tn,tdata,xdata,LE_loc,...
        BC_x_0,BC_x_1,A_pos,A_neg,modelno);
    
    %a vector full of the weights that correspond to points
    weight_matrix = zeros(numel(data),1);
    
    %calculate residuals
    res = model (:) - data(:);
    
    for j = 1:11
            %where is model between different densities
            model_subset = (model(:) >= .1*(j-1))&(model(:) < .1*(j)); 
            %just set the weight
            weight_matrix(model_subset) = weight_f(j);
    end

    %     now calculate the likelihood function
    J = sum((res./weight_matrix).^2); %sum(log(sqrt(WLS_SV)*weight_matrix) + 


end