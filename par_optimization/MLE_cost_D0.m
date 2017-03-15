function [J,WLS_SV,WLS_weight,weight_matrix,res,model] = MLE_cost_D0(data,q,x,dx,xn,x_int,xbd_0,...
    xbd_1,t,dt,tn,tdata,xdata,LE_loc,BC_x_0,BC_x_1,A_pos,A_neg,modelno)

    weight_f = zeros(11,1);
    N = numel(data);
    
    model = diff_conv_f_vel_spline_OLS_D0(q,x,dx,xn,x_int,xbd_0,xbd_1,t,dt,tn,tdata,xdata,LE_loc,...
        BC_x_0,BC_x_1,A_pos,A_neg,modelno);
    
    %keeps track of which data entries correspond to which weights
    %(1,2,...,11)
    model_subset_ind = zeros(numel(data),1);
    %a vector full of the weights that correspond to points
    weight_matrix = zeros(numel(data),1);
    
    %calculate residuals
    res = model (:) - data(:);
    
    OLS_SV = 1/N*sum(res.^2); %slightly biased sample variance
    
    for j = 1:11
            %where is model between different densities
            model_subset = (model(:) >= .1*(j-1))&(model(:) < .1*(j)); 
            %take corresponding subset of res
            res_subset = res(model_subset); 
            %calculate weight for that section
            weight_f(j) = sqrt((1/numel(res_subset)*sum(res_subset.^2))/OLS_SV); 
            
            %keep track of which points correspond to which weights
            model_subset_ind(model_subset) = j;
            %put weights into corresponding location in matrix
            weight_matrix(model_subset) = weight_f(j);
    end

    
    %Now we have OLS data, do 10 iterations to estimate WLS data
    WLS_weight = zeros(11,1);
    for i = 1:10
        
        %calculate WLS sample variance
         WLS_SV = 1/N*sum((res./weight_matrix).^2);
 
        %now estimate each weight
        for j = 1:11
            res_subset = res(model_subset_ind==j);

            WLS_weight(j) = sqrt((1/numel(res_subset)*sum(res_subset.^2))/WLS_SV); 

            %update weight vector
            weight_matrix(model_subset_ind == j) = WLS_weight(j);

        end
        
    end

    %     now calculate the likelihood function
    J = sum(log(sqrt(WLS_SV)*weight_matrix) + (res./weight_matrix).^2/(2*WLS_SV));


end