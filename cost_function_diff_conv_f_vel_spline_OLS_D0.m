function cost = cost_function_diff_conv_f_vel_spline_OLS_D0(data,q,x,dx,xn,x_int,xbd_0,xbd_1,t,dt,tn,tdata,xdata,LE_loc,...
        BC_x_0,BC_x_1,A_pos,A_neg,modelno)

    
    %simulate model
    model = diff_conv_f_vel_spline_OLS_D0(q,x,dx,xn,x_int,xbd_0,xbd_1,t,dt,tn,tdata,xdata,LE_loc,...
        BC_x_0,BC_x_1,A_pos,A_neg,modelno);

    cost = sum(sum((model-data).^2));  
    
    

end