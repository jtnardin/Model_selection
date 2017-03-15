%compute AIC
%N is the number of data points, J_all is the cost vector, and q_length is
%a vector giving the number of parameters estimated in each model. 
function AIC = AIC_compute_WLS_f(N,J_all,WLS_weights,res,q_length)

    AIC = zeros(length(J_all),1);

    
    
    for i = 1:length(J_all)

        Jwls = 1/N*sum((res{i}./WLS_weights{i}).^2);
        
        AIC(i) = N*(1+log(2*pi)) + 2*N*log(Jwls)+ 2*sum(WLS_weights{i}) + 2*(q_length(i)+12);

    end

end