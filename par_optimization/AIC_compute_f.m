%compute AIC
%N is the number of data points, J_all is the cost vector, and q_length is
%a vector giving the number of parameters estimated in each model. 
function AIC = AIC_compute_f(N,J_all,q_length)

    AIC = zeros(length(J_all),1);

    for i = 1:length(J_all)

        AIC(i) = N*((1+log(2*pi)) + 2*log(J_all(i)/N)) + 2*q_length(i);

    end

end