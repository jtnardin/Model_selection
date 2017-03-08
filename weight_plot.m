
clear all; clc

welllet = 'F';
well = 5;

Dtype = '0';

%load in data and best-fits
load(['MLE_EST_' welllet num2str(well) '_WLS_D' Dtype '.mat'])

figure
hold on

y = 0:.1:1;

for i = 1:6
    
    plot(y,WLS_weight_all{i})
        
end

legend('v, xn = 25','v(t), xn = 25','v, xn = 50','v(t), xn = 50','v, xn = 100','v(t), xn = 100','location','northwest') 

xlabel('f(t,q)')
ylabel('$w_j$','interpreter','latex')

title(['Weights for all model sims, data = '  welllet num2str(well)])


    exportfig(gcf,[welllet num2str(well) '_weights_D' Dtype '.eps'],'color','rgb')
    saveas(gcf,[welllet num2str(well) '_weights_D' Dtype '.fig'])
    