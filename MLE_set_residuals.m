%script MLE_residuals.m written 2-28-17 by JTN
%to look at residuals for MLE estimation

clear all; clc

welllet = 'F';
well = 3;

Dtype = 'V';

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
            [J,WLS_SV_all(count),WLS_weight_all{count},weight_matrix_all{count},res_all{count},model_sims{count}] ...
                = MLE_cost_D0(data,q,x,dx,xn,x_int,xbd_0,xbd_1,...
                t,dt,tn,tdata,xdata,LE_loc,BC_x_0,BC_x_1,A_pos,A_neg,j);
       elseif strcmp(Dtype,'V')
           [J,WLS_SV_all(count),WLS_weight_all{count},weight_matrix_all{count},res_all{count},model_sims{count}] ...
                = MLE_cost_DV(data,q,x,dx,xn,x_int,xbd_0,xbd_1,...
                t,dt,tn,tdata,xdata,LE_loc,BC_x_0,BC_x_1,A_pos,A_neg,j);
       end
       
       
       
       toc     

       count = count + 1;
   
    end
end



%now get weights and residuals


modelnum = [1 2 1 2 1 2 ];
xn_num = [25 25 50 50 100 100];

figure(1)
set(gcf,'units','normalized','outerposition',[0 0 1 1])
figure(2)
set(gcf,'units','normalized','outerposition',[0 0 1 1])

    for i = 1:6

        %get model and residuals   
          
        
        err = 'WLS';
        C = weight_matrix_all{i};
        res_WLS = (res_all{i}./C)/sqrt(2*WLS_SV_all(i));
        model = model_sims{i};
        
        
        figure(1)
        subplot(3,2,i)

        %figure 1 -- rj v. fj
        plot(model(:),res_WLS(:),'b.')

        xlabel('$f(t_j,\hat{q})$','interpreter','latex','fontsize',30)
        ylabel('$r_j$','interpreter','latex','fontsize',30)

        set(gca,'fontsize',20)

        title([welllet num2str(well) ' ' err ' Residuals, Model ' num2str(modelnum(i)) ', xn = ' num2str(xn_num(i))],'fontsize',30)
         
        axis([0 max(max(model_sims{i})) -5 5])
        

        figure(2)
        subplot(3,2,i)
        [X,T] = meshgrid(xdata,tdata);

        %generate 2d contours of rj vs of xj and tj

        surf(X,T,reshape(res_WLS,length(tdata),length(xdata)),'edgecolor','none')

        xlabel('x','fontsize',30)
        ylabel('t','fontsize',30)

        set(gca,'fontsize',20)

        title([welllet num2str(well) ' ' err ' Residuals, Model ' num2str(modelnum(i)) ', xn = ' num2str(xn_num(i))],'fontsize',30)
         axis([0 1 0 tdata(end)])

        caxis([-2 2])
        colorbar


    end

    
    figure(1)

    exportfig(gcf,[welllet num2str(well) '_MLE_Res_' err '_D' Dtype '.eps'],'color','rgb')
    saveas(gcf,[welllet num2str(well) '_MLE_Res_' err '_D' Dtype '.fig'])


    figure(2)

    exportfig(gcf,[welllet num2str(well) '_MLE_Res_contour_' err '_iterative_D' Dtype '.eps'],'color','rgb')
    saveas(gcf,[welllet num2str(well) '_MLE_Res_contour_' err '_iterative_D' Dtype '.fig'])


    
    
    %plot fits against data
    
    
count = 1;
        
for i = 30:30:90
    figure; hold on
    
%     nonzerodataWLS = dataWLS(i,:)~=0;
    
    
    plot(xdata,data(i,:),'.','color','k')
%     plot(xdata(nonzerodataWLS),data(i,nonzerodataWLS),'k.')
    plot(xdata,model_sims{1}(i,:),'b')
    plot(xdata,model_sims{2}(i,:),'color',[0 .5 0])
    plot(xdata,model_sims{3}(i,:),'r')
    plot(xdata,model_sims{4}(i,:),'c')
    plot(xdata,model_sims{5}(i,:),'y')
    plot(xdata,model_sims{6}(i,:),'m')
    
       
    
    if count == 1
       legend('data','v, xn = 25','v(t), xn = 25','v, xn = 50','v(t), xn = 50','v, xn = 100','v(t), xn = 100','location','southwest') 
    end
    
    xlabel('x')
    ylabel('u')
    title(['MLE WLS Model fits ' welllet num2str(well) ', t = ' num2str(round(tdata(i)))])
    
    exportfig(gcf,['F' num2str(well) '_fits_to_data_WLS_MLE_D' Dtype '_' num2str(count) '.eps'],'color','rgb')
    saveas(gcf,['F' num2str(well) '_fits_to_data_WLS_MLE_D' Dtype '_' num2str(count) '.fig'])
    
    count = count + 1;
end    
    
%plot velocities
figure
hold on


%plot constant velocities
if strcmp(Dtype,'0')
    plot([tdata(1) tdata(end)],[q_all{1}(1) q_all{1}(1)],'b')
    plot([tdata(1) tdata(end)],[q_all{3}(1) q_all{3}(1)],'r')
    plot([tdata(1) tdata(end)],[q_all{5}(1) q_all{5}(1)],'y')
elseif strcmp(Dtype,'V')
    plot([tdata(1) tdata(end)],[q_all{1}(2) q_all{1}(2)],'b')
    plot([tdata(1) tdata(end)],[q_all{3}(2) q_all{3}(2)],'r')
    plot([tdata(1) tdata(end)],[q_all{5}(2) q_all{5}(2)],'y')
end

n = 5;

tsamp = augknt([tdata(1) tdata(end) tdata(round(linspace(1,length(tdata),n)))],2);

colors = 'gcm';

for i = 1:3

    if strcmp(Dtype,'0')
        v_spline = spmak(tsamp,q_all{2*i}(1:end-1));
    elseif strcmp(Dtype,'V')
        v_spline = spmak(tsamp,q_all{2*i}(2:end-1));
    end
    
    V = @(t) fnval(v_spline,t);
    
    plot(tdata,V(tdata),colors(i))    
    
end

xlabel('t')
ylabel('v(t)')

title(['Velocities for fit models to data ' welllet num2str(well)])

legend('v, xn = 25','v, xn = 50','v, xn = 100','v(t), xn = 25','v(t), xn = 50','v(t), xn = 100','location','southeast') 


exportfig(gcf,['F' num2str(well) '_vel_MLE_WLS_D' Dtype '_' num2str(count) '.eps'],'color','rgb')
saveas(gcf,['F' num2str(well) '_vel_MLE_WLS_D' Dtype '_' num2str(count) '.fig'])
   
    
    
