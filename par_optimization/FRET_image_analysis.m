%FRET_image_analysis written 3-7-17 by JTN
% to quantify the FRET RAtios from experimental
% data and compare with best-fit velocities

if ~exist('EGF_wt_cell')
    load('FRET_imaging.mat','EGF_wt_cell')
end

load('cell_data_1d_struct_mod.mat')

welllet = 'F';

FRET_mean = cell(3,4);
FRET_mean_behind_LE = cell(3,4);

x = linspace(0,1,540);


%generate FRET cell from all experiments.
for i = 1:3
    for j = 1:4
    
        %loop through each time point and take mean over y-direction to get
        %a sense of  1D FRET wave . Then get average cellular FRET level
        for k = 1:144
            FRET_mean{i,j}(:,k) = mean(EGF_wt_cell{i,j}(:,:,k),2);
            %calculate LE location
            LE = leading_edge_calc(FRET_mean{i,j}(:,k),x,0.5,1);
            %find data points behind LE
           behind_LE = (x <= LE);
           %find mean value of data behind LE
            FRET_mean_behind_LE{i,j}(k) = mean(FRET_mean{i,j}(behind_LE,k));
            
        end
    end
end

%now plot + compare to best-fit v(t) or something

figure
set(gcf,'units','normalized','outerposition',[0 0 1 1])

cell_density = [1700 2500 3000 4000];


%loop through the different densities
for i = 1:4

        load(['MLE_EST_' welllet num2str(i+1) '_WLS_D0.mat'])

        %number of velocity points estimated
        q = q_all{6};
        n = length(q)-1;
    
        %data range considered for this data set
        tdata = 0:1/3:1/3*(size(cell_data_1d_mod{5,i},2)-1);
    
        %create v-spline function
        tsamp = augknt([tdata(1) tdata(end) tdata(round(linspace(1,length(tdata),n)))],2);
        v_spline = spmak(tsamp,q(1:end-1));
        V = @(t) fnval(v_spline,t);
        
        subplot(2,2,i)
        
        plot(tdata,V(tdata)/max(V(tdata)))
        hold on
        plot(tdata,FRET_mean_behind_LE{2,i}(1:length(tdata))/max(FRET_mean_behind_LE{2,i}))
        
        xlabel('Time')
        ylabel('Relative speed / FRET')
        title(['FRET and v(t), density = ' num2str(cell_density(i)) ' $cells/mm^2$'],'interpreter','latex')

end

exportfig(gcf,'FRET_v_speed.eps','color','rgb')
saveas(gcf,'FRET_v_speed.fig')

figure
set(gcf,'units','normalized','outerposition',[0 0 1 1])
figure
set(gcf,'units','normalized','outerposition',[0 0 1 1])

cell_density = [1700 2500 3000 4000];


%loop through the different densities
for i = 1:4

        load(['MLE_EST_' welllet num2str(i+1) '_WLS_D0.mat'])

        %number of velocity points estimated
        q = q_all{6};
        n = length(q)-1;
    
        %data range considered for this data set
        tdata = 0:1/3:1/3*(size(cell_data_1d_mod{5,i},2)-1);
    
        %create v-spline function
        tsamp = augknt([tdata(1) tdata(end) tdata(round(linspace(1,length(tdata),n)))],2);
        v_spline = spmak(tsamp,q(1:end-1));
        V = @(t) fnval(v_spline,t);
        
        %plot mean FRET
        figure(1)
        subplot(2,2,i)
        
        plot(tdata,V(tdata)/max(V(tdata)))
        hold on
        plot(tdata,FRET_mean_behind_LE{2,i}(1:length(tdata))/max(FRET_mean_behind_LE{2,i}))
        
        xlabel('Time')
        ylabel('Relative speed / FRET')
        title(['FRET and v(t), density = ' num2str(cell_density(i)) ' $cells/mm^2$'],'interpreter','latex')

        figure(2)
        subplot(2,2,i)
        
        
         surf(tdata,x,FRET_mean{2,i}(:,1:length(tdata)),'edgecolor','none')
         view(2)
         colorbar
         xlabel('Time')
         ylabel('Space')
         
         title(['FRET Ratio for density = ' num2str(cell_density(i)) ' $cells/mm^2$'],'interpreter','latex')

end

% figure(1)
% exportfig(gcf,'FRET_v_speed.eps','color','rgb')
% saveas(gcf,'FRET_v_speed.fig')

figure(2)
exportfig(gcf,'FRET.eps','color','rgb')
saveas(gcf,'FRET.fig')