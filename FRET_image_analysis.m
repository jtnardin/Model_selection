%FRET_image_analysis written 3-7-17 by JTN
% to quantify the FRET RAtios from experimental
% data and compare with best-fit velocities

if ~exist('EGF_wt_cell')
    load('FRET_imaging.mat','EGF_wt_cell')
end

welllet = 'F';
well = 2;


FRET_data = EGF_wt_cell{double(welllet)-68,well-1};
FRET_mean = cell(3,4);
FRET_mean_behind_LE = cell(3,4);

x = linspace(0,1,540);

%generate FRET cell from all experiments.
for i = 1:3
    for j = 1:4
    
        %loop through each time point and take mean over y-direction to get
        %a sense of  1D FRET wave 
        for k = 1:144
            FRET_mean{i,j}(:,k) = mean(EGF_wt_cell{i,j}(:,:,k),2);
        end
    end
end


for i = 1:3
    for j = 1:4
            
        for k = 1:144
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