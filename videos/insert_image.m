%insert_image.m written 4-5-17 by JTN to import data from experimental
%images and create 1d cell profiles. From these, we then find the average
%cell profile from triplicate data


clear all; clc


%initialize data cell
cell_data_1d_mod = cell(6,8);

x = linspace(0,1,540);

%to determine regions of interest -- calculated by eye from experimental
%videos to see what regions are relevant, and for what time period.
tend = 144*ones(6,8);
xstart = 1*ones(6,8);
xend = 540*ones(6,8);

yend = 540*ones(6,8);

xstart(4,1) = 320;
xstart(4,2) = 311;
xstart(4,4) = 148;
xstart(5,2) = 230;
xstart(5,3) = 153;
xstart(6,2) = 37;
xstart(6,3) = 115;
xstart(6,4) = 39;

xend(4,2) = 505;
xend(4,4) = 503;

yend(4,2) = 500;
yend(4,3) = 472;
yend(5,4) = 490;

tend(4,1) = 125;
tend(4,3) = 112;
tend(5,3) = 131;
tend(5,4) = 104;
tend(6,2) = 123;
tend(6,3) = 112;
tend(6,4) = 112;

%loop through each video and get relevant data
for i = 1:6
    for j = 1:8
        %load in corresponding video
        fname = ['WL_2_' char(65+i) '0' num2str(j+1) '.tif'];

        if i == 5 && j == 3 %had to modify this video -- lots of noise in wound
            fname = ['WL_2_' char(65+i) '0' num2str(j+1) '_concat.tif'];
        end
            
        %retrive info of image
        info = imfinfo(fname);

        %number of images
        num_images = tend(i,j);

        %xnum, ynum correspond to height and width of video resp.
        xnum = info.Width;
        ynum = xend(i,j) - xstart(i,j) + 1;

        %initialize matrices for videos
        A = zeros(ynum,xnum,num_images);

        B = zeros(540,540);
        
        for k = 1:num_images
            B = imread(fname,k,'info',info);
            A(:,:,k) = B(xstart(i,j):xend(i,j),:);
        end
        
        %ignore errant cells in wound space
        if yend(i,j)~=540
            A(:,yend(i,j):540,:)=0; 
        end
                
        %normalizing factor gotten from data in back of sheet
        img_mean = squeeze(mean(mean(A(:,1:100,:))));

        
        %now put 1d info into big vector
        %initialize entry in data cell
        cell_data_1d_mod{i,j} = zeros(xnum,num_images);
        
        
        %remove any noise past the LE
        % do so by defining cutoff density
        if i == 5 && j ==3
            LE_dens = 0.06; %this sample is particularly noisy
        else
            LE_dens = 0.02;
        end

        for k = 1:num_images
            %normalized cell density
            cell_data_1d_mod{i,j}(:,k) = mean(A(:,:,k))/img_mean(k);
            %find where cell dens small
            LE = leading_edge_calc(cell_data_1d_mod{i,j}(:,k),x,LE_dens,0);
            %set those points to zero
            cell_data_1d_mod{i,j}(x>LE,k) = 0; 
        end
        
    end
end

%Now create mean cell profiles
mean_cell_data = cell(8,2);
mean_cell_data_sv = cell(8,2);

for i = 1:8
    for j = 1:2
        
        %figure out minimum duration for triplicate data
        tend_mean = min(tend((j-1)*3+1:j*3,i));
        
        LE_init = zeros(3,1);
        
        %find three initial LE locations
        for k = 1:3
            LE_init(k) = leading_edge_calc(cell_data_1d_mod{(j-1)*3+k,i}(:,1),1:540,0.5,0);
        end
        
        %sort by initial LE location for alignment
        [LE_init,I]  =sort(LE_init);
        %align back of sheets to smallest initial LE
        cell_data_1d_mod{(j-1)*3+I(2),i}(1:LE_init(2) - LE_init(1) ,:) = [];
        cell_data_1d_mod{(j-1)*3+I(3),i}(1:LE_init(3) - LE_init(1) ,:) = [];
        %align front of sheets to largest initial LE
        cell_data_1d_mod{(j-1)*3+I(2),i}(end-(LE_init(3) - LE_init(2)) + 1:end,:) = [];
        cell_data_1d_mod{(j-1)*3+I(1),i}(end-(LE_init(3) - LE_init(1)) + 1:end,:) = [];
        
        %how many spatial points do we have now
        xend = size(cell_data_1d_mod{(j-1)*3+1,i},1);
        
        %initialize array for taking the mean
        temp_data = zeros(3,xend);
        
        %initialize mean data entry
        mean_cell_data{i,j} = zeros(xend,tend_mean);
        mean_cell_data_sv{i,j} = zeros(xend,tend_mean);
        
        %loop through each time point and find mean cell profile
        for k = 1:tend_mean
            for l = 1:3
                temp_data(l,:) = smooth(cell_data_1d_mod{(j-1)*3+l,i}(:,k))';
            end
            
            mean_cell_data{i,j}(:,k) = mean(temp_data);
            mean_cell_data_sv{i,j}(:,k) = sqrt(var(temp_data)); %standard deviation
            
        end
        
    end
end

%save
save('mean_cell_prof_data.mat','mean_cell_data','mean_cell_data_sv')

