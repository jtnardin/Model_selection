%written 4-5-17 by JTN to look at mean cell profiles \pm their SD's

i=4;
j=2;


[xn,tn] = size(mean_cell_data{i,j});
X = [1:xn,fliplr(1:xn)];

for k = 1:tn

    figure(1)
    
    Y = [mean_cell_data{i,j}(:,k)+mean_cell_data_sv{i,j}(:,k);flipud(mean_cell_data{i,j}(:,k)-mean_cell_data_sv{i,j}(:,k))];
    
    fill(X,Y,'b')
    hold on
    plot(mean_cell_data{i,j}(:,k))
    
    hold off

        
    figure(2)
    plot(cell_data_1d_mod{(j-1)*3+1,i}(:,k))
    hold on
    plot(cell_data_1d_mod{(j-1)*3+2,i}(:,k))
    plot(cell_data_1d_mod{(j-1)*3+3,i}(:,k))
    hold off

    pause(.125)

   
end
