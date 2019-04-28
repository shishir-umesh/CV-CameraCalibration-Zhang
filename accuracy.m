%We can get all the world coordinates as follows
grid_coordinates = [];

for i = 0:width
    for j = 0:height
        %We append each calculated grid coordinate to the grid_coordinates
        grid_coordinates = [grid_coordinates; i*30, j*30, 1];
    end
end
%We also have all the homographies for each image

for i = 1:4
    H = eval(['H' num2str(i)]);
    p_approx = H*grid_coordinates';
    
    %We normalise the points as follows
    for j = 1:length(p_approx)
        p_approx(:,j) = p_approx(:,j) / p_approx(3,j);
    end
    
    img = eval(['img' num2str(i)]);
    figure(), imshow(img)
    hold on
    title(['Figure 1 : Projected grid corners fro >> ' files(i)])
    plot(p_approx(1,:),p_approx(2,:),'ro', 'MarkerSize',3); 
    hold off
end