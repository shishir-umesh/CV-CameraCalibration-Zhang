grid_coordinates = [];

for i = 0:width
    for j = 0:height
        %We append each calculated grid coordinate to the grid_coordinates
        grid_coordinates = [grid_coordinates; i*30, j*30, 1];
    end
end
%We also have all the homographies for each image

for i = 1:4
    
    %Obtain the corresponding homography of the image
    H = eval(['H' num2str(i)]);
    p_approx = H*grid_coordinates';
    
    %We normalise the points as follows
    for j = 1:length(p_approx)
        p_approx(:,j) = p_approx(:,j) / p_approx(3,j);
    end
    
    %Plot the approximate projected grid points onto the image
    img = eval(['img' num2str(i)]);
    figure(), imshow(img)
    hold on
    title(['Figure 1 : Projected grid corners for >> ' files(i)])
    plot(p_approx(1,:),p_approx(2,:),'ro', 'MarkerSize',3); 
    hold off
    
    %Obtain the harris corners for the image.
    sigma = 2;
    thresh = 500;
    radius = 2;
    [cim,r,c,rsubp,csubp]=harris(rgb2gray(img),sigma,thresh,radius,1);
    title(["Figure 2 : Harris corners for >> " files(1)])
    
    %Computing the closest points to the Harris corners of the image
    D = dist2(p_approx(1:2,:)',[csubp, rsubp]);
    [D_sorted, D_index] = sort(D, 2);
    p_correct(:,:,i) = [csubp(D_index(:,1)),rsubp(D_index(:,1)),ones((width+1)*(height+1),1)];
    figure(), imshow(img), title(["Figure3 : Grid Points for >> " files(i)])
    hold on
    plot(p_correct(:,1,i),p_correct(:,2,i),'g+')
    hold off

    %Use the newly found closest harris corners to compute a new homography
    H_new(:,:,i) = homography2d(grid_coordinates',p_correct(:,:,i)');
    H_new(:,:,i) = H_new(:,:,i)/H_new(3,3,i);
    disp(["New Homography H for >> " files(i)])
    disp(H_new(:,:,i))
end

Hnew1 = Hnew(:,:,1);
Hnew2 = Hnew(:,:,2);
Hnew3 = Hnew(:,:,3);
Hnew4 = Hnew(:,:,4);
