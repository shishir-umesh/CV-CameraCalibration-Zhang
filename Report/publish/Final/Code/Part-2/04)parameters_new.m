% We now use our new homographies to estimate our intrinsic and extrinsic
% parameters.

V = [];
for i = 1:4
    H_temp = eval(['Hnew' num2str(i)]);
    h1 = H_temp(:,1);
    h2 = H_temp(:,2);
    h3 = H_temp(:,3);
    
    v11 = [h1(1)*h1(1), h1(1)*h1(2)+h1(2)*h1(1), h1(2)*h1(2), h1(3)*h1(1)+h1(1)*h1(3), h1(3)*h1(2)+h1(2)*h1(3), h1(3)*h1(3)]';
    v12 = [h1(1)*h2(1), h1(1)*h2(2)+h1(2)*h2(1), h1(2)*h2(2), h1(3)*h2(1)+h1(1)*h2(3), h1(3)*h2(2)+h1(2)*h2(3), h1(3)*h2(3)]';
    v22 = [h2(1)*h2(1), h2(1)*h2(2)+h2(2)*h2(1), h2(2)*h2(2), h2(3)*h2(1)+h2(1)*h2(3), h2(3)*h2(2)+h2(2)*h2(3), h2(3)*h2(3)]';
    
    V = [V; v12'; (v11-v22)'];
end

[U, Sigma, V_transpose] = svd(V);

b = V_transpose(:,end);

B11 = b(1);
B12 = b(2);
B22 = b(3);
B13 = b(4);
B23 = b(5);
B33 = b(6);

B = [B11, B12, B13; B12, B22, B23; B13, B23, B33];

v0 = (B12*B13 - B11*B23)/(B11*B22 - B12^2);
lambda = B33 - (B13^2 + v0*(B12*B13-B11*B23))/B11;
alpha = sqrt(lambda/B11);
beta = sqrt(lambda*B11/(B11*B22-B12^2));
gamma = -B12*alpha^2*beta/lambda;
u0 = gamma*v0/alpha - B13*alpha^2/lambda;

%Therefore our intrinsic matrix A can be defined as follows,
A = [alpha, gamma, u0; 0, beta, v0; 0, 0, 1];

%Now we calculate the extrinsic parameters and store them for future use
for i = 1:4
    H_temp = eval(['Hnew' num2str(i)]);
    h1 = H_temp(:,1);
    h2 = H_temp(:,2);
    h3 = H_temp(:,3);
    
    lambda_r = 1/ norm(A\h1);
    r1 = lambda_r*(A\h1);
    r2 = lambda_r*(A\h2);
    r3 = cross(r1,r2);
    t(:,i)  = lambda_r*(A\h3);
    
    R = [r1, r2, r3];
    
    [U,S,Vprime] = svd(R);
    Rotation(:,:,i) = U*Vprime;
    
    disp(["Rotation matrix R for images" files(i)])
    disp(Rotation(:,:,i))
    disp(["Translation vector for images" files(i)])
    disp(t(:,i))
    
    %We now need to compute the Reprojection Error between the points in
    %p_correct and the points we get by projecting grid corners to the
    %image using the new homography
    x1 = p_correct(:,1,i);
    y1 = p_correct(:,2,i);
    
    H = eval(['Hnew' num2str(i)]);
    points_projection = H*grid_coordinates';
    for j=1:length(points_projection)
        points_projection(:,j) = points_projection(:,j) /points_projection(3,j);
    end
    points_projection = points_projection';
    
    x2 = points_projection(:,1);
    y2 = points_projection(:,2);
    
    disp(["New Homography Reprojection error for >> " files(i)])
    total_err_reprojection = sum(sqrt((x1(:)-x2(:)).^2 + (x1(:)-x2(:)).^2));
    disp(["Total Reprojection Error (as Euclidean Distance) >> " total_err_reprojection]);
    disp(["Average Reprojection Error per point >> " total_err_reprojection/80]);
    
    H = eval(['H' num2str(i)]);
    points_projection_2 = H*grid_coordinates';
    for j=1:length(points_projection_2)
        points_projection_2(:,j) = points_projection_2(:,j) /points_projection_2(3,j);
    end
    points_projection_2 = points_projection_2';
    
    x2 = points_projection_2(:,1);
    y2 = points_projection_2(:,2);
    
    disp(["Part 2 Homography Reprojection error for >> " files(i)]);
    total_err_reprojection = sum(sqrt((x1(:)-x2(:)).^2 + (x1(:)-x2(:)).^2));
    disp(["Total Reprojection Error (as Euclidean Distance) >> " total_err_reprojection]);
    disp(["Average Reprojection Error per point >> " total_err_reprojection/80]);
end

    
