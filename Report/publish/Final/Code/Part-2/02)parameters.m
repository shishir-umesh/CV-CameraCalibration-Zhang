%We follow the instructions in section 2.4.4 of the Camera Calibration
%chapter by Zhengyou Zhang

%Since we have 4 computed homographies, we can stack n equations to get Vb=0
V = [];
for i = 1:4
    %Assigns the i(th) Homography to the variable H
    H = eval(['H' num2str(i)]);
    
    %Let the i-th column vector of H(homography) be h(i) =
    %[h(i1),h(i2),h(i3)]'
    h1 = H(:,1);
    h2 = H(:,2);
    h3 = H(:,3);
    
    %Using v(ij) = [hi(1)*hj(1), hi(1)*hj(2)+hi(2)*hj(1), hi(2)*hj(2), hi(3)*hj(1)+hi(1)*hj(3), hi(3)*hj(2)+hi(2)*hj(3),
    %hi(3)*hj(3)]'
    v11 = [h1(1)*h1(1), h1(1)*h1(2)+h1(2)*h1(1), h1(2)*h1(2), h1(3)*h1(1)+h1(1)*h1(3), h1(3)*h1(2)+h1(2)*h1(3), h1(3)*h1(3)]';
    v12 = [h1(1)*h2(1), h1(1)*h2(2)+h1(2)*h2(1), h1(2)*h2(2), h1(3)*h2(1)+h1(1)*h2(3), h1(3)*h2(2)+h1(2)*h2(3), h1(3)*h2(3)]';
    v22 = [h2(1)*h2(1), h2(1)*h2(2)+h2(2)*h2(1), h2(2)*h2(2), h2(3)*h2(1)+h2(1)*h2(3), h2(3)*h2(2)+h2(2)*h2(3), h2(3)*h2(3)]';
    
    %We stack the 'n' equations for 'n' images to obtain V
    V = [V; v12'; (v11-v22)'];
end

%We now need to solve for Vb = 0
%The solution for Vb = 0 from the chapter is said to be well known as the
%eigenvector of V'.V associated with the smallest eigenvalue(equivalently,
%the right singular vector of V associated with the smallest singular
%value).
%We acheive this by by applying the singular value decomposition on V.
[U, Sigma, V_transpose] = svd(V);

b = V_transpose(:,end);

%Once 'b' is estimated, we can now compute all the camera intrinsic
%parameters by first describing the matrix "B".

%We know that B is symmetric and is defined by a 6D vector as follows,
% b = [B11, B12, B22, B13, B23, B33]'
% Where,
% B = [ B11,    B12,    B13]
%     [ B12,    B22,    B23]
%     [ B13,    B23,    B33]
B11 = b(1);
B12 = b(2);
B22 = b(3);
B13 = b(4);
B23 = b(5);
B33 = b(6);

B = [B11, B12, B13; B12, B22, B23; B13, B23, B33];

%Display the computed matrix B
disp("Matrix B >>");
disp(B);

%We can now calculate the intrinsic parameters using the equations
%mentioned in Page 21 of Camera Calibration chapter.
v0 = (B12*B13 - B11*B23)/(B11*B22 - B12^2);
lambda = B33 - (B13^2 + v0*(B12*B13-B11*B23))/B11;
alpha = sqrt(lambda/B11);
beta = sqrt(lambda*B11/(B11*B22-B12^2));
gamma = -B12*alpha^2*beta/lambda;
u0 = gamma*v0/alpha - B13*alpha^2/lambda;

%Therefore our intrinsic matrix A can be defined as follows,
A = [alpha, gamma, u0; 0, beta, v0; 0, 0, 1];

%Display the computed intrinsic parameters matrix K
disp("Intrinsic Parameters matrix >>");
disp(A);

%Now once we have computed the matrix B and the intrinsic parameters of the 
%camera.We can now compute the R and t values for each image.
%We follow the steps mentioned in page 21.

files = ["images2", "images9", "images12", "images20"];

for i = 1:4
    %Assigns the i(th) Homography to the variable H
    H = eval(['H' num2str(i)]);
    
    %Let the i-th column vector of H(homography) be h(i) =
    %[h(i1),h(i2),h(i3)]'
    h1 = H(:,1);
    h2 = H(:,2);
    h3 = H(:,3);
    
    %INV(A)*b can be slower and lass accurate than A\b. Consider using A\b for 
    %INV(A)*b or b\A for b*INV(A).
    lambda_r = 1/ norm(A\h1);
    r1 = lambda_r*(A\h1);
    r2 = lambda_r*(A\h2);
    r3 = cross(r1,r2);
    t  = lambda_r*(A\h3);
    
    R = [r1, r2, r3];
    
    disp(["Rotation Matrix for >> " files(i)]);
    disp(R);
    disp(["Translation vector for >> " files(i)]);
    disp(t);
    
    %Now we need to confirm that our rotation matrix is in fact a rotation matrix.
    %We print R'*R to check if it is an identity matrix.
    disp(["Transpose(R)*R for >> " files(i)]);
    disp(R'*R);
    
    %We notice that it is not an identity matrix, hence we can enforce R to be
    %a rotation matrix by SVD decomposition of R by setting the singular values to ones
    % U*Sigma*V_T = svd(R), new rotation matrix will be U*V_T
    [U, Sigma, V_transpose] = svd(R);
    R_new = U*V_transpose;
    
    disp(["New rotation matrix for >> " files(i)]);
    disp(R_new);
    disp(["Transpose(R_new)*R_new for >> " files(i)]);
    disp(R_new'*R_new);
    
end



