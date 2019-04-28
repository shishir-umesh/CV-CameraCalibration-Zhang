%We follow the instructions in section 2.4.4 of the Camera Calibration
%chapter by Zhengyou Zhang

%Since we have 4 computed homographies, we can stack n equations to get Vb=0
V = [];
for i = 1:4
    %Assigns the i(th) Homography to the variable H
    H = eval(['H' num2str(i)]);
    
    %Let the i-th coloumn vector of H(homography) be h(i) =
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

%Once 'b' is estimated, we can now compute all the camera instrinsic
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

%We also know that B is defined as B = (lambda)*(inverse(A))'*A
