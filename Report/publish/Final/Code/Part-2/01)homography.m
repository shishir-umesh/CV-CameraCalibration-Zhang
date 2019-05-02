%Obtain the world coordinates of the chessboard pattern
%Width and height square size of chessboard pattern as specified in assignment
width = 9;
height = 7;
squareSize = 30;

%Get the grid corners from the chessboard pattern from the world
%coordinates system
%We start from the origin(bottom left corner) and move clockwise
world_coordinates = [0,0,1; 0,height*squareSize,1; width*squareSize,height*squareSize,1;width*squareSize,0,1];

% To convert from Nx3 into a 3xN matrix
world_coordinates=world_coordinates';
disp(world_coordinates);

%Read the 4 images provided of the grid
img1 = imread("images2.png");
img2 = imread("images9.png");
img3 = imread("images12.png");
img4 = imread("images20.png");

%Manually accept the grid corners from the user
imshow(img1), [x1,y1] = ginput(4);
imshow(img2), [x2,y2] = ginput(4);
imshow(img3), [x3,y3] = ginput(4);
imshow(img4), [x4,y4] = ginput(4);

%Display the grid corners 
disp([x1,y1]);
disp([x2,y2]);
disp([x3,y3]);
disp([x4,y4]);

ones = [1,1,1,1];

%Obtain the 3xN form of the image coordinates --> [X,Y,1]
image_coordinates_1 = [x1(:), y1(:), ones(:)]';
image_coordinates_2 = [x2(:), y2(:), ones(:)]';
image_coordinates_3 = [x3(:), y3(:), ones(:)]';
image_coordinates_4 = [x4(:), y4(:), ones(:)]';

%Compute the homography2d(x1, x2) where x1-->x2
H1 = homography2d(world_coordinates, image_coordinates_1);
H2 = homography2d(world_coordinates, image_coordinates_2);
H3 = homography2d(world_coordinates, image_coordinates_3);
H4 = homography2d(world_coordinates, image_coordinates_4);

%Display the computed homographies
disp("Homography of --> images2.png")
disp(H1)
disp("Homography of --> images9.png")
disp(H2)
disp("Homography of --> images12.png")
disp(H3)
disp("Homography of --> images20.png")
disp(H4)