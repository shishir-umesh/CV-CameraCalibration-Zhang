n1 = imread("images2.png");
imshow(n1)
[x1,y1] = ginput(4);
hold on;
plot(x1, y1, 'o', 'LineWidth', 2, 'MarkerSize', 15);
z1 = [1,1,1,1];
h1 = [x1(:), y1(:), z1(:)];
h1 = transpose(h1);


n2 = imread("images9.png");
imshow(n2)
[x2,y2] = ginput(4);
hold on;
plot(x2, y2, 'o', 'LineWidth', 2, 'MarkerSize', 15);
z2 = [1,1,1,1];
h2 = [x2(:), y2(:), z2(:)];
h2 = transpose(h2);

H = homography2d(h1, h2)
