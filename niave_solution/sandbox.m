clc, clear,
addpath('.\flowColorCode')

im1 = imread('..\imagePairs\Army\frame10.png');
im2 = imread('..\imagePairs\Army\frame11.png');
figure(1), clf
%imshow(im1);

centerpoint = [495, 195];% for Army
%centerpoint = [316, 53]; % for Teddy
%centerpoint = [476, 235]; % for Mequon
winSize = 90;
winRange = floor(centerpoint-winSize/2);
winRange = [winRange, winRange+winSize-1];
im1 = im1(winRange(2):winRange(4),winRange(1):winRange(3),:);
imshow(im1);

%%
im1_CLAB = rgb2lab(im1);
ci = im1_CLAB(45,45,:);
sigma_c_sq = 15^2;
sigma_x_sq = 15^2;
w = zeros(size(im1_CLAB,1),size(im1_CLAB,2));
for r = 1:size(im1_CLAB,1)
    for c = 1:size(im1_CLAB,2)
        cj = im1_CLAB(r,c,:);
        c_norm = (ci(2)-cj(2))^2+(ci(3)-cj(3))^2;
        p_norm = ((45-r)^2 + (45-c)^2);
        w(r,c) = exp(-c_norm/(2*sigma_c_sq)-p_norm/(2*sigma_x_sq));
    end
end
imshow(w);

%%
im1 = imresize(im1,[100,100]);
im2 = imresize(im2,[100,100]);


alpha=1;
ite=100;
    
[u, v] = HS(im1, im2, alpha ,ite, [], [], 0);

uv(:,:,1) = u;
uv(:,:,2) = v;

imshow(flowToColor(uv));