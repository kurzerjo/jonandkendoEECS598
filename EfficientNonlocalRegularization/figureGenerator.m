clc, clear,
addpath('..\flowColorCode')

im1 = imread('..\imagePairs\Mequon\frame10.png');
im2 = imread('..\imagePairs\Mequon\frame11.png');
if isinteger(im1);
    im1 = double(im1);
    im2 = double(im2);
end
%% Create Structure Texture Decomp image 
im1_tex = structure_texture_decomposition_rof(im1);

%% Show an original image and the structure decomp visions in same frame
quHei = floor(size(im1_tex,1)/4);
tmp = vertcat(im1(quHei:3*quHei,:,:),im1_tex(quHei:3*quHei,:,:));
tmp([1,2,end/2+1,end/2,end-1,end],:,:) = 0;
tmp(:,[1,2,end-1,end],:) = 0;
imshow(uint8(tmp))

%% Show pyramids of images in decreasing coarseness
clc
pyramid_spacing = 2;
pyramid_levels  =  1 + floor( log(min(size(im1_tex, 1), size(im1_tex,2))/16) / log(pyramid_spacing) );

% Construct image pyramid, using setting in Bruhn et al in  "Lucas/Kanade.." (IJCV2005') page 218
factor            = sqrt(2);  % sqrt(3)
smooth_sigma      = sqrt(pyramid_spacing)/factor;   % or sqrt(3) recommended by Manuel Werlberger   
f                 = fspecial('gaussian', 2*round(1.5*smooth_sigma) +1, smooth_sigma); 
pyramid_images    = compute_image_pyramid(im1_tex, f, pyramid_levels, 1/pyramid_spacing);

tmp = pyramid_images{1};
for i = 2:length(pyramid_images)
    tmp2 = pyramid_images{i};
    tmp2(end,:,:) = 0; tmp2(:,end,:) = 0;
    tmp(1:size(tmp2,1), 1:size(tmp2,2),:) = tmp2;
end
tmp([end-1,end],:,:) = 0; tmp(:,[end-1,end],:) = 0;
tmp([1,2],:,:) = 0; tmp(:,[1,2],:) = 0;
imshow(uint8(tmp))


