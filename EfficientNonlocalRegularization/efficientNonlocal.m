clc, clear,
addpath('..\flowColorCode')

im1 = imread('..\imagePairs\Army\frame10.png');
im2 = imread('..\imagePairs\Army\frame11.png');
%figure(1), clf
%imshow(im1);
if isinteger(im1);
    im1_8 = double(im1);
    im2_8 = double(im2);
end
%%
images  = cat(length(size(im1))+1, im1_8, im2_8);
size(images)
for i = 1:size(images,3)
    images(:,:,i,:)  = structure_texture_decomposition_rof( squeeze(images(:,:,i,:)));
end
%%
size(images)
tmp = images(:,:,:,1);
%imshow(uint8(tmp))
imshow(uint8(tmp))

%%
pyramid_spacing = 2;
pyramid_levels  =  1 + floor( log(min(size(images, 1), size(images,2))/16) / log(pyramid_spacing) );

% Construct image pyramid, using setting in Bruhn et al in  "Lucas/Kanade.." (IJCV2005') page 218

factor            = sqrt(2);  % sqrt(3)
smooth_sigma      = sqrt(pyramid_spacing)/factor;   % or sqrt(3) recommended by Manuel Werlberger   
f                 = fspecial('gaussian', 2*round(1.5*smooth_sigma) +1, smooth_sigma);    

pyramid_images    = compute_image_pyramid(images, f, pyramid_levels, 1/pyramid_spacing);


%%
centerpoint = [495, 195];% for Army
%centerpoint = [316, 53]; % for Teddy
%centerpoint = [476, 235]; % for Mequon
winSize = 90;
winRange = floor(centerpoint-winSize/2);
winRange = [winRange, winRange+winSize-1];
im1_w = im1(winRange(2):winRange(4),winRange(1):winRange(3),:);
%imshow(im1_w);

%%
im1_CLAB = rgb2lab(im1_w);
ci = im1_CLAB(45,45,:);
sigma_c_sq = 8^2;
sigma_x_sq = 15^2;
w = zeros(size(im1_CLAB,1),size(im1_CLAB,2));
for r = 1:size(im1_CLAB,1)
    for c = 1:size(im1_CLAB,2)
        cj = im1_CLAB(r,c,:);
        c_norm = (ci(2)-cj(2))^2+(ci(3)-cj(3))^2;
        p_norm = ((45-r)^2 + (45-c)^2);
%        p_norm = 0;
        w(r,c) = exp(-c_norm/(2*sigma_c_sq)-p_norm/(2*sigma_x_sq));
    end
end

surf(w);

%%
clc, clear
T = 50;
x = linspace(-200,200,10000);
ux = zeros(size(x));
%px = x.^2;                      % quadratic
%px = abs(x);                    % L1 norm
%px = (T/10)*log(1+(x.^2)/(2*1.5^2));   % Lorentzian
%px = sqrt((x.^2 + (0.001)^2));  % Charbonnier
%px = (x.^2 + (0.001)^2).^.45;   % Generalized Charbonnier
px = sqrt(abs(x))*5 ;           % something else

px = min(px,T);    %L1 norm
omega_n = [50 4];
sigma_n = [.5 .1];
data_n = [30 15 5;
          20 20 20];
      
K = size(data_n,2);
lb = zeros(size(data_n));
ub = ones( size(data_n)) * 30;
fun = @(data_n) (T - data_n(1,1) * exp(-(x.^2)/(2*data_n(2,1).^2)) ...
                   - data_n(1,2) * exp(-(x.^2)/(2*data_n(2,2).^2)) ...
                   - data_n(1,3) * exp(-(x.^2)/(2*data_n(2,3).^2)) ...
                   - px).*(0.1+exp(-(x.^2)/(2*(T/2).^2)));
fun2 = @(data_n) T - data_n(1,1) * exp(-(x.^2)/(2*data_n(2,1).^2)) ...
                   - data_n(1,2) * exp(-(x.^2)/(2*data_n(2,2).^2)) ...
                   - data_n(1,3) * exp(-(x.^2)/(2*data_n(2,3).^2));
options.Algorithm = 'levenberg-marquardt';
data_n = lsqnonlin(fun,data_n,[],[],options);
%data_n = lsqnonlin(fun,data_n,lb,ub);
data_n

plot(x,px,x,fun2(data_n))
axis([-100 100 0 T+10])




%imshow(flowToColor(uv));