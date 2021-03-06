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

%% Color consanctancy figure
im1 = imread('..\imagePairs\Army\frame10.png');
centerpoint = [495, 195];% for Army
winSize = 90;
winRange = floor(centerpoint-winSize/2);
winRange = [winRange, winRange+winSize-1];

im1_w = im1(winRange(2):winRange(4),winRange(1):winRange(3),:);


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
        w(r,c) = exp(-c_norm/(2*sigma_c_sq)-p_norm/(2*sigma_x_sq));
    end
end
w = uint8(w/max(w(:))*255);
w = cat(3, w, w, w);
im = horzcat(im1_w, w);
im = insertShape(im, 'circle', [45 45 5], 'LineWidth', 2,'Color','red');

imshow(im,'InitialMagnification','fit');

%% Penalty function estimation
clc, clear
T = 50;
x = linspace(-200,200,10000);
ux = zeros(size(x));
%px = x.^2;                      % quadratic
%px = abs(x);                    % L1 norm
%px = (T/10)*log(1+(x.^2)/(2*1.5^2));   % Lorentzian
%px = sqrt((x.^2 + (0.001)^2));  % Charbonnier
px = (x.^2 + (0.001)^2).^.45;   % Generalized Charbonnier
%px = sqrt(abs(x))*5 ;           % something else

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

h = plot(x,fun2(data_n),'r-.',x,px);
set(h(2),'LineWidth',2);
set(h(1),'LineWidth',4);
axis([-150 150 0 T+10])

