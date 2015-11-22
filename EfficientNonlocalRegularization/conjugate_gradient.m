% conjugate gradient
function [x] = conjugate_gradient(images, uv, B, wB, it)
nl_weight = 1.5;
npixels = size(wB,1)/2;
jointImg = images(:,:,:,1);

sigma_s = 9;   %sigma_x
sigma_r = 8;    %sigma_c

sz = [size(uv,1) size(uv,2)];

[term1,~] = adaptive_manifold_filter(ones(sz), sigma_s, sigma_r, [],jointImg);

% compute wA and b
wA = zeros(npixels*2,1);

for i = 1:2
    [term2,~] = adaptive_manifold_filter(uv(:,:,i), sigma_s, sigma_r, [],jointImg);
    tmp = term1.*uv(:,:,i) - term2;
    wA(npixels*(i-1)+1:npixels*i) = tmp(:);
end



b = nl_weight*wA + wB;
% initial guess
x = zeros(2*npixels,1);
Ax = x;
%initial residual vector
r = b - Ax;
p = r;

res_old = r'*r;
iteration = 0;

while iteration < it
    
    % HDGF
    for i = 1:2
        half_p = reshape(p(npixels*(i-1)+1:npixels*i), sz);
        [term2,~] = adaptive_manifold_filter(half_p, sigma_s, sigma_r, [],jointImg);
        
        tmp = term1.*half_p - term2;
        Ap(npixels*(i-1)+1:npixels*i) = tmp(:);
    end
    Ap = reshape(Ap,length(Ap(:)),1);
    Ap = nl_weight*Ap + B*p;
    
    alpha = res_old/(p'*Ap);
    x = x + alpha*p;
    r = r - alpha*Ap;
    res_new = r'*r;
    
    if(norm(r) < 1e-6)
        break;
    end
    
    p = r + (res_new/res_old)*p;
    res_old = res_new;
    
    iteration = iteration + 1;
end


end