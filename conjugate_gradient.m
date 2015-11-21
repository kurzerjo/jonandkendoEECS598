% conjugate gradient
function [x] = conjugate_gradient(A, b, x)

% using HDGF
Ax = A*x;

%initial residual vector
r = b - Ax;
p = r;

res_old = r'*r;

iteration = 0;

while iteration < 50
    
    % HDGF
    Ap =A*p;
    
    alpha = res_old/(p'*Ap);
    x = x + alpha*p;
    r = r - alpha*Ap;
    res_new = r'*r;
    
    if(sqrt(res_new) < 0.00001)
        break;
    end
    
    p = r + (res_new/res_old)*p;
    res_old = res_new;
    
    iteration = iteration + 1;
end


end