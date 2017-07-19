function [ history ] = stvnnr_mrf_recon(x_ref, im_init, y_zf,  param, gt_data)
%iter_admm ADMM algorithm for solving locally low rank regularization:
%
% \min_x 0.5 * || Ax - y ||_2^2 + lambda_1 * |x|_{TV} + lambda_2 * |x|_{*} ) (1)
%   where |x|_{TV} denotes spatial 2D or spati-temporal 3D total variation
%   and |x|_{*} denotes the global low-rank term.
%
% Inputs:
%  x_ref = pointer to solution to (1), with result stored in x_ref.data 
%  in_init = initial estimate of the the recon images
%  y_zf =  k-space measurements
%  param = struct keeping the all the parameters and operations used in the algorithm 

% Optional:
%  gt_data = fully sampled image data to be used for comparison in
%  retropspective sampling
%
% Outputs:
%  history = struct of history/statistics from the iterations of the algorithm 

if nargin < 5
    use_gt = false;
else
    use_gt = true;
end

%% Operators
A_for = param.A_for;
A_adj = param.A_adj;

%% Parameters and finite-difference matrices
lambda_1 = param.lambda1;
lambda_2 = param.lambda2;
L = param.L;
maxIter = param.IterNo;
D = param.D;

%% Initialization
[Nx,Ny,Nt] = size(im_init);
im_init = double(im_init);
X_est = im_init;
Y_est = im_init;
t1 = 1; t2 =1;

lambda_shrink = (t1*lambda_2)./(1+t1*L);

psnr_res = zeros(1,maxIter);
obj_fun = zeros(1,maxIter);
fprintf('\n STVNNR reconstruction is starting..!');

tic
for k = 1:maxIter
    
    fprintf('\n Iteration %i',k);
    
    %% Compute gradient
    f1 = A_for(X_est);
    fg = (t1/(1+t1*L))*double(A_adj(f1-y_zf));
    ft = ((t1*lambda_1)/(1+t1*L))*compute_TV(Y_est,D,'b');
    X_p = X_est-fg-ft;
  
    X_p = reshape(X_p,[Nx*Ny, Nt]);
    
    %% Evaluate matrix shrinkage operator (X subproblem)
    X_new = SVD_shrinkage(X_p, lambda_shrink);
    
    X_new = reshape(X_new, Nx, Ny, Nt);
    %% Compute Y subproblem
    Y_p = Y_est + t2*lambda_1*compute_TV((2*X_new-X_est), D, 'f');
    
    %% Project Y onto l-inf unit ball
    Y_new = sign(Y_p).*min(abs(Y_p),1);
    
    %% Assign the previous values
    Y_est = Y_new;
    X_est = X_new;
    
    %% Cost calculation
    X_orig = X_new;
    
    e = A_for(X_orig) - y_zf;     
    [~,sigmaq,~] = givefastSVD(reshape((X_orig),[Nx*Ny,Nt]));
    V1 = sum(sum(sum(compute_TV(X_orig, D, 'f'))));
    obj_fun(k) = sum(abs(e(:)).^2) + lambda_1*V1 + lambda_2*sum(abs(sigmaq(:)));
    
    if use_gt
        psnr_val=PSNR(abs(X_orig),abs(gt_data));
        fprintf(': Cost function = %.3f, PSNR = %.3f dB...',obj_fun(k), psnr_val);
        psnr_res(k) = psnr_val ;
        
        if (k>1) && (psnr_res(k)-psnr_res(k-1)<0)
            break ;
        end
    else
        if (k>1) && (obj_fun(k)-obj_fun(k-1)>0)
            break ;
        end
        
         fprintf(': Cost function = %.3f..', obj_fun(k));
    end
    
    x_ref.data = X_orig;
end

t2 = toc;

history.nitr = k;
history.run_time = t2;
history.psnr = psnr_res;
history.obj_func = obj_fun;


end


function [tv_sum] = compute_TV(input, D, direction)
% Inputs:
%    input = 3D input matrix in the form of Nx X Ny X Nt
%    D = A struct keeping the matrices to calculate finite-differences
%        along X,Y and time dimension (if it is for 3D)
%    direction = 'f' for forward TV, 'b' is backward (transpose) TV
%                 operation
  
% Output:
%    tv_sum = sum (L1-norm) of TV term along all available dimensions


if ~isfield(D,'Dt')
    type = '2D';
else
    type = '3D';
end

[Nx,Ny,Nt] = size(input);
Dx=D.Dx; Dy=D.Dy; Dt=D.Dt;

if strcmp(type,'2D')
    input = reshape(input,[Nx*Ny,Nt]);
    ss_tv = zeros(size(input));

    if direction == 'f'
        for n=1:size(input,2)
            D1 = abs(Dx*input(:,n))+ eps;
            D2 = abs(Dy*input(:,n))+ eps;
            ss_tv(:,n) = D1+D2;
        end
        
    elseif direction == 'b'
        for n=1:size(input,2)
            D1 = abs(Dx'*input(:,n))+ eps;
            D2 = abs(Dy'*input(:,n))+ eps;
            ss_tv(:,n) = D1+D2;
        end
    end
    
elseif strcmp(type,'3D')
    if direction == 'f'
        D1 = abs(Dx*input(:))+ eps;
        D2 = abs(Dy*input(:))+ eps;
        D3 = abs(Dt*input(:))+ eps;
        ss_tv =  D1+D2+D3;
        
    elseif direction == 'b'
        D1 = abs(Dx'*input(:))+ eps;
        D2 = abs(Dy'*input(:))+ eps;
        D3 = abs(Dt'*input(:))+ eps;
        ss_tv =  D1+D2+D3;    
    end
end
  
tv_sum = reshape(ss_tv, [Nx Ny Nt]);

end



function  [data_shrinked] = SVD_shrinkage (data, lambda)
% Inputs:
%    data = 2D input matrix in the form of (Nx X Ny) X Nt
%    lambda = scalar value denoting the threshold for singular values
% Output:
%    data_shringked = output data obtained by shrinking the input data

[u,sigma,v] = givefastSVD(data);
s=diag(sigma);

thres = lambda;

s=(s-thres);
s = s.*(s>0);
data_shrinked=u*(diag(s))*v';

end


