function [out,par] = admm_recon(data,par)
% % Temporal subpsace reconstruction based on admm 
% =========================================================================
% [out,par] = cs_recon(data,par)
%   Input
%       data.Y          [Ncoils, Full K space, NTimepoints]  Viewshared/zero-filled raw data             
%       data.Y          [Nx,Ny,(Nz), T timepoints]           Undersampled data              
%       data.X          [Nx,Ny,(Nz), T timepoints]           Initial image estimate
%  
%       par.            Anything that somehow controls the recon pipeline goes into this struct
% 
%   Output
%       out.X           [Nx,Ny,(Nz),  T timepoints]          2D+T Reconstructed image
%       out.qmap        [Nx,Ny,(Nz),  Q maps]                2D+Q Estimated parametric maps
%   Function
%       Temporal subspace reconstruction based on Tamir et al.
% =========================================================================
% v1.1: 19.09.16
% =========================================================================
% ADMM code from J. Tamir, adapted by P. Gomez
% =========================================================================
 
%% get operators
A_for = par.oper.A_for;
AHA = par.oper.AHA;
if par.f.match
    D_match = par.oper.D_match;
end

%% ADMM
iter_ops.max_iter = par.recon.admm_max_iter;
iter_ops.rho = par.recon.admm_rho;
iter_ops.objfun = @(a, sv, lam) 0.5*norm_mat(data.Y - A_for(a))^2 + lam*sum(sv(:));

llr_ops.lambda = par.recon.admm_lambda;
llr_ops.block_dim = [par.recon.llr_block_dim par.recon.llr_block_dim];

lsqr_ops.max_iter = par.recon.lsqr_max_iter;
lsqr_ops.tol = par.recon.lsqr_tol;

alpha_ref = RefValue;
if par.f.apply_tempsubspace
    if par.ind.Nslices > 1
        alpha_ref.data = zeros(par.ind.ix, par.ind.iy, par.ind.iz, par.ind.temp_coeff);
    else
        alpha_ref.data = zeros(par.ind.ix, par.ind.iy, par.ind.temp_coeff);
    end
else
    if par.ind.Nslices > 1
        alpha_ref.data = zeros(par.ind.ix, par.ind.iy, par.ind.iz, par.ind.NTimepoints);
    else
        alpha_ref.data = zeros(par.ind.ix, par.ind.iy, par.ind.NTimepoints);
    end
    
end
par.recon.history = iter_admm(alpha_ref, iter_ops, llr_ops, lsqr_ops, AHA, data.X, @admm_callback);

%% Match to dictionary + out
if par.f.compress
    V = par.V(:,1:par.ind.Tv);
    data.X = temporal_forward(alpha_ref.data,V); %forward temporal on alpha vals
elseif par.f.apply_tempsubspace
    V = par.V(:,1:par.ind.temp_coeff);
    data.X = temporal_forward(alpha_ref.data,V); %forward temporal on alpha vals
end
out = D_match(data); %match to dictionary to get qmaps 
if par.f.alphaout
    out.alpha = alpha_ref.data; %can be used to calculate params in projected space
end

end
