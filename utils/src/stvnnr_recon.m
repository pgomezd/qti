function [out,par] = stvnnr_recon(data,par)
% % Fast Spatio-temporal TV + Global LR reconstruction based on primal-dual splitting 
% =========================================================================
% [out,par] = fb_recon(data,par)
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
%       Spatio-temporal TV + Global LR minimiation based reconstruction
% =========================================================================
% v1.1: 14.07.17
% =========================================================================
% Code adapted by C. Ulas
% =========================================================================
 
%% get operators
A_for = par.oper.A_for;
A_adj = par.oper.A_adj;

if par.f.match
    D_match = par.oper.D_match;
end

% load finite-difference matrices for TV
load dif_Q.mat;
D.Dx = Q1;
D.Dy = Q2;
D.Dt = Q3;

param = par.recon.STVNNR;
param.A_for = A_for;
param.A_adj = A_adj;
param.D = D;

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
par.recon.history = stvnnr_mrf_recon(alpha_ref, data.X, data.Y, param);

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