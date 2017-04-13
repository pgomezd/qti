function [out,par] = get_sens(data,par)
% % get coil kernel and estimate of data 
% =========================================================================
% [out,par] = get_sens(raw,par);
%   Input
%       Either:
%           data.Y      [Nz,Ncoils, ksamp, T timepoints]  raw data             
%       Or:
%           data.Y      [SpokeLength,Nspokes, C coils, T timepoints] raw data
%
%       par.            Anything that somehow controls the recon pipeline goes into this struct
% 
%   Output
%       out.X           [Nx,Ny,(Nz),  T timepoints]          2/3D+T Reconstructed image
%       par.recon.sens  [Nx,Ny,(Nz),  C coils]               Coil kernel
%   Function
%       1. Reconstructs coil images with nufft/gridding operator
%       2. Estimates coil kernel with function from G. Buonincontri
% =========================================================================
% v1.1: 15.09.16
% =========================================================================
% (c) P. Gomez
% =========================================================================

if par.f.verbose; fprintf('Getting initial image estimate \n'); end;
%% apply dcf on raw data
switch par.recon.psd
    case 'qti'
        switch par.recon.sp_method
            case 'grid'
            case 'grid_cart'
            case 'nufft'
                raw_in = bsxfun(@times,data.Y,par.recon.dcf); 
            case 'nufft_cart'
        end      
    case 'rufis'
    case 'cart'
    case 'sim'
end


%% recon coil image
switch par.recon.psd
    case 'qti'
        switch par.recon.sp_method
            case 'grid'
            case 'grid_cart'
            case 'nufft'
                [out_coils,par] = adj_nufft(raw_in,par); 
            case 'nufft_cart'
        end
    case 'rufis'
    case 'cart'
    case 'sim'
end


%% get coil kernel and estimate image
switch par.recon.psd
    case 'qti'
        %% No parallel imaging in this version, simply output the coils
        out.X = out_coils;
    case 'rufis'    
end

%% out
if par.f.apply_tempsubspace && ~(strcmp(par.recon.method,'ADMM') || strcmp(par.recon.method,'IPA'))
    out.X = temporal_forward(out.X,par.V(:,1:par.ind.temp_coeff)); %project back to full temp space for dictionary matching
end
out.Y = data.Y; %output raw k-space as well (with no temporal projection or density compensation)

%% Use data to get mask
if par.f.calcmask
    out.mask = get_mask(out.X,par); %review mask thresh
end

end