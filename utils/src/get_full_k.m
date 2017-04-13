function [out,par] = get_full_k(raw,par)
% % get full k space data 
% =========================================================================
% [out] = get_full_k(data,par)
%   Input
%       raw       [Ncoils, K sampled points, NTimepoints]  raw data - single slice               
%
%       par.        Anything that somehow controls the recon pipeline goes into this struct
% 
%   Output
%       out       [Ncoils, Full K space, NTimepoints]  vieshared/zero-filled raw data 
%   Function
%       1. Zero-fills or does viewsharing on raw data to get full kspace
%       2. Calculates sampling mask
%       3. Compresses data into V singular values (if par.f.compress)
% =========================================================================
% v1.1: 15.09.16
% =========================================================================
% (c) P. Gomez
% =========================================================================
%% Adjust for single coil
if ismatrix(raw)
    raw = reshape(raw,[1 size(raw)]);
end
%% Get variables
Ncoils = size(raw,1);
dsamp = size(raw,2);
NTimepoints = size(raw,3);
nint = par.ind.nint;

%% raw method: zero-filling/viewsharing/collapsing into single k
switch par.recon.raw_method
    case 'zerofill'
        craw = zeros(Ncoils,dsamp*nint,NTimepoints,'single'); %entire k-space as zero
        for nc = 1:nint %loop over interleaves - each interleave has a unique k-space location
            crepind = nc:nint:NTimepoints; %index interleaves with the same k-space location
            ni=mod((nc-1),nint)+1; %find k-space location in dsamp*nintl matrix
            craw(:,(ni-1)*dsamp+1:ni*dsamp,crepind) = raw(:,:,crepind);
        end
    case 'viewshare'
    case 'collapse'
end

%% out
out = craw;

end