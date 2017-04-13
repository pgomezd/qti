function [out,par] = compress_data(data,par)
%% Compress data in temporal dimension
% =========================================================================
% [out,par] = compress_data(data,par)
%   Input
%       data         [N datadims, T timeoints] data              
%
%       par.        Anything that somehow controls the recon pipeline goes into this struct
%
% %
%   Output
%       out          [N datadims  Tv singular values] compressed data        
%
% %   
%   Function
%       Performs SVD compression of data in the temporal domain       
% =========================================================================
% v1.1: 31.07.16
% =========================================================================
% P. Gomez
% =========================================================================

%% Load/create V singular vectors
load(par.dir.dict,'V');
if exist('V','var')
    par.ind.NTimepoints = size(V,2); %#ok<NODEF>
else
    warning('Singular vectors V not found, compressing and saving dictionary');
    load(par.dir.dict,'D');
    if ~isfield(par.ind,'Tv')
        error('par.ind.Tv required for dicionary compression');
    end
    par.ind.NTimepoints = par.ind.Tv;
    [~,~,V]=svd(D,0); %#ok<NODEF>
    V=V(:,1:par.ind.NTimepoints);
    D=D*V;  %#ok<NASGU>
    save(par.dir.dict,'D','V','-v7.3');
end

%% Compress data
datadims = size(data);
data = reshape(data,[numel(data)/datadims(end) datadims(end)]); %concatenate all dims but the last one
data = data(:,:)*V;  % compress
out = reshape(data,[datadims(1:end-1) par.ind.NTimepoints]); %reshape back
end