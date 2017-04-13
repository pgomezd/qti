function [x,par] = adj_nufft(y,par)
% Adjunt nufft operator from k space y to image space x
% In
%     y: [Nz,Ncoils, ksamp, T] sampled k-space data
% Out
%     x: [Nx, Ny, (Nz),C, T] 2D coil imgs 
%   par.V data-driven SVD

% variables
Nslices = size(y,1);
Ncoils  = size(y,2);
if par.f.apply_tempsubspace
    NTimepoints = par.ind.temp_coeff;
else
    NTimepoints = par.ind.NTimepoints;
end
mtx_reco = par.ind.mtx_reco;
Gn = par.recon.Gn; %nufft operator

x = zeros(mtx_reco,mtx_reco,Nslices,Ncoils,NTimepoints,'single');
if par.f.verbose_iter; fprintf('Reconstructing image with nufft operator \n'); end;
for s=1:Nslices
    if par.f.verbose_iter; fprintf('Reconstructing slice %d/%d \n',s,Nslices); end;
    craw = get_full_k(squeeze(y(s,:,:,:)),par); %do viesharing/zero-filling + compression    
    %% apply tempsubspace
    if par.f.apply_tempsubspace || par.f.compress
        if par.f.verbose_iter; fprintf('Applying temporal subspace compression on data \n'); end;
        if par.f.compress_with_data && ~isfield(par,'V')
            fprintf('Computing subspace from single data slice \n');
            tmp=permute(y,[4 1 2 3]);
            [~,~,par.V]=svd(tmp(:,:).',0); 
            par.V=par.V(:,1:max(par.ind.Tv,par.ind.temp_coeff));
        end
            
        if par.f.compress 
            V = par.V(:,1:par.ind.Tv);
        else
            V = par.V(:,1:par.ind.temp_coeff);
        end
        craw = temporal_adjoint(craw,V);
    end
    for c = 1:Ncoils
        if par.f.verbose_iter; fprintf('Reconstructing temporal frames of coil %d/%d \n',c,Ncoils); end;
        for t = 1:NTimepoints %recon each image individually
            x(:,:,s,c,t) = (Gn'*squeeze(craw(c,:,t))); %recon with nufft, dcf previously applied
        end
    end
end
x = squeeze(x);