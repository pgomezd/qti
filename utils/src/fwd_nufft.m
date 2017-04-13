function y = fwd_nufft(x,par)
% Forward nufft operator from image space x, to k space y
% In
%     x: [Nx, Ny,(Nz) C, T] 2D coil imgs
% Out
%     y: [Nslices,Ncoils, ksamp, T] 

% adjust for coils
if par.ind.Ncoils == 1
    x = reshape(x,size(x,1),size(x,2),1,size(x,3));
end
% variables
dsamp = par.ind.dsamp;
nint = par.ind.nint;
Nslices = par.ind.Nslices;
Gn = par.recon.Gn;
NTimepoints = par.ind.NTimepoints;

if Nslices>1
    Ncoils = size(x,4);
    NTimepoints_x = size(x,5);
else
    Ncoils = size(x,3);
    NTimepoints_x = size(x,4);
end

if par.f.apply_tempsubspace || par.f.compress
    if par.f.compress
        V = par.V(:,1:par.ind.Tv);
    else
        V = par.V(:,1:par.ind.temp_coeff);
    end
end

y = zeros(Nslices,Ncoils,dsamp,NTimepoints,'single'); %review
if par.f.verbose_iter; fprintf('Applying forward nufft operator on coil data \n'); end;
for s=1:Nslices
    if par.f.verbose_iter; fprintf('Reconstructing slice %d/%d \n',s,Nslices); end;
    for c = 1:Ncoils
    cy = zeros(dsamp*nint,NTimepoints_x,'single');
    if par.f.verbose_iter; fprintf('Applying forward nufft operator to temporal frames of coil %d/%d \n',c,Ncoils); end;
        for t = 1:NTimepoints_x
            if Nslices>1
                cy(:,t) = (Gn * squeeze(x(:,:,s,c,t))); %forward model with nufft
            else
                cy(:,t) = (Gn * squeeze(x(:,:,c,t))); %forward model with nufft
            end 
        end

        %% project back to original space with tempsubpsace
        if par.f.apply_tempsubspace || par.f.compress
            if par.f.verbose_iter; fprintf('Projecting back to full temporal space \n'); end;
            cy = temporal_forward(cy,V);
        end

        %% select only sampled locations
        for t=1:NTimepoints
            ni  = mod((t-1),nint)+1; %find interleave index
            y(s,c,:,t) = cy((ni-1)*dsamp+1:ni*dsamp,t); %select only sampled location in dsamp*nint kspace grid
        end
    end
end
           
end