function [out] = get_masked_voxels(data,par)
%% Get masked voxels
% =========================================================================
% Fit only voxels inside mask
% =========================================================================
% v1.1: 15.08.16
% v2.1: 30.03.17 - added mask adjustment 2x2 patches. Tag: pgd.2x2
% =========================================================================
% P. Gomez
% =========================================================================

%% get variables
if strcmp(par.dict.type,'atlas') %atlas-based dictionary pgd.atlas
    w_side = par.ind.w_side;
    S = par.ind.S;
    L = par.ind.L;
end

if ~isfield(data,'mask')
   % --- 30.03.2017 pgd.2x2 --- %
   data.mask = ones(par.ind.N,1)>0; %if no mask, fit all voxels with boolean mask
end
%% mask data
if strcmp(par.dict.type,'cluster') %cluster-based dictionary
    d.x = single(data.Xp(data.mask>0,:));
    % --- 31.07.2016 pgd.atomind --- %
    par.atomind = [];
    % --- 31.07.2016 pgd.atomind --- %
% --- pgd.atlas --- %
elseif strcmp(par.dict.type,'atlas') %atlas-based dictionary pgd.atlas
    d.x = single(data.Xp(data.mask>0,:));
    if isfield(data,'qmap')
        d.q = single(data.Qp(data.mask>0,:));
    end
    Nd = size(d.x,1); 
    if par.f.windowpatches % windowpatches == get patches for the entire size of the search window
        atomind = zeros(L,Nd); %L x Nd
        mask_ind = single(1:numel(data.mask)*S);
        mask_ind = bsxfun(@times,reshape(mask_ind,[size(data.mask),S]),data.mask);

        % Get mask_ind window patches
        if length(datadims)==3 %2D + Q data
            Xmask = padarray(mask_ind,[w_dist(1) w_dist(2)]); %zero-padding
            Xmask_p = zeros(N,W,S);
            ws1 = w_side(1);
            ws2 = w_side(2);
            if par.f.use_parallel
                parfor s = 1:S
                    Xmask_c =  Xmask(:,:,s);
                    Xmask_p(:,:,s) = im2col(Xmask_c,[ws1 ws2])'; %mask
                end
            else
                for s = 1:S
                    Xmask_c =  Xmask(:,:,s);
                    Xmask_p(:,:,s) = im2col(Xmask_c,[ws1 ws2])'; %mask
                end
            end
            clear Xmask_c;
        elseif length(datadims)==4 %3D + Q data
            Xmask_p = vol2col(mask_ind,[w_side(1) w_side(2) w_side(3)]); %slow - should be improved
        end

        Xmask_p = Xmask_p(:,:);
        Xmask_p = Xmask_p(data.mask>0,:)'; 

        %Select indexes inside atlas mask (full ST + atlas) and subject mask
        for n = 1:size(d.x,1)
            atomind(par.amask(Xmask_p(Xmask_p(:,n)>0,n)),n)=1; %slow - should be improved
        end
        clear Xmask_p;
        clear mask_ind;
        par.atomind = atomind;
        clear atomind;
    else %no window patches == iterate on every voxel
        if par.f.verbose; fprintf('Voxel-wise indexing \n'); end;
        par.atomind = get_atomind(data.mask,par.amask,w_dist,L,Nd,datadims,S);            
    end
else
    error('unknown dictionary type, par.dict.type should be either cluster or atlas');
end
% --- pgd.atlas --- % 

%% out
out = d;
out.mask = data.mask;
out.X = data.X;
out.Xp = data.Xp;
if par.f.Yout && isfield(data,'Y')  
   out.Y = data.Y; 
end
end