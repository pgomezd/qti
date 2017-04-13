%% col2imn: columns to 2D image 
% =========================================================================
% out = col2imn(col,im_dim,p_side)
% %
%   Input
%       
%      col            [Nx*Ny, Pi*Pj, (Q)] 2D image patches as column vectors per voxel in the image 
%      im_dim         [Nx,Ny]             image dimensions
%      patch_side     [Pi]                patch side
%
%   Output
%      out            [Nx,Ny,(Q)]        2D image volume (+ Q parameters - optional)
%
%   Function
%         This function undoes im2col.m by creating a 2D + Q volume. 
%               - Only dense patch extraction implemented.
%               - Uses zero-padding to treat boundary conditions
%               - Updates the volume by averaging all contributions of
%               every patch
% =========================================================================
% v1.1: 31.10.15
% =========================================================================
% P. Gomez
% GE Global Research
% =========================================================================
function out = col2imn(col,im_dim,patch_side)
%% setup and defaults
if numel(im_dim)~=2
    error('Only 2D data dimensions supported');
end

if numel(patch_side)~=1
    error('Only patches with same side supported');
end
patch_dist = floor(patch_side/2); %distance to center patch
P = patch_side^numel(im_dim);
col_dim = size(col);
if numel(col_dim)>2
   col = col(:,:,:); %concatenate last dimension
   Q = size(col,3); %get number of parameters
else 
    Q = 1; 
end

%% Reshape and shift
col = reshape(col,[im_dim,P,Q]);
col = padarray(col,[patch_dist patch_dist]);

for p = 1:P 
        col(:,:,p,:) = circshift(col(:,:,p,:),p-(ceil(p/patch_side)-1)*patch_side-patch_side+patch_dist,1); %shift in x 
        col(:,:,p,:) = circshift(col(:,:,p,:),ceil(p/(patch_side))-patch_side+patch_dist,2); %shift in y 
end

%Update considering boundary conditions
imdim = numel(im_dim);
p=1:P;
ind_patch_x = (p-(ceil(p/patch_side)-1)*patch_side-patch_side+patch_dist);
ind_patch_y = ceil(p/(patch_side))-patch_side+patch_dist;

out = zeros([size(col,1),size(col,2),Q]);

%inds {dims,conds}
inds = cell(imdim,3);
for d=1:imdim
    inds{d,1} = patch_dist+1; %-1
    inds{d,2} = patch_dist+1+1:size(col,d)-patch_dist-1; %0
    inds{d,3} = size(col,d)-patch_dist; %1
end

voldim = 3; %smallest neihgboring dimensions
for b = 1:imdim^3  %3 conditions per dimension -1,0,1
    
    cond_x = b-(ceil(b/voldim)-1)*voldim-voldim+1;
    cond_y = ceil(b/(voldim))-voldim+1;
    
    ind_x = inds{1,cond_x+2};
    ind_y = inds{2,cond_y+2};
    
    p_x = [find(ind_patch_x == cond_x),find(ind_patch_x == 0)]; %find indexes that match boundary or are not 0
    p_y = [find(ind_patch_y == cond_y),find(ind_patch_y == 0)]; 
    
    ind_patch = sort(p_x(ismember(p_x,p_y)));
    
    out(ind_x,ind_y,:) = ...
        squeeze(mean(col(ind_x,ind_y,ind_patch,:),3)); %mean == patch-wise update
end
out=out(patch_dist+1:end-patch_dist,patch_dist+1:end-patch_dist,:);
end
