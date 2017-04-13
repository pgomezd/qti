%% vol2col: 3D volume to columns 
% =========================================================================
% out = vol2col(vol,patch)
% %
%   Input
%         
%       vol             [Nx,Ny,Nz,(P)]  3D image volume (+ Q parameters - optional)
%       patch_dim       [Pi,Pj,Pk]  3D patch size
%
%   Output
%       out             [Nx*Ny*Nz, Pi*Pj*Pk, (Q)] 3D image patches as row vectors per voxel in the image 
%
%   Function

%         Extracts 3D image patches from a 3D volume.
%               - Only dense patch extraction implemented.
%               - Uses zero-padding to treat boundary conditions
% =========================================================================
% v1.1: 29.10.15
% v1.2: 31.10.15 - extension to N-D datasets
% =========================================================================
% P. Gomez
% GE Global Research
% =========================================================================
%%
function out = vol2col(vol,patch_dim)
%% setup and defaults
dimVol = size(vol);
if numel(dimVol)<3
    error('data should be at least 3D');
end

if numel(patch_dim)~=3
    error('Only 3D patch dimensions supported');
end
patch_dist = floor(patch_dim/2); %distance to center patch

%% 3D patch extraction
if numel(dimVol) == 3    
    out = zeros(numel(vol),patch_dim(1)*patch_dim(2)*patch_dim(3));
    vol = padarray(vol,patch_dist);
    for z=1+patch_dist(3):size(vol,3)-patch_dist(3) %loop through all slices w/o padding
        for pz=1:patch_dim(3) %at every slice, get all patches in z
           ind1 = (z-patch_dist(3)-1)*dimVol(1)*dimVol(2)+1:(z-patch_dist(3))*dimVol(1)*dimVol(2);
           ind2 = (pz-1)*patch_dim(1)*patch_dim(2)+1:pz*patch_dim(1)*patch_dim(2);
           out(ind1,ind2) = im2col(vol(:,:,z-patch_dist(3)+pz-1),[patch_dim(1) patch_dim(2)])';
        end
    end
else
   vol = vol(:,:,:,:); %concatenate last dimensions
   out = zeros(dimVol(1)*dimVol(2)*dimVol(3),patch_dim(1)*patch_dim(2)*patch_dim(3),dimVol(end));
   vol = padarray(vol,patch_dist);
   for p=1:dimVol(4)
        for z=1+patch_dist(3):size(vol,3)-patch_dist(3) %loop through all slices w/o padding
            for pz=1:patch_dim(3) %at every slice, get all patches in z
               ind1 = (z-patch_dist(3)-1)*dimVol(1)*dimVol(2)+1:(z-patch_dist(3))*dimVol(1)*dimVol(2);
               ind2 = (pz-1)*patch_dim(1)*patch_dim(2)+1:pz*patch_dim(1)*patch_dim(2);
               out(ind1,ind2,p) = im2col(vol(:,:,z-patch_dist(3)+pz-1,p),[patch_dim(1) patch_dim(2)])';
            end
        end 
   end
end

end