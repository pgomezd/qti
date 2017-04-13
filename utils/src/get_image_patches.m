function [out] = get_image_patches(data,par)
%% Get image patches
% =========================================================================
% Get patches for matching to dictionary
% =========================================================================
% v1.1: 15.08.16
% v2.1: 30.03.17 - added 2x2 patches. Tag: pgd.2x2
% =========================================================================
% P. Gomez
% =========================================================================

%% get variables
datadims = par.ind.datadims;
P = par.ind.P;
Q = par.ind.Q;
T = par.ind.NTimepoints;
ix = par.ind.ix;
iy = par.ind.iy;
if par.ind.Nslices > 1    
    iz = par.ind.Nslices;
end
% --- 30.03.2017 pgd.2x2 --- %
if P==4
   par.ind.N = (mod(ix,2)+ix)*(mod(iy,2)+iy)/P*iz;
end
% --- 30.03.2017 pgd.2x2 --- %
N = par.ind.N;

p_side = par.ind.p_side;
p_dist = par.ind.p_dist;

%% extract patches/reformat data for matching
if P == 1
    if length(datadims)==3 %2D + Q data
        if par.f.verbose; fprintf('Reformatting 2D data for dictionary matching \n'); end;
    elseif length(datadims)==4 %3D + Q data
        if par.f.verbose; fprintf('Reformatting 3D data for dictionary matching \n'); end; 
    end
    Xp = reshape(data.X,[N,T]);
    if isfield(data,'qmap')
        Qp = reshape(data.qmap,[N,Q]);
    end
% --- 30.03.2017 pgd.2x2 --- %
elseif P == 4
   if length(datadims)==3 %2D + Q data
        if par.f.verbose; fprintf('Extracting multiparmetric 2x2+T distinct image patches \n'); end;
        Xp = zeros((mod(ix,2)+ix)*(mod(iy,2)+iy)/P,P,T);
        if par.f.use_parallel
            Xi = data.X;
            parfor t = 1:T
               Xp(:,:,t) = im2col(squeeze(Xi(:,:,t)),[p_side p_side],'distinct')'; 
            end
        else
            for t = 1:T
               Xp(:,:,t) = im2col(squeeze(data.X(:,:,t)),[p_side p_side],'distinct')'; 
            end
        end

        if isfield(data,'qmap')
            Qp = zeros((mod(ix,2)+ix)*(mod(iy,2)+iy)/P,P,Q);
            qmap = data.qmap;
            if par.f.use_parallel
               parfor q = 1:Q
                    Qp(:,:,q) = im2col(squeeze(qmap(:,:,q)),[p_side p_side],'distinct')'; 
               end
            else
               for q = 1:Q
                    Qp(:,:,q) = im2col(squeeze(qmap(:,:,q)),[p_side p_side],'distinct')'; 
               end
            end
            clear qmap;
        end
        if isfield(data,'mask')
            data.mask = im2col(data.mask,[p_side p_side],'distinct')';
            data.mask = sum(data.mask,2) == P;
        end
   elseif length(datadims)==4 %3D + Q data
       if par.f.kernel3D
            error('Only slice-wise processing with 2x2 patches');
       else
            if par.f.verbose; fprintf('Extracting multiparmetric 2x2+T distinct image patches from 3D volume \n'); end;
            Xp = zeros((mod(ix,2)+ix)*(mod(iy,2)+iy)/P,iz,P,T);
            if par.f.use_parallel
                Xi = data.X;
                parfor t = 1:T
                    for z=1:iz
                        Xp(:,z,:,t) = im2col(squeeze(Xi(:,:,z,t)),[p_side p_side],'distinct')';
                    end
                end
            else
                for t = 1:T
                    for z=1:iz
                        Xp(:,z,:,t) = im2col(squeeze(data.X(:,:,z,t)),[p_side p_side],'distinct')';
                    end
                end
            end
            Xp = reshape(Xp,[(mod(ix,2)+ix)*(mod(iy,2)+iy)/P*iz,P,T]);

            if isfield(data,'qmap')
                Qp = zeros((mod(ix,2)+ix)*(mod(iy,2)+iy)/P,iz,P,Q);
                qmap = data.qmap;
                if par.f.use_parallel
                    parfor q = 1:Q
                       for z = 1:iz
                           Qp(:,z,:,q) = im2col(squeeze(qmap(:,:,z,q)),[p_side p_side],'distinct')';
                       end
                    end
                else
                    for q = 1:Q
                       for z = 1:iz
                            Qp(:,z,:,q) = im2col(squeeze(qmap(:,:,z,q)),[p_side p_side],'distinct')'; 
                       end
                    end
                end
                clear qmap;
                Qp = reshape(Qp,[(mod(ix,2)+ix)*(mod(iy,2)+iy)/P*iz,P,Q]);
            end
            if isfield(data,'mask')
                msk = zeros((mod(ix,2)+ix)*(mod(iy,2)+iy)/P,iz,P);
                cmsk = data.mask;
                if par.f.use_parallel
                    parfor z = 1:iz
                        msk(:,z,:) = im2col(squeeze(cmsk(:,:,z)),[p_side p_side],'distinct')';
                    end
                else
                    for z = 1:iz
                        msk(:,z,:) = im2col(squeeze(cmsk(:,:,z)),[p_side p_side],'distinct')';
                    end 
                end
                msk = reshape(msk,[(mod(ix,2)+ix)*(mod(iy,2)+iy)/P*iz,P]);
                data.mask = sum(msk,2) == P;
                clear msk;
            end
       end
   end 
else
% --- 30.03.2017 pgd.2x2 --- %
    if length(datadims)==3 %2D + Q data
        if par.f.verbose; fprintf('Extracting multiparmetric 2D+T image patches \n'); end;
        Xp = zeros(N,P,T);
        if par.f.use_parallel
            Xi = data.X;
            parfor t = 1:T
               Xp(:,:,t) = im2col(padarray(squeeze(Xi(:,:,t)),[p_dist p_dist]),[p_side p_side])'; 
            end
        else
            for t = 1:T
               Xp(:,:,t) = im2col(padarray(squeeze(data.X(:,:,t)),[p_dist p_dist]),[p_side p_side])'; 
            end
        end

        if isfield(data,'qmap')
            Qp = zeros(N,P,Q);
            qmap = data.qmap;
            if par.f.use_parallel
               parfor q = 1:Q
                    Qp(:,:,q) = im2col(padarray(squeeze(qmap(:,:,q)),[p_dist p_dist]),[p_side p_side])'; 
               end
            else
               for q = 1:Q
                    Qp(:,:,q) = im2col(padarray(squeeze(qmap(:,:,q)),[p_dist p_dist]),[p_side p_side])'; 
               end
            end
            clear qmap;
        end
    elseif length(datadims)==4 %3D + Q data
        % --- 30.05.2016 pgd.3Dkernel --- %
        if par.f.kernel3D
            if par.f.verbose; fprintf('Extracting multiparamteric 3D+T image patches \n'); end;
            Xp =  vol2col(data.X,[p_side p_side p_side]);

            if isfield(data,'qmap')
               Qp  = vol2col(data.qmap,[p_side p_side p_side]);
            end
        else %use a 2D kernel on 3D data %pgd.3Dkernel
            if P == 1
                if par.f.verbose; fprintf('Reformatting 2D slice from 3D data for dictionary matching \n'); end;
            else
                if par.f.verbose; fprintf('Extracting multiparmetric 2D+T image patches from 3D volume \n'); end;
            end
            Xp = zeros(ix*iy,iz,P,T);
            if par.f.use_parallel
                Xi = data.X;
                parfor t = 1:T
                    for z=1:iz
                        Xp(:,z,:,t) = im2col(padarray(squeeze(Xi(:,:,z,t)),[p_dist p_dist]),[p_side p_side])';
                    end
                end
            else
                for t = 1:T
                    for z=1:iz
                        Xp(:,z,:,t) = im2col(padarray(squeeze(data.X(:,:,z,t)),[p_dist p_dist]),[p_side p_side])';
                    end
                end
            end
            Xp = reshape(Xp,[N,P,T]);

            if isfield(data,'qmap')
                Qp = zeros(ix*iy,iz,P,Q);
                qmap = data.qmap;
                if par.f.use_parallel
                    parfor q = 1:Q
                       for z = 1:iz
                            Qp(:,z,:,q) = im2col(padarray(squeeze(qmap(:,:,z,q)),[p_dist p_dist]),[p_side p_side])'; 
                       end
                    end
                else
                    for q = 1:Q
                       for z = 1:iz
                            Qp(:,z,:,q) = im2col(padarray(squeeze(qmap(:,:,z,q)),[p_dist p_dist]),[p_side p_side])'; 
                       end
                    end
                end
                clear qmap;
                Qp = reshape(Qp,[N,P,Q]);
            end 
        end
        % --- 30.05.2016 pgd.3Dkernel --- %
    end
end
%% out
out.Xp = Xp(:,:);
out.X = data.X;

if isfield(data,'qmap')
    Qp = Qp(:,:);
    out.Qp = Qp;
end
if isfield(data,'mask')
    out.mask = data.mask;
end
if par.f.Yout && isfield(data,'Y')  
   out.Y = data.Y; 
end

end