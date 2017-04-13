function out = qti_dtm_fit(dict,data,par)
%% qti fit with DTM
% =========================================================================
% out = qti_dtm_fit(dict,data,par)
% %
%   Input
%         
%       dict
%           .D          [K entries, PT elements]   dictionary
%           .normD      [K entries, 1]             norm of dictionary
%           .lut        [K entries, PQ biological parameter patches]  
%           .lutN       [K entries, PQ biological parameter patches] normalized (optional)
%           .V          [K entries, V singular vectors] temporal subspace of dict (optional) 
%       data
%           .x          [N voxels,  PT elements] image data
%           .q          [N voxels,  PQ params]   estimated parameters (optional)
%       par
%           .rp.alpha   [1, PT elements]    image space weights
%           .rp.beta    [1, PQ params]      parameter space weights
%       	.rp.gamma   {0,1} Regularization between image (gamma = 0) and parameter space (gamma = 1) 
%
%   Output
%       out.X           [N voxels, PT elements]            fit of data Xf=norm(Xi)/norm(D)*D
%       out.qmap        [N voxels, PQ params]              estimated biological parameters
%       out.mt          [N voxels, Kn nearest neighbors]   matching result (voxel-atom correlation)
%       out.dm          [N voxels, Kn nearest neighbors]   selected dictionary entries in L
%       out.pd          [N voxels, 1]                      proton density  
%
%   Function
%       Finds the highest correlated atom in the dictionary to the data
%       by projection max <Di,x>
% =========================================================================
% Log
% v1.1: 20.10.15
% v1.2: 31.10.15 - replaced patch extraction and update with vol2col and col2vol.
% v1.3: 22.11.15 - code integration into 1 script.
% v2.1: 08.02.16 - included option for atlas-based dictionaries. Tag: pgd.atlas
% v2.1: 16.02.16 - extended to k-nearest neighbor matching. Tag: pgd.knn
% v2.2: 21.02.16 - added a second iteration to match in parameter space. Tag: pgd.pmiter
% v2.3: 07.03.16 - modifications to run matching on GPU. Tag: pgd.pmiter
% v2.4: 05.04.16 - optimization and general debugging.
% v2.5: 31.05.16 - added flag for 2D kernels on 3D data. Tag: pgd.3Dkernel
% v3.1: 15.05.16 - removed bpar and migrated all to par
% v4.1: 15.09.16 - split code into separate files, use operator instead
% v4.2: 26.09.16 - re-adjusted PD calculations. Tag: pgd.pdcalc
% v4.5: 30.09.16 - moved lut and lutN into dict struct, removed lut centering
% =========================================================================
% P. Gomez
% =========================================================================

%% Set-up
x = data.x;

%% Normalization
D = dict.D;
if isfield(dict,'normD'); normD = dict.normD; end;    
if ~isfield(dict,'lut')
    error('dict.lut required for parameter estimation'); 
end
if ~exist('normD','var')
    warning('normD not found, normalizing dictionary \n');
    normD = zeros(1,size(D,1));
    if par.f.use_parallel
        parfor l = 1:size(D,1)
            normD(l)=norm(D(l,:));
            D(l,:)=D(l,:)/normD(l);
        end
    else
        for l = 1:size(D,1)
            normD(l)=norm(D(l,:));
            D(l,:)=D(l,:)/normD(l);
        end
    end
    D(isnan(D))=0;
end

if (par.f.normLUT || ~isfield(par,'lutN')) && isfield(data,'q') && par.rp.gamma > 0  %normalize the LUT weighted by par.rp.beta
    if par.f.normLUT
        if par.f.verbose; fprintf('Normalizing lut \n'); end;
    elseif ~isfield(dict,'lutN')
        warning('dict.lutN not found, normalizing lut');
    end
    lutN = zeros(size(dict.lut));
    if ~isempty(par.rp.beta)
        beta = par.rp.beta;
    else
        beta = [];
    end
    lut = dict.lut;
    if par.f.use_parallel
        parfor l = 1:size(D,1)
            if isempty(beta)
               lutN(l,:) = lut(l,:)/(norm(lutN(l,:))); %normalize
            else
               lutN(l,:) = lut(l,:)/(sqrt(sum(beta.*(lutN(l,:).^2)))); %normalize with weights
            end
        end
    else
        for l = 1:size(D,1)
             if isempty(beta)
               lutN(l,:) = lut(l,:)/(norm(lutN(l,:))); %normalize
            else
               lutN(l,:) = lut(l,:)/(sqrt(sum(beta.*(lutN(l,:).^2)))); %normalize with weights
            end
        end
    end
    dict.lutN = lutN;
    clear lut*;
    clear beta;
end

    
%% DTM Fit
blockSize = par.fp.blockSize;
iter = ceil(size(x,1)/blockSize);
if par.f.use_gpu
    warning('PD calculation not available on GPU');
    if isfield (data,'q') %runs on GPU
        [mt,dm] = qti_dtm_gpu(D,x,par,data.q,dict.lutN);
    else
        [mt,dm] = qti_dtm_gpu(D,x,par); 
    end
else
    mt = zeros(size(x,1),par.ind.knn);
    dm = zeros(size(x,1),par.ind.knn);
    % --- 17.08.16-05.09. pgd.pdcalc --- %
    pd = zeros(size(x,1),1);
    X = zeros(size(x));
    % --- 17.08.16-05.09. pgd.pdcalc --- %
    if par.f.verbose; fprintf('Matching data \n'); end;
    for i = 1:iter
        if par.f.verbose_iter; fprintf(['Matching block ', num2str(i) '/' num2str(iter), '\n']); end;
        if i<iter
            cind = (i-1)*blockSize+1:i*blockSize;    
        else
            cind = (i-1)*blockSize+1:size(x,1); 
        end
        % --- 21.02.16 pgd.pmiter --- %
        % --- 31.07.2016 pgd.atomind --- %
        if ~isempty(par.atomind)
            ip=(D*ctranspose((bsxfun(@times,x(cind,:),par.rp.alpha)))).*par.atomind(:,cind);      
        else
            ip=D*ctranspose(x(cind,:)); 
        end
        % --- 31.07.2016 pgd.atomind --- %
            if isfield(data,'q') && par.rp.gamma > 0
                ipq = (dict.lutN*(bsxfun(@times,data.q(cind,:),par.rp.beta)).').*par.atomind(:,cind);
                if par.ind.knn > 1
                    [match,dictentry] = sort((1-par.rp.gamma)*abs(ip)+par.rp.gamma*ipq,'descend'); %review: sorting is slow - needs better solution
                    % --- 18.02.16 pgd.knn --- %
                    mt(cind,:) = match(1:par.ind.knn,:)';
                    dm(cind,:) = dictentry(1:par.ind.knn,:)';
                    % --- pgd.knn --- %
                else
                    % --- 31.07.2016 pgd.atomind --- %
                    if ~isempty(par.rp.gamma)
                        [mt(cind),dm(cind)] = max((1-par.rp.gamma)*abs(ip)+par.rp.gamma*ipq,[],1);
                    else
                        [mt(cind),dm(cind)] = max(abs(ip),[],1);
                    end
                    % --- 31.07.2016 pgd.atomind --- %
                end
            else
                % --- 18.02.16 pgd.knn --- %
                if par.ind.knn > 1
                    [match,dictentry] = sort(abs(ip),'descend'); %review: sorting is slow - needs better solution
                    mt(cind,:) = match(1:par.ind.knn,:)';
                    dm(cind,:) = dictentry(1:par.ind.knn,:)';
                else
                % --- 18.02.16 pgd.knn --- % 
                    [mt(cind),dm(cind)] = max(abs(ip),[],1);
                end
            end
            % --- 21.02.16 pgd.pmiter --- %
            % --- 17.08.16-26.09. pgd.pdcalc --- %
            for ic = 1:length(cind)
                pd(cind(ic)) = ip(dm(cind(ic),1),ic);
                X(cind(ic),:) = pd(cind(ic)).*D(dm(cind(ic),1),:); %re-scale X as in Davies et al.
                pd(cind(ic)) = pd(cind(ic))./normD(dm(cind(ic),1)); %get PD after signal scaling. 
            end       
    end
    % --- 17.08.16-26.09. pgd.pdcalc --- %
end            
% clear ip*;
clear match;
clear dictentry;

%% Qmap
% --- 18.02.16 pgd.knn --- %
qmap = zeros(size(dm,1),size(dict.lut,2),par.ind.knn);
for kn = 1:par.ind.knn
    qmap(:,:,kn) = bsxfun(@times,mt(:,kn),dict.lut(dm(:,kn),:));
end
% --- 18.02.16 pgd.knn --- %
qmap = bsxfun(@rdivide,sum(qmap,3),sum(mt,2)); %weighted average of selected params

%% out
out.Xfit = X;
out.X = data.X;
out.Xp = data.Xp;
out.mask = data.mask;

if par.f.qout; out.qmapfit = qmap; end;
if par.f.mtout; out.mtfit = mt; end;
if par.f.dmout; out.dmfit = dm; end;
if par.f.pdout; out.pdfit = pd; end;

if par.f.Yout && isfield(data,'Y')  
   out.Y = data.Y; 
end

end


function [mt,dm] = qti_dtm_gpu(D,x,par,varargin)
%% qti fit with DTM on GPU
% =========================================================================
% Run qti fit on GPU device
% =========================================================================
% v1.1: 02.11.15
% =========================================================================
% P. Gomez
% =========================================================================

%% Get qmaps and LUT if available
if length(varargin)==1
    error('q and lutN required for data-driven updates');
elseif length(varargin)==2
    par.f.qmapUpdate = 1;  
    q = varargin{1};
    lutN = varargin{2};
end

%% Memory allocation and problem slicing to fit to GPU memory
maxMem = par.fp.maxMem;
memAtom = size(D,2)*8; %8 bytes complex single
L = floor(maxMem/memAtom*.45); %45% of the memory goes to the Dictionary
memL = L*8; %8 bytes complex single 
N = floor(maxMem/(memAtom+memL)*.45); %define N s.t. memory allows for ip = d*x';

%define number of segments in L and N
segL = ceil(size(D,1)/L);
segN = ceil(size(x,1)/N);

%% DTM Fit
% perform fit on each individual segment
d_mt = gpuArray(zeros(segL,size(x,1)));
d_dm = gpuArray(zeros(segL,size(x,1)));
for n=1:segN
%         fprintf(['segment n = ' num2str(n) '/' num2str(segN) '\n']);
    if n<segN
        nind = (n-1)*N+1:n*N;
    else
        nind = (n-1)*N+1:size(x,1);
    end
    d_x = gpuArray(x(nind,:));
    % --- 07.03.16 pgd.pmiter --- %
    % --- 31.07.2016 pgd.atomind --- %
    if ~isempty(par.rp.alpha)
        d_a = gpuArray(par.rp.alpha);
        if  par.f.qmapUpdate
                d_q = gpuArray(q(nind,:));
                d_b = gpuArray(par.rp.beta);
        end
    end
    % --- 31.07.2016 pgd.atomind --- %
    % --- 07.03.16 pgd.pmiter --- %
    for l = 1:segL 
        if l<segL
            lind = (l-1)*L+1:l*L;
        else
            lind = (l-1)*L+1:size(D,1);
        end
        d_D = gpuArray(D(lind,:));
        % --- 31.07.2016 pgd.atomind --- %
        if ~isempty(par.atomind)
            d_atind = gpuArray(par.atomind(lind,nind));
            % --- 07.03.16 pgd.pmiter --- %
            d_ip = (d_D*ctranspose((bsxfun(@times,d_x,d_a)))).*d_atind;
            if  par.f.qmapUpdate
                d_lutN = gpuArray(lutN(lind,:));
                d_ipq = (d_lutN*(bsxfun(@times,d_q,d_b))').*d_atind;
                d_gamma = gpuArray(par.rp.gamma);
                if par.ind.knn > 1
                    [d_match,d_dictentry] = sort((1-d_gamma)*abs(d_ip)+d_gamma*d_ipq,'descend'); %sorting is slow - needs better solution
                    d_mt(l,nind,:) = d_match(1:par.ind.knn,:)';
                    d_dm(l,nind,:) = lind(d_dictentry(1:par.ind.knn,:)');
                else
                    [d_mt(l,nind),d_cind] = max((1-d_gamma)*abs(d_ip)+d_gamma*d_ipq,[],1);
                    d_dm(l,nind) = lind(d_cind);
                end
                clear d_lutN;
                clear d_ipq;
                clear d_gamma;
            else
            % --- 07.03.16 pgd.pmiter --- %
                % --- 18.02.16 pgd.knn --- %
                if par.ind.knn > 1
                    [d_match,d_dictentry] = sort(abs(d_ip),'descend'); %sorting is slow - needs better solution
                    d_mt(l,nind,:) = d_match(1:par.ind.knn,:)';
                    d_dm(l,nind,:) = lind(d_dictentry(1:par.ind.knn,:)');
                % --- 18.02.16 pgd.knn --- %
                else
                    [d_mt(l,nind),d_cind] = max(abs(d_ip),[],1);
                    d_dm(l,nind) = lind(d_cind);
                end
            end
        else
            d_ip=d_D*ctranspose(d_x);
            [d_mt(l,nind),d_cind] = max(abs(d_ip),[],1);
            d_dm(l,nind) = lind(d_cind);
        end
        % --- 31.07.2016 pgd.atomind --- %

        clear d_ip;
        clear d_D;
        clear d_atind;
    end
    clear d_x;
end

%% Gather data and get max per segmentS
[d_mseg,d_indseg] = max(d_mt,[],1);
mt = gather(d_mseg)';
dmg = gather(d_dm);
indseg = gather(d_indseg);
dm = zeros(size(x,1),1);
for cn = 1:size(x,1)  
    dm(cn) = dmg(indseg(cn),cn);
end

end