function par = setup_qti_par(dict,data,par)
% % Setup qti par for recon
% =========================================================================
% Default values in par structure for recon
% =========================================================================
% v1.1: 30.05.16
% =========================================================================
% (c) P. Gomez 2016
% =========================================================================

%% Flags
if ~isfield(par,'f'); par.f.f = true; end;  
if ~isfield(par.f,'compress'); par.f.compress = false; end; %compress the dictionary
if ~isfield(par.f,'calcmask'); par.f.calcmask = true; end; %automatic mask calculation based on data
if ~isfield(par.f,'match'); par.f.match = false; end; %dictionary matching
if ~isfield(par.f,'use_parallel'); par.f.use_parallel = false; end; %uses parallel toolbox when possible
if ~isfield(par.f,'use_gpu'); par.f.use_gpu = false; end; %uses gpu device when possible
if ~isfield(par.f,'debug'); par.f.debug = false; end; %debugging, outputs intermediate recon results (raw data, coil kernels,...)
if ~isfield(par.f,'verbose'); par.f.verbose = true; end; %verbose - output comments
if ~isfield(par.f,'verbose_iter'); par.f.verbose_iter = false; end; %verbose - output iterative comments
if ~isfield(par.f,'scaleK'); par.f.scaleK = false; end %if true, scales raw k-space
if ~isfield(par.f,'apply_tempsubspace'); par.f.apply_tempsubspace = false; end %if true, applies temporal subspace projection in recon operator
        
%% Recon 
if ~isfield(par,'recon'); par.recon.f = true; end;  
if ~isfield(par.recon,'spks'); par.recon.spks = false; end; %spike removal
if ~isfield(par.recon,'mask_thresh'); par.recon.mask_thresh = 70; end; %percentil threshold for automatically masking data [0-100]
if ~isfield(par.recon,'method'); par.recon.method = 'TM'; end % reconstruction method on data.  'NN' no recon, 'TM' (template matching), 'IPA' (compressed sensing w/IPA), 'ADMM' (temporal subspace with LLR regularization based on ADMM) 
if ~isfield(par.recon,'sp_method'); par.recon.sp_method = 'nufft'; end; %spiral recon method. options: 'grid', 'nufft'
if ~isfield(par.recon,'raw_method'); par.recon.raw_method = 'zerofill'; end; %method for filling the k-space. options: 'zerofill', 'viewshare', 'collapse'
if ~isfield(par.recon,'psd'); par.recon.psd = 'qti'; end; % for data pre-processing. options: 'qti', 'rufis','sim'
if ~isfield(par.recon,'scale_K'); par.recon.scale_K = 1; end; % k-space trajectory scaling - artificial zooming in/out of image.
if ~isfield(par.recon,'nechoes'); par.recon.nechoes = 5; end; % number of virtual echoes to bin data into

%% Indexes 
if ~isfield(par,'ind'); par.ind.f = true; end;  
if ~isfield(par.ind,'mtx_reco'); par.ind.mtx_reco = 128; end; %matrix dims for recon
if ~isfield(par.ind,'Nx'); par.ind.Nx = par.ind.mtx_reco; end; 
if ~isfield(par.ind,'Ny'); par.ind.Ny = par.ind.mtx_reco; end;
if ~isfield(par.ind,'Tv'); par.ind.Tv = 10; end; %temporal dictionary compression with SVD
if ~isfield(par.ind,'temp_coeff'); par.ind.temp_coeff = 4; end; % number of subspace coefficients

%% Dictionary
if ~isfield(par,'dict'); par.dict.f = true; end;  
if ~isfield(par.dict,'type'); par.dict.type = 'cluster'; end; %dictionary type, options: cluster, atlas

%% Method-specific flags/params/indexes/check
switch par.recon.method   
    case 'NN'

    case 'TM'
        if ~par.f.match
            warning('par.recon.method = TM requires par.f.match, setting flag')
            par.f.match = true;
        end
    case 'IPA'
        if ~isfield(par.f,'adaptiveStepSearch'); par.f.adaptiveStepSearch = false; end %if true, searches for the best step size in every iteration
        if ~isfield(par.recon,'iter');  par.recon.iter = 10; end; % CS iterations
        if ~isfield(par.recon,'kappa'); par.recon.kappa = 0.99; end; % recon threshold
        if ~isfield(par.recon,'mu');    par.recon.mu = 32; end; % step size 
        if ~isfield(par.recon,'tol');   par.recon.tol = 1e-4; end; % iteration tolerance diff
        
        if ~par.f.match
            warning('par.recon.method = IP requires par.f.match, setting flag')
            par.f.match = true;
        end
    case 'ADMM'
        if ~isfield(par.recon,'llr_block_dim'); par.recon.llr_block_dim = 8; end; % block dims for LLR regularization (must be divisible by par.ind.mtx_reco)
        if ~isfield(par.recon,'admm_max_iter'); par.recon.admm_max_iter = 50; end; % ADMM iters
        if ~isfield(par.recon,'admm_rho'); par.recon.admm_rho = 0.1; end; % ADMM rho for data consitency term
        if ~isfield(par.recon,'admm_lambda'); par.recon.admm_lambda = 0.04; end; % ADMM lambda for regularization term
        if ~isfield(par.recon,'lsqr_max_iter'); par.recon.lsqr_max_iter = 10; end; % LSQR iters
        if ~isfield(par.recon,'lsqr_tol'); par.recon.lsqr_tol = 1e-4; end; % LSQR tol
        
    case 'STVNNR'
       if ~isfield(par.recon.STVNNR,'lambda1'); par.recon.STVNNR.lambda1 = 5*10^(-4); end;
       if ~isfield(par.recon.STVNNR,'lambda2'); par.recon.STVNNR.lambda2 = 5*10^(-3); end;
       if ~isfield(par.recon.STVNNR,'L'); par.recon.STVNNR.L= 1/L; end;
       if ~isfield(par.recon.STVNNR,'IterNo'); par.recon.STVNNR.IterNo = 400; end;
    otherwise
        error('Unknown par.recon.method. Options: NN, TM, IPA, ADMM, STVNNR')
end

switch par.recon.sp_method
    case 'grid'
        
    case 'grid_cart'
        
    case 'nufft'
        if ~isfield(par.recon,'Jk'); par.recon.Jk = [6 6]; end; % interpolation kernel for nuFTT
        
    case 'nufft_cart'
        if ~isfield(par.recon,'Jk'); par.recon.Jk = [6 6]; end; % interpolation kernel for nuFTT
        
    otherwise
        error('Unknown par.recon.sp_method. Options: grid, grid_cart, nufft, nufft_cart')   
end

switch par.recon.raw_method
    case 'zerofill'
                
    case 'viewshare'
        if ~isfield(par.recon,'vs_boundary'); par.recon.vs_boundary = 'closest_frame_int_gb'; end % viesharing boundary method. options: no_share, closest_int, closest_frame, closest_frame_int_gb

    case 'collapse'
        
    otherwise
        error('Unknown par.recon.raw_method. Options: zerofill, viewshare, collapse') 
end

switch par.recon.psd
    case 'qti'
        if ~isfield(par,'grad'); par.grad.f = true; end;  
        if ~isfield(par.grad,'delay'); par.grad.delay = 0; end;
    case 'rufis'
        if ~isfield(par.recon,'oversamp'); par.recon.oversamp = 2; end % data oversampling factor
        if ~isfield(par.f,'sort'); par.f.sort = false; end % radial data sorting
    case 'cart' %cartesian data

    case 'sim' %simulated data - Brainweb phantom

    otherwise
        error('Unknown par.recon.psd. Options: qti, rufis, cart, sim')         
end
                
%% Fitting params
if ~isfield(par,'fp'); par.fp.f = true; end;
if par.f.use_gpu
    if ~isfield(par.fp,'memMax'); par.fp.par.memMax = 0.8; end; %memory to be occupied when fitting data
end

%% Dictionary specific params
% --- pgd.atlas --- %
if strcmp(par.dict.type,'atlas')
    % Inputs
    if ~isfield(par,'amask'); error('par.amask must be given as an input if par.dict.type = atlas');end
    if ~isfield(par.ind,'w_side'); error('par.ind.w_side should be given as an input if par.dict.type = atlas'); end
    
    %Search window
    par.ind.w_dist = floor(par.ind.w_side/2);
    if length(par.ind.w_side) ~=length(datadims)-1
        error('par.ind.w_side should be an array with %d elements',length(datadims)-1);
    end
    if numel(par.ind.w_side)==2
        par.ind.W = par.ind.w_side(1)*par.ind.w_side(2);
    elseif numel(par.ind.w_side)==3
       par.ind.W = par.ind.w_side(1)*par.ind.w_side(2)*par.ind.w_side(3);
    end
end
% --- pgd.atlas --- %

%% Dictionary matching params
if par.f.match
    if ~isfield(dict,'lut') 
        error('dict.lut required as an input if par.f.match = true'); 
    end
    
    if ~isfield(par.f,'Yout'); par.f.Yout = false; end; %output flags on data: Y = k-space space
    if ~isfield(par.f,'Xout'); par.f.Xout = true; end; %output flags on data: X = image space
    if ~isfield(par.f,'qout'); par.f.qout = true; end; %output flags on data: q = parameter maps
    if ~isfield(par.f,'mtout'); par.f.mtout = true; end; %output flags on data: mt = match to dictionary: voxel - atom correlation
    if ~isfield(par.f,'dmout'); par.f.dmout = true; end; %output flags on data: dm = selected dictionary entry at every voxel
    if ~isfield(par.f,'alphaout'); par.f.alphaout = false; end; %output flags on data: alpha = reconstructed temporal coefficients
    if ~isfield(par.f,'pdout'); par.f.pdout = true; end; %output flags on data: pd = proton density
    if ~isfield(par.f,'normX'); par.f.normX = false; end; %normalize the input signal
    if ~isfield(par.f,'normLUT'); par.f.normLUT = false; end;  %normalize the parameter look-up table
    if ~isfield(par.f,'kernel3D'); par.f.kernel3D = true; end; % use a 3D kernel for the spatiotemporal dictionary. pgd.3Dkernel
    if ~isfield(par.f,'windowpatches'); par.f.windowpatches = false; end %get patch-based representation of the entire window (if par.dict.type = atlas)
    if ~isfield(par.recon,'update'); par.recon.update = 'voxel'; end % update step for spatiotemporal dictonary. Options: patch, voxel
    
    % Memory allocation to divide data into blocks
    if ~isfield(par.fp,'blockSize')
        if ispc && par.f.use_gpu %set number of data blocks based on available memory
            mem = memory;
            maxAllowedArray = mem.MaxPossibleArrayBytes*par.fp.memMax; %use certain % of max possible memory
            if isa(data.X,'single')
                    dataBytes = 4*2; %complex single data
            elseif isa(data.X,'double')
                dataBytes = 8*2; %complex double data
            else
                error('data must be single or double precision');
            end
            maxVoxels = floor(maxAllowedArray/(size(par.lut,1)*dataBytes)); %

            if maxVoxels > size(data.X,1)
                par.fp.blockSize = size(data.X,1); %block size covers all of data 
            else
                par.fp.blockSize = maxVoxels; %block size constrained to max voxels - given by % of max possible mem.
            end
        else
            par.fp.blockSize = par.ind.mtx_reco*par.ind.mtx_reco; %constrain blocks to 1 slice
        end
    end
    
    %Index setting: lut, q and patches
    if ~isfield(par.ind,'p_side'); par.ind.p_side = 1; end; %size of patch side for spatiotemporal dictionary matching
    if ~isfield(par.ind,'knn'); par.ind.knn = 1; end % number of nearest neighbors to find in dictionary
    par.ind.p_dist = floor(par.ind.p_side/2);
    par.ind.L = size(dict.lut,1);
    if par.ind.p_side > 1
        if isfield(par.ind,'datadims')
            if par.f.kernel3D
                par.ind.P = par.ind.p_side^(length(par.ind.datadims)-1); 
            else
                par.ind.P = par.ind.p_side^2; %2D kernel on 3D data
            end
        else
            warning('unkwown par.ind.datadims, setting 2D kernel for spatiotemporal dictionary matching')
            par.f.kernel3D = false;
            par.ind.P = par.ind.p_side^2;
        end
    else
       par.ind.P = 1;  
    end
    if mod(size(dict.lut,2),par.ind.P) ==0
        par.ind.Q = size(dict.lut,2)/par.ind.P;
    else
        error('size(dict.lut,2) should be == par.ind.P*par.ind.Q'); 
    end
   
    % Dictionary specific params
    if strcmp(par.dict.type,'atlas')
        % --- pgd.atomind --- %
        par.rp.alpha = ones([1 size(dict.lut,2)]); %equal weights to all
        par.rp.beta =  ones([1 size(dict.lut,2)]); %equal weights to all
        par.rp.gamma = 0;
        % --- pgd.atomind --- %
    elseif strcmp(par.dict.type,'cluster')
        % --- pgd.atomind --- %
        par.atomind =  [];
        par.rp.alpha = []; %equal weights to all
        par.rp.beta =  []; %equal weights to all
        par.rp.gamma = [];
        % --- pgd.atomind --- %
    end
end

%% Parameter check
if (par.recon.mask_thresh < 60 && par.f.calcmask)
   warning('par.recon.mask_thresh recommended between 65 - 85');
end

if par.f.compress && par.f.apply_tempsubspace
    error('par.f.compress incompatible with par.f.apply_tempsubspace, choose one');
end

if par.f.compress
    if isempty(dict)
        par.f.compress_with_data = true; %use data to project onto subspace
    elseif ~isfield(dict,'V') 
        error('dict.V needed for dictionary compression');
    else
        par.f.compress_with_data = false;
        if par.ind.temp_coeff > size(dict.V,2)
           warning('max par.ind.temp_coeff is %d',size(dict.V,2));
           par.ind.temp_coeff = size(dict.V,2);
        end
        par.V = dict.V; %write to par for later use % review
    end
end

if par.f.apply_tempsubspace
    if isempty(dict)
        par.f.compress_with_data = true; %use data to project onto subspace
    elseif ~isfield(dict,'V') 
        error('dict.V needed for temporal subspace compression');
    else
        par.f.compress_with_data = false;
        if par.ind.temp_coeff > size(dict.V,2)
           warning('max par.ind.temp_coeff is %d',size(dict.V,2));
           par.ind.temp_coeff = size(dict.V,2);
        end
        par.V = dict.V; %write to par for later use % review
    end
end


%% Checks for release version
if par.f.use_gpu
    error('GPU implementation not available in this version')
end

switch par.recon.method
    case 'IPA'
        error('par.recon.method=IPA not available in this version')
end

switch par.recon.psd  
    case 'qti'
        
    otherwise
        error('Only par.recon.psd = qti available in this version')
end
switch par.recon.raw_method
    case 'zerofill'
        
    otherwise
        error('Only par.recon.raw_method=zerofill available in this version')
end

switch par.recon.sp_method
    case 'nufft'
        
    otherwise
       error('Only par.recon.sp_method=nufft available in this version') 
end

    