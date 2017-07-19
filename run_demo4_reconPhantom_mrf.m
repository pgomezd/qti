%% Demo 3: Recon Brainweb phantom with MRF recon
% =========================================================================
addpath(genpath('utils'))

%% Setup par
par.f.calcmask = false; %automatic masking
par.f.use_parallel = true;
par.f.match = true; %dictionary matching
par.f.Yout = true; %output raw k-space as well
par.f.pdout = true;
par.f.mtout = false;
par.f.dmout = false;
par.f.qout = true; 
par.f.Xout = true;
par.f.alphaout = false;
par.f.compress = false;
par.recon.psd = 'qti';
par.f.match = true; %match to dictionary
par.ind.mtx_reco = 256; %matrix size
par.f.use_parallel = true; 
par.fp.blockSize = floor(256*256/10);
par.recon.scale_K = 256/25                 ; %ksp scaling (artifical zoom in/out)
par.recon.sp_method = 'nufft'; 
par.f.scaleK = true; %scales k by 1/maximum
par.recon.lsqr_max_iter = 10;
par.recon.admm_max_iter = 250;
par.recon.lsqr_tol = 1e-4;
par.recon.llr_block_dim = 8;
par.recon.admm_lambda = 0.00005;
par.recon.admm_rho = 0.9; 
par.ind.temp_coeff = 8;
par.recon.method = 'STVNNR';
par.f.apply_tempsubspace = true; %compresses data by applying temp subspace

% parameters for STVNNR
par.recon.STVNNR.lambda1 = 9*10^(-4); % parameter for TV term
par.recon.STVNNR.lambda2 = 7*10^(-3); % parameter for low-rank term
par.recon.STVNNR.L = 1/82;   % Step size parameter, should be close to zero. 1/acc_factor is a good choice.
par.recon.STVNNR.IterNo = 1000;  % Maximum number of iterations for the algorithm 


%% Get files
dict =  (['data' filesep 'dict.mat']);
pfile = (['data' filesep 'raw.mat']);
wave =  (['data' filesep 'wave.mat']);

%% Recon
out = recon_qti(dict,pfile,wave,par);

%% Load phantom
load(['data' filesep 'ref.mat'])

%% Compare
close all;
figure;
set(gcf,'units','normalized','outerposition',[0.2 0.2 0.9 0.9]);
titles = {'Ref', 'Recon', 'Difference'};
ylabs = {'T1 (ms)','T2 (ms)'};
for in = 1:6
    subplot(2,3,in)
    switch in
        case 1
            cimg = ref.qmap(:,:,1)*1e3;
        case 2
            cimg = out.qmap(:,:,1)*1e3;
        case 3
            cimg = (ref.qmap(:,:,1)-out.qmap(:,:,1))*1e3;
        case 4
            cimg = ref.qmap(:,:,2)*1e3;
        case 5
            cimg = out.qmap(:,:,2)*1e3;
        case 6
            cimg = (ref.qmap(:,:,2)-out.qmap(:,:,2))*1e3;   
    end
    imagesc(cimg)
    axis image
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    
    colormap(ax,hot);
    cb = colorbar(ax,'south');
    cb.Position(1) = cb.Position(1)+0.01; %move
    cb.Position(3) = cb.Position(3)-0.02; %reduce height
    cb.Position(4) = cb.Position(4)/4; %reduce width
    cb.Color = 'w';
    cb.AxisLocation = 'in';

    if in <= 3
        caxis([0 4000])  
        title(titles{in})
    else
       caxis([0 350])        
    end
    
    if in == 1
        ylabel(ylabs{1})
    elseif in == 4
        ylabel(ylabs{2})
    end
           
end






