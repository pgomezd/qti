function [out,par] = recon_qti(dict,pfile,varargin)
%% QTI recon
% =========================================================================
% [out,par] = recon_qti(dict,pfile,wave,par)
% [out,par] = recon_qti(dict,data,par)
% Input
%       dict                    variable or path to file with the following variables
%           .D                  [L atoms, T(P) elements]            Dictionary
%           .lut                [L atoms, Q(P) maps]                Look-up table (optional, if par.f.match)
%           .V                  [T timepoints, V singular vectors]  Singular vectors (optional, if par.f.compress)
%           .normD              [L atoms, 1]                        norm of dictionary (optional)
%       Either:
%           pfile               [T timeoints, K kspace points, 1,1, Nslices, Ncoils]  Raw Pfile as obtained from a GE scanner (path to file)              
%           wave                                                                               
%               .k              [1, K kspace points]        Complex k-space k = x + iy 
%               .dcf            [1, K kspace points]        Density compensation
%               .npix                                       Nominal matrix resolution of acquisition
%
%       Or:
%               data.Y         [Nx,Ny,(Nz), T timepoints]   Undersampled data 
%               data.U         [Nx,Ny,(Nz)]                 Undersampling mask
%               data.mask      [Nx,Ny,(Nz)]                 Fitting mask, optional
%
%       Or:
%               data.X         [Nx,Ny,(Nz), T timepoints]   Reconstructed image
%               data.qmap      [Nx,Ny,(Nz), Q maps]         Estimated quantitative maps (optional)
%
%       par.        Anything that somehow controls the recon pipeline goes into this struct
%
% %
% Output
%       out.qmap       [Nx,Ny,(Nz),  Q maps]                Quantitative maps  
%       out.X          [Nx,Ny,(Nz),  T timepoints]          Reconstructed image
%       out.Xfit       [Nx,Ny,(Nz),  T timepoints]          Dictionary fit
%       out.pd         [Nx, Ny,(Nz), 1]                     Proton density
%       out.mt         [Nx, Ny,(Nz), Kn nearest neighbors]  Voxel - atom correlation
%       out.dm         [Nx, Ny,(Nz), Kn nearest neighbors]  Selected dictionary entry at every voxel
%       out.mask       [Nx, Ny,(Nz)]                        Calculated mask (if par.f.calcmask)
%
% %   
% Function
%   Recons QTI data for multiple inputs with different methods and implementations
% =========================================================================
% Available implementations:
%   Reconstruction (par.recon.method):
%       par.recon.method        = 'NN'          Match to dictionary with initial image estimate 
%       par.recon.method        = 'TM'          Template matching: fft/nufft recon on undersampled data
%       par.recon.method        = 'ADMM'        CS recon with temporal subpsace + local low rank regularization
%   Dictionary matching (par.dict.type)
%       par.dict.type           = 'cluster'     L x T(P) dictionary from simulated/clustered data
%       par.dict.type           = 'atlas'       L x T(P) dictionary from an atlas of co-registered training subjects
%   Other implementations
%       par.f.compress          = true          SVD Compression of k-space data in the temporal domain
%       par.f.use_parallel      = true          Use parallel toolbox when possible
%   For more flags/recon parameters see setup_mrf_par.m
% =========================================================================
% External tools:
%   J. Fessler's nuFFT toolbox.           /utils/lrt
%   M. Lustig's wrapper for the nuFFT.    /utils/lrt/@NUFFT  
%   J. Tamir's subspace recon.            /utils/t2shuffle
%   B. Hargreaves' EPG simulations.       /utils/epg
% =========================================================================
% External references:
%   MR fingerprinting                     Ma et al. Magnetic Resonance Fingerprinting. Nature 2013.
%   SVD compression for MRF               McGivney et al. SVD Compression for Magnetic Resonance Fingerprinting in the Time Domain. IEEE TMI 2014.
%   T2-shuffling and subspace recon       Tamir et al. T2 Shuffling: Sharp, Multiconstrast, Volumetric Fast Spin-Echo Imaging. MRM 2016.   
%
% References:
%   Spatiotemporal dictionary matching.   Gomez et al. Learning a Spatiotemporal Dictionary for Magnetic Resonance Fingerprinting with Compressed Sensing. MICCAI-PMI 2015.
%   Dictionary clustering.                Gomez et al. 3D Magnetic Resonance Fingerprinting with a Clustered Spatiotemporal Dictionary. IMSMRM 2016.
%   Atlas-based dictionary matching.      Gomez et al. Simultaneous Parameter Mapping, Modality Synthesis, and Anatomical Labeling of the Brain with MR Fingerprinting. MICCAI 2016.
%   ADMM recon for MRF/QTI                Gomez et al. Accelerated parameter mapping with compressed sensing: an alternative to MR fingerprinting. ISMRM 2017
%   Quantitative transient-state imaging. Gomez et al. Ultrafast magnetic resonance imaging and parametric mapping with efficient transient-state encoding. Submitted.
% =========================================================================
% v1.1: 02.06.16
% v1.2: 24.07.16 - moved recon functions to individual files and general optimization
% v2.1: 15.08.16 - merged file with recon_qti file, left functions as aux functions
% v3.1: 15.09.16 - changed forward and backward operations to use operators
% v3.2: 04.10.16 - moved zero-filling/viewsharing and sense estimation into get_sens
% v3.4: 21.11.16 - general debugging, incorporated cartesian iterations
% v4.1: 12.04.17 - removed functions not specific to QTI publication
% =========================================================================
% (c) P. Gomez 2017
% =========================================================================

%% Defaults and setup
if isa(pfile,'char') %input: pfile path
    if isempty(varargin)
       error('Path to waveform needed as third argument if input is path to raw pfile');
    elseif length(varargin)==1
        wave = varargin{1};
        par.f.defaultSetup = false;
    elseif length(varargin)==2
        wave = varargin{1};
        par = varargin{2};
        par.f.defaultSetup = false;
    else
        error('Only path to wave and par can given as additional inputs')
    end
    par.recon.data_in = 'pfile';
    par.dir.pfile = pfile;
    data = pfile;
    clear pfile; %rename to data
elseif isa(pfile,'struct') %input: data struct
    data = pfile;
    clear pfile; %rename to data
    if (~isfield(data,'Y') && ~isfield(data,'U')) 
        if ~isfield(data,'X')
            error('data.X or data.Y/data.U required as input parameters');
        end
    end
    if isempty(varargin)
        par.f.defaultSetup = true; 
    elseif length(varargin)==1
        par = varargin{1};
        par.f.defaultSetup = false; 
    else
       error('Only par can given as additional input if input is data struct') 
    end 
    par.recon.data_in = 'struct';
else
    error('The second argument of the function must contain either the path to a raw pfile or a data struct');
end
if isa(dict,'char') %input: dictionary file to path
    load(dict); %if char, load dict to have file
end

% set-up variables
par = setup_qti_par(dict,data,par);

% prep data
[data,par]=prep_MR_rawdata(data,wave,par);

% get initial image estimate and sens maps
[data,par] = get_sens(data,par);

% get operators
par = get_qti_operators(dict,par);
        
%% Recon with different methods
switch par.recon.method
    case 'NN'
        if par.f.match
            out = par.oper.D_match(data); %projection/match to dictionary
        else
            out = data; %output data estimate with no matching
        end
    case 'TM'
        data.X = par.oper.A_adj(data.Y); %recon image
        if par.f.match
            out = par.oper.D_match(data); %projection/match to dictionary
        end
    case 'IPA'
    case 'ADMM'
        out = admm_recon(data,par);
end

end

