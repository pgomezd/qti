% % Prep_MR_rawdata
% =========================================================================
% [out,par] = prep_MR_rawdata(pfile,wave,par)
%   Input
%       pfile        [T timepoints (Echo), K kspace points (Readout), 1,1, Nslices, Ncoils] raw data as obtained from a GE scanner                
%       wave
%        Either:
%            .k      [1, K kspace points]           Complex k-space k = x + iy
%            .dcf    [1, K kspace points]           Density compensation
%            .npix                                  Nominal matrix resolution of acquisition
%        Or:
%            wave    [3, spokeLength, numSpokes]    3D k-space coordinates for par.recon.psd = rufis
%
%       par.        Anything that somehow controls the recon pipeline goes into this struct
% 
%   Output
%       Either:
%           data.Y      [Nslices,Ncoils, Full K space, NTimepoints] Raw data             
%       Or:
%           data.Y      [Nx,Ny,(Nz), C coils, T timepoints]         Raw data on Carteisan grid
%   Function
%       Preprocces raw data
%               case 'nufft'
%                   Saves nuFFT operator for iterations
%                   Outputs raw data
% =========================================================================
% v1.1: 02.06.16
% v1.1: 29.07.16 - leaves data in raw format
% v2.1: 05.09.16 - applies dcf on raw data for nuFFT operator
% v3.1: 03.10.16 - switch case for general gridding
% v3.2: 04.10.16 - flag for processing of fidall + rufis data
% v3.3: 13.10.16 - switch for par.recon.psd with individual cases
% =========================================================================
% (c) P. Gomez 2016
% =========================================================================
function [out,par] = prep_MR_rawdata(pfile,wave,par)

%% Load raw data
rawdata=load(pfile);
raw = rawdata.raw;

switch par.recon.psd   
    case 'qti'
        par.ind.NTimepoints = size(raw,1); %number of frames / temporal points
        %% Get waveform variables
        wf = load(wave);          % load waveform

        %% Pre-process raw data
        k = wf.k;
        ind = wf.ind;
        if ismatrix(raw)
            nexc = size(raw,1);
            Nslices = 1;
            Ncoils = 1;
        else
            [nexc,~,~,~,Nslices,Ncoils] = size(raw,1); 
        end
        nintl = size(ind,1);     % #of interleaves for encoding one image
        ns1 = sum(ind(:));
        tmp = zeros(1,size(k,2),2);
        tmp(1,:,1) = real(k); tmp(1,:,2) = imag(k);
        k = tmp;
        
        %% selecting indexed raw data
        nrep = ceil(nexc/nintl); % pgd: select always entire interleaves
        dsamp=sum(ind(1,:));    % get number of actual samples per interleave
        NTimepoints = nexc; 
        dtmp = zeros(nrep,ns1,1,1,Nslices,Ncoils);
        for lr=1:nrep
            for lsli=1:Nslices
                for lc=1:Ncoils
                    intl_ind = (lr-1)*nintl+(1:nintl);
                    intl_ind = intl_ind(intl_ind <=nexc); %pgd: make sure intl are inside nexc
                    tmp = raw(intl_ind,:,1,1,lsli,lc).';
                    dtmp(lr,1:length(intl_ind)*dsamp,1,1,lsli,lc) = tmp(ind(1:length(intl_ind),:).');
                end
            end
        end
        raw = dtmp;
        nexc = nrep;
        if any(isnan(raw(:))) || isempty(raw) || any(isinf(raw(:))) 
            error('raw contains nan,inf or empty'); 
        end

        %% pgd: Reformat raw
        if par.f.verbose; fprintf('Reformatting raw data \n'); end;
        switch par.recon.raw_method
            case 'collapse'
            otherwise
                % raw in format: [Nslices, Ncoils, K sampled points, NTimepoints] 
                cdims = size(raw);
                raw = reshape(raw,[nexc dsamp nintl cdims(3:end)]); %separate dsamp from nintl
                raw = permute(raw,[6 7 2 3 1 4 5]); %get rid on singleton dims
                raw = raw(:,:,:,:); %concatenate nexc with nintl
                raw = raw(:,:,:,1:NTimepoints); %ger rid of extra timepoints
        end
        if any(isnan(raw(:))), error('raw raw contains NaN'); end
        if any(isinf(raw(:))), error('raw raw contains inf'); end

        %% pgd: re-scale k-space
        if par.f.scaleK
            if par.f.verbose; fprintf('Scaling raw data \n'); end;
            par.recon.scale = max(raw(:)); %scale by max;
            raw = raw/par.recon.scale; 
        end
        
        %% Get output indexes
        par.ind.nint = nintl;
        [par.ind.Nslices, par.ind.Ncoils, par.ind.dsamp,par.ind.NTimepoints] = size(raw);  
        if par.ind.Nslices == 1
            par.ind.datadims = [par.ind.mtx_reco par.ind.mtx_reco par.ind.NTimepoints];
            par.ind.ix = par.ind.mtx_reco;
            par.ind.iy = par.ind.mtx_reco;
            par.ind.N = par.ind.mtx_reco^2;
        else
            par.ind.datadims = [par.ind.mtx_reco par.ind.mtx_reco par.ind.Nslices par.ind.NTimepoints];
            par.ind.N = par.ind.mtx_reco^2*par.ind.Nslices;
            par.ind.ix = par.ind.mtx_reco;
            par.ind.iy = par.ind.mtx_reco;
            par.ind.iz = par.ind.Nslices;
        end

        %% get dcf
        % review here
        switch par.recon.raw_method
           case 'zerofill'
               par.recon.dcf = reshape(wf.dcf(1:dsamp),[1 1 dsamp]);
           case 'viewshare'
               error('viesharing not available in this version');
           case 'collapse'
                error('collapsing not available in this version');
        end
      
        % adjust for compression
        if par.f.compress 
            par.ind.NTimepoints = par.ind.Tv;
        end
        
        
        %% pgd: out operators
        switch par.recon.sp_method 
            case 'nufft'
                if par.f.verbose; fprintf('Getting nufft operator \n'); end;
                % Get and store nufft operator
                ck = k(1,:,1) + 1i*k(1,:,2); %get complex k-space
                par.recon.Gn = NUFFT(ck,1,1,[0 0],[par.ind.mtx_reco par.ind.mtx_reco],2,par.recon.Jk);  %get nufft operator, no dcf: NUFFT(ck,1,...), dcf: NUFFT(ck,wf.dcf,...)             
        end
        
        %% pgd: out raw 
        raw(isnan(raw))=0;
        raw(isinf(raw))=0;
        par.recon.kmask = raw ~= 0; %undersampling mask
        out.Y = raw;        

    case 'rufis'
        
    case 'cart'
        
    case 'sim'
end

