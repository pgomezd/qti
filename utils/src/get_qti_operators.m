function par = get_qti_operators(dict,par)
% % Get qti operators
% =========================================================================
% Get operators for forward/inverse modelling
% =========================================================================
% v1.1: 15.09.16
% =========================================================================
% (c) P. Gomez 2016
% =========================================================================

%% Sampling mask
P_for = @(y) y.*par.recon.kmask; % undersampling mask must be of same dims as y 
par.oper.P_for = P_for;

%% DCF 
switch par.recon.psd
    case 'qti'
        switch par.recon.sp_method
            case 'grid'
            case 'grid_cart'
             case 'nufft'
                D_adj = @(y) bsxfun(@times,par.recon.dcf,y);
            case 'nufft_cart'
        end     
    case 'rufis'
    case 'cart'
    case 'sim'
end
par.oper.D_adj = D_adj;    

%% Fourier + temporal subpsace
switch par.recon.psd
    case 'qti'
        switch par.recon.sp_method
            case 'grid'
            case 'grid_cart'
            case 'nufft'
                F_for = @(x) fwd_nufft(x,par);
                F_adj = @(y) adj_nufft(y,par);
            case 'nufft_cart'
        end     
    case 'rufis'
    case 'cart'
    case 'sim'
 end

%% Sensitivities - not implemented in this version
S_for = @(x) x;
S_adj = @(c) c;

%% Dictionary matching
if par.f.match
    if par.f.compress && size(dict.D,2) > par.ind.Tv*par.ind.P
       V = par.V(:,1:par.ind.Tv);
       dict.D = dict.D*V; 
    end
    
    D_patch =   @(d) get_image_patches(d,par); % Reformat data for fitting
    D_mask  =   @(d) get_masked_voxels(d,par); % Mask data to fit
    D_tm    =   @(d) qti_dtm_fit(dict,d,par); % Match patch-based image matrix to the dictionary
    D_out   =   @(d) update_out_data(d,par); % Reshape, update, and out
    par.oper.D_match =   @(d) D_out(D_tm(D_mask(D_patch(d)))); %full dictionary matching operator
end

%% Full forward/adjoint operators
par.oper.A_for = @(a) P_for(F_for(S_for(a)));
par.oper.A_adj = @(y) S_adj(F_adj(D_adj(P_for(y))));
par.oper.AHA =   @(a) S_adj(F_adj(D_adj(P_for(F_for(S_for(a))))));

end