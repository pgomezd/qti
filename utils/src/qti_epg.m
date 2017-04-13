function dict = qti_epg(T1,T2,Diff,seq,par)
%% Simulate qti signal evolutions with EPG
% =========================================================================
% D = qti_epg(T1,T2,seq,par)
% %
%   Input
%       T1                  [L x 1]     T1 relaxation [s]
%       T2                  [L x 1]     T2 relaxation [s]
%       Diff                [L x 6]     Diffusion tensor components (Dii, Dij) [m^2/s] 
%       seq     Sequence details
%                   .FA     [T x 1] flip angles [deg]
%                   .TR     [T x 1] repetition times [s]
%                   .TE     [T x 1] echo times [s]
%               Optional sequence details
%                   .Tg1    [T x 1] Time of first diff encoding gradient (pre-readout) [s]
%                   .T2_Tg1 [T x 1] Time to Tg1 [s]
%                   .kg1    [T x 3] k-space traversal due to Tg1 [rad/m in x,y,z]
%                   .Tg2    [T x 1] Time of second diff encoding gradient (post-readout) [s]
%                   .T2_Tg2 [T x 1] Time to Tg2 
%                   .kg2    [T x 3] k-space traversal due to Tg2 [rad/m in x,y,z]
%                   .TI             Duration of inversion pulse (if par.seq.inversion)
%                   .inversion      Inversion pulse flag
%                   .phase_cycle    Phase cycling flag
%   Output
%       dict
%           .D_nn       [L combinations x T timepoints]     Unnormalized dictionary 
%           .D          [L combinations x T timepoints]     Normalized dictionary of simulated signals
%           .V          [T timepoints x V singular values]  Singular vectors  
%           .lut        [L timepoints x Q parameters]       Parameter look-up table
% =========================================================================
% Log
% v1.1: 24.02.17
% v1.2: 03.03.17: added diffusion encoding gradients and diffusion tensor
% v2.1: 13.04.17: removed tensor model, only T1/T2 relaxation
% =========================================================================
% (c) P. Gomez
% =========================================================================
%% Defaults and set-up
if ~isfield(par,'f'); par.f.f = true; end;  
if ~isfield(par.f,'normalize_dict'); par.f.normalize_dict = true; end
if ~isfield(par.f,'get_svd'); par.f.get_svd = true; end
if ~isfield(par.f,'use_parallel'); par.f.use_parallel = false; end; 
if ~isfield(par.f,'save_bv'); par.f.save_bv = false; end; 
if ~isfield(par,'ind'); par.ind.f = true; end;  
if ~isfield(par.ind,'Tv'); par.ind.Tv = 20; end;  

%% Review variables
if ~isfield(seq,'FA'); error('seq.FA needed as an input'); end
if ~isfield(seq,'TR'); error('seq.TR needed as an input'); end
if ~isfield(seq,'TE'); error('seq.TE needed as an input'); end
if ~isfield(seq,'inversion');   seq.inversion = true;  end
if (seq.inversion && ~isfield(seq,'TI'))
    error('seq.TI needed as input if seq.inversion = true');
end
if ~isfield(seq,'phase_cycle'); seq.phase_cycle = false; end
if ~isfield(seq,'order'); seq.order = 8; end %number of higher order components to get rid off
if numel(seq.TR)==1
    seq.TR = repmat(seq.TR,[1 length(seq.FA)]);
end
if numel(seq.TE)==1
    seq.TE = repmat(seq.TE,[1 length(seq.FA)]);
end

%% Diffusion gradients
if ~isfield(seq,'Tg1');             seq.Tg1 = 0;                      end %default is no gradient
if ~isfield(seq,'T2_Tg1');          seq.T2_Tg1 = 0;                   end 
if ~isfield(seq,'kg1');             seq.kg1 = 0;                      end 
if ~isfield(seq,'Tg2');             seq.Tg2 = 0;                      end
if ~isfield(seq,'T2_Tg2');          seq.T2_Tg2 = 0;                   end 
if ~isfield(seq,'kg2');             seq.kg2 = 0;                      end 
if ~isfield(seq,'accumulate_bvals');   seq.accumulate_bvals = false;  end 

if numel(seq.Tg1)==1;   seq.Tg1 = repmat(seq.Tg1,[length(seq.FA) 1]);       end
if numel(seq.T2_Tg1)==1;seq.T2_Tg1 = repmat(seq.T2_Tg1,[length(seq.FA) 1]); end
if numel(seq.Tg2)==1;   seq.Tg2 = repmat(seq.Tg2,[length(seq.FA) 1]);       end
if numel(seq.T2_Tg2)==1;seq.T2_Tg2 = repmat(seq.T2_Tg2,[length(seq.FA) 1]); end
if size(Diff,1)==1;      Diff = repmat(Diff,[length(T1) 1]);                end
if size(seq.kg1,1)==1;   seq.kg1 = repmat(seq.kg1,[length(seq.FA) 1]);      end
if size(seq.kg2,1)==1;   seq.kg2 = repmat(seq.kg2,[length(seq.FA) 1]);      end

%% Review tensor
if max(size(seq.kg1,2),size(seq.kg2,2)) > 3
    error('Directional gradient must have maximum 3 components');
elseif max(size(seq.kg1,2),size(seq.kg2,2)) > 1 && size(Diff,2) ~=6
   error('Directional gradients, D should be a vector with six components'); 
elseif  max(size(seq.kg1,2),size(seq.kg2,2)) == 1 && size(Diff,2) == 6
    error('Non-directional gradients, D should be a scalar')
end

%% No tensor model
if size(Diff,2) > 1
    error('Tensor model not available in this version');
end

%% Simulate signal
T = length(seq.FA);
L = length(T1);
D = zeros(L,T,'single');
if par.f.use_parallel
    parfor l=1:L
        D(l,:) = qti_epg_single(T1(l),T2(l),squeeze(Diff(l,:)).',seq);
    end
else
    for l=1:L
        D(l,:) = qti_epg_single(T1(l),T2(l),squeeze(Diff(l,:)).',seq);
    end
end

%% Get singular values
if par.f.get_svd 
    [~,~,V]=svd(D,0); 
    dict.V=V(:,1:par.ind.Tv);
end

%% Normalize
D_nn = D; %save no norm D for subspace projection;    
if par.f.normalize_dict
    normD = zeros(L,1,'single');
    if par.f.use_parallel
        parfor l = 1:size(D,1)
            normD(l) = norm(D(l,:));
            D(l,:)=D(l,:)/normD(l);
        end
    else
        for l = 1:size(D,1)
            normD(l) = norm(D(l,:));
            D(l,:)=D(l,:)/normD(l);
        end
    end
    dict.D = D;
    dict.normD = normD;
end
      
%% Output    
dict.D_nn = D_nn;
if sum(sum(Diff)) == 0
    dict.lut  = ([T1(:),T2(:)]); 
else
    dict.lut  = ([T1(:),T2(:),Diff]);
end
end
