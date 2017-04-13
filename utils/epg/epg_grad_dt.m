%
%function [FpFmZ] = epg_grad(FpFmZ,noadd)
%	Propagate EPG states through a "unit" gradient.
%	
%	INPUT:
%		FpFmZ = 3xN vector of F+, F- and Z states.
%		noadd = 1 to NOT add any higher-order states - assume
%			that they just go to zero.  Be careful - this
%			speeds up simulations, but may compromise accuracy!
%
%	OUTPUT:
%		Updated FpFmZ state.
%
%       SEE ALSO:
%               epg_grad, epg_grelax
%
%       B.Hargreaves.
% =========================================================================
% Modified by P. Gomez to shift the bvals as well
% =========================================================================
% Log
% v1.1: 10.03.17
% =========================================================================
% P. Gomez
% =========================================================================

function [FpFmZ,bvals] = epg_grad_dt(FpFmZ,noadd,bvals)

if (nargin < 2) noadd=0; end;	% Add by default.  

% Gradient does not affect the Z states.

if (noadd==0)
  FpFmZ = [FpFmZ [0;0;0]];	% Add higher dephased state.
  bvals(end+1,:,:) = 0; %add higher dephased state
end;

FpFmZ(1,:) = circshift(FpFmZ(1,:),[0 1]);	% Shift Fp states.
FpFmZ(2,:) = circshift(FpFmZ(2,:),[0 -1]);	% Shift Fm states.
FpFmZ(2,end)=0;					% Zero highest Fm state.
FpFmZ(1,1) = conj(FpFmZ(2,1));			% Fill in lowest Fp state.


bvals(:,:,1) = circshift(bvals(:,:,1),[1 0]); %F+
bvals(:,:,2) = circshift(bvals(:,:,2),[-1 0]); %F-
bvals(1,:,1) = 0; %first F+ state to 0
bvals(end,:,2) = 0; %highest F- state to 0
