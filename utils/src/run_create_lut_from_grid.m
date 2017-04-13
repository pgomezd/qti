%% Create look up table (lut) from grid
% Script that creates combinations of biological parameters that match the
% simulated FISP 
% =========================================================================
% Script 
%   1. Uses nd grid to create bpar.lut: [L entries, Q parameters]
% =========================================================================
% v1.1: 19.10.15 
% =========================================================================
% P. Gomez
% GE Global Research
% =========================================================================

%Create par.lut to make fit easier
if par.ind.Q == 3
   [insout1, insout2, insout3] = ndgrid(par.bp.T1,par.bp.T2,par.bp.B1);
   bpar.lut=zeros(numel(insout1),par.ind.Q);
   par.lut(:,1) =  insout1(:);
   par.lut(:,2) =  insout2(:);
   par.lut(:,3) =  insout3(:);
  
else
    par.lut = zeros(size(D,1),par.ind.Q);
    clear insout;
    clear insin;
    insout = '[';
    insin = ' = ndgrid(';

    for q=1:par.ind.Q
         insout = [insout,'out',num2str(q),','];
         insin = [insin,'par.bp.(par.str.bp{',num2str(q),'}),'];
    end
    insout(end)=']';
    insin(end)=')';
    insin = [insin,';'];

    eval([insout,insin]); %evaluate instruction nd grid for all Q parameters
    for q=1:par.ind.Q
       par.lut(:,q)= eval(['out',num2str(q),'(:)']); %vectorized version of each out grid
       eval(['clear out',num2str(q)]); %clear out vector
    end 
end
clear ins*;