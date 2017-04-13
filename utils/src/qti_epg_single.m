function S = qtip_epg_single(T1,T2,Diff,seq)
%% set-up
F = [0;0;1];	  % Equilibrium Magnetization.	
T = length(seq.FA);
S = zeros(T,1);
if seq.order == 0
    seq.order = T; %as many orders as repetitions
end
F(1,seq.order)=0; % Allocate states
noadd = 1;	% Don't add higher-order states

  
%% inversion pulse
if seq.inversion
    F = epg_rf(F,pi,0); % RF Rotation: configuration state/flip angle/phase cycling
    F = epg_grelax(F,T1,T2,seq.TI,0,Diff,0,noadd);	%T1,T2/time/diff grad/diff/gradon/add
end

%% main sequence
for t=1:T
    % RF Rotation: configuration state/flip angle/phase cycling
    if seq.phase_cycle
      F = epg_rf(F,seq.FA(t)*pi/180,mod(t+1,2)); 
    else
      F = epg_rf(F,seq.FA(t)*pi/180,0);    
    end
    % First block Relax-Gradient-Relax
    if seq.Tg1(t)>0
        F = epg_grelax(F,T1,T2,seq.T2_Tg1(t),seq.kg1(t,:),Diff,0,noadd);	%T1,T2/time/diff grad/diff/gradon/add
        F = epg_grelax(F,T1,T2,seq.Tg1(t),seq.kg1(t,:),Diff,1,noadd);
    end
    F = epg_grelax(F,T1,T2,seq.TE(t)-seq.Tg1(t)-seq.T2_Tg1(t),seq.kg1(t,:),Diff,0,noadd);
   
    % Record signal
    S(t) = F(1,1);		
    
    % Second block Relax-Gradient-Relax
    if seq.Tg2(t)>0
        F = epg_grelax(F,T1,T2,seq.T2_Tg2(t),seq.kg2(t,:),Diff,0,noadd);	
        F = epg_grelax(F,T1,T2,seq.Tg2(t),seq.kg2(t,:),Diff,1,noadd);
    end
        F = epg_grelax(F,T1,T2,seq.TR(t)-seq.TE(t)-seq.Tg2(t)-seq.T2_Tg2(t),0,Diff,0,noadd); 
end
end