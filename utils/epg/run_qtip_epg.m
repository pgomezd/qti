%% epg_qtip
% simulate qtip signal evolution with EPG
% 24.02.17 P. Gomez

%% set-up
T1 = 0.685; % seconds
T2 = 0.065; % seconds
T = 500;	 % Number of Sequence repetitions.
TR = 8e-3;	 % 8ms
TE = 2e-3;   % 2ms
FA_min = 1;  
FA_max = 70; 
TI = 18e-3; %duration of inversion pulse before pulse
FA = linspace(FA_min,FA_max,T); %linear flip angle ramp 
F = [0;0;1];	% Equilibrium Magnetization.	
                % [F+; F-; Z],  all longitudinal in Z0 state.
invert = 1;
                
%% Simulate                
S = zeros(T,1);
phase_cycle = 0;

% inversion pulse
F = epg_rf(F,pi,0); % RF Rotation: Configuration state / Flip angle / phase cycling
F = epg_grelax(F,T1,T2,TI,0,0,0,0);	% T1,T2 relaxation / diff grad / diffusion / spoiling / add
% main sequence    
for t=1:T
    if phase_cycle
      F = epg_rf(F,FA(t)*pi/180,mod(t+1,2)); % RF Rotation: Configuration state / Flip angle / phase cycling    s
    else
      F = epg_rf(F,FA(t)*pi/180,0); % RF Rotation: Configuration state / Flip angle / phase cycling    
    end
    F = epg_grelax(F,T1,T2,TE,0,0,0,0);	% T1,T2 relaxation / diff grad / diffusion / spoiling / add
    S(t) = F(1,1);		% Record signal
    F = epg_grelax(F,T1,T2,TR-TE,1,0,1,0);	% T1,T2 relaxation / diff grad / diffusion / spoiling / add
end

%% Simulate                
S2 = zeros(T,1);
phase_cycle = 1;
% inversion pulse
F = [0;0;1];	% Equilibrium Magnetization.	
F = epg_rf(F,pi,0); % RF Rotation: Configuration state / Flip angle / phase cycling
F = epg_grelax(F,T1,T2,TI,0,0,0,0);	% T1,T2 relaxation / diff grad / diffusion / spoiling / add
% main sequence    
for t=1:T
    if phase_cycle
      F = epg_rf(F,FA(t)*pi/180,mod(t+1,2)); % RF Rotation: Configuration state / Flip angle / phase cycling    s
    else
      F = epg_rf(F,FA(t)*pi/180,0); % RF Rotation: Configuration state / Flip angle / phase cycling    
    end
    F = epg_grelax(F,T1,T2,TE,0,0,0,0);	% T1,T2 relaxation / diff grad / diffusion / spoiling / add
    S2(t) = F(1,1);		% Record signal
    F = epg_grelax(F,T1,T2,TR-TE,1,0,1,0);	% T1,T2 relaxation / diff grad / diffusion / spoiling / add
end

%% plot
close all
plot(abs(S))
hold on
plot(abs(S2))
