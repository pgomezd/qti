%% Demo 2: Simulate transient-state signals
% =========================================================================
addpath(genpath('utils'))
%% get mean values of tissue classes
T1 = [1700,685,4000,1900]/1e3;  %GM, WM, CSF, BV
T2 = [95,65,1500,275]/1e3;


%% Set-up acquisition
clear seq;
par.f.use_parallel = false;
par.f.normalize_dict = false;
Nreps = 500;
FA_min = 1;
FA_max = 70;
seq.FA = linspace(FA_min,FA_max,Nreps);
seq.TR = 8e-3;
seq.TE = 2e-3;
seq.TI = 18e-3;
seq.inversion = true;
gamma = 4258*2*pi;	% gamma, [Rad/s/G]
G = 18*(1e-3)*(10e3); % gradient amplitude, [mT/m]*[T/mT]*[G/T] = [G/m]
seq.T2_Tg2  = 3e-3; %3 ms reading time during spiral
seq.Tg2     = 1e-3; %1 ms gradient dephasing 

%% Simulate signals
dict = qti_epg(T1,T2,0,seq,par);

%% Get acquisition time
TR = seq.TR*1e3; % in ms
qti_acq = zeros(Nreps,1);
qti_acq(1) = 20; %inversion time before (ms)
for t=1:length(qti_acq)
    if t==1
        qti_acq(t) = TR+qti_acq(t);
    else
        if t <= Nreps
            qti_acq(t) = qti_acq(t-1) + TR;
        end      
    end
end

%% plot signals
close all;
cols = {[0.6 0.6 0.6],[0 0 0],[0.1 0.4 .95],[1 0 0]}; %gray,black,blue, red
legs = {'GM','WM','CSF','BV'}; 
wdt = 1.5;
stls = {':','-','-.','--'};
ax_sig = 0.15;

for l=1:size(dict.D_nn,1)
    hold on
    plot(qti_acq/1e3,abs(dict.D_nn(l,:)),'LineWidth',wdt,'Color',cols{l},'LineStyle',stls{l})
end
axis([0 qti_acq(end)/1e3 0 ax_sig])          
ylabel('Signal (au)')
xlabel('Acquisition time (s)')
lgd = legend(legs);
lgd.Box = 'off';


