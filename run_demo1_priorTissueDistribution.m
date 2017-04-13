%% Demo 1: Prior distribution
% =========================================================================
%% get grid and priors
T1_step = 0.01*2; %in s
T2_step = 0.005*2; %in s
T1_r = 0.01:T1_step:6;
T2_r = 0.005:T2_step:3;
[T1_grd,T2_grd] = ndgrid(T1_r, T2_r);
parameter_space = [T2_grd(:),T1_grd(:)];
mu = [ ...
  0.095 1.7; ... %GM 
  0.065 .685 ; ... %WM 
  1.5   4.0; ... % CSF 
  0.275 1.9 ; ...% BV 
 ];

sigma = cat(3,...
 [0.010 0;0 0.200], ... % GM STD
 [0.010 0;0 0.090], ... % WM STD
 [0.200 0;0 0.300], ... % CSF STD
 [0.050 0;0 0.200] ... % BV STD
 );   

%% simualate 
sigma = sigma*1.2; %increase prior by 20% to cover more space
prior = zeros([size(T1_grd),size(mu,1)]);
for p = 1:size(mu,1)
    cp = reshape(mvnpdf(parameter_space,mu(p,:),sigma(:,:,p)),size(T2_grd));
    prior(:,:,p) = cp/sum(cp(:));
end
prior = (1/size(prior,3))*sum(prior,3);

%% plot image
close all;
imagesc(T2_grd(:)*1e3,T1_grd(:)*1e3,prior)
axis xy
ylabel('T1 (ms)')
xlabel('T2 (ms)')


