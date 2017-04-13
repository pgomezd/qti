%% Demo 3: Create raw data with trjaectory waveform
% =========================================================================
addpath(genpath('utils'))
%% Load phantom and select slice
im_vol = loadminc(['data' filesep 'phantom.mnc']); % volume is  217   181   181 
mtx_size = 256;
img = zeros(mtx_size,mtx_size);
slice_ind = 90;
img(21:237,37:217) = im_vol(:,:,slice_ind);
img(img==8) = 3; %set glial matter to WM
img = img.*(img<4);  % keep only WM/GM/CSF
img = img(end:-1:1,:);
%% create reference qmaps
% Tissue        T1 (ms)      T2 (ms)     phantom index
% Background    -           -           0
% CSF           4000        1500        1
% GM            1700        95          2
% WM            685         65          3    
Q = 2;
T1 = [0;4000;1700;685]/1e3;
T2 = [0;1500;95;65]/1e3;
PD = [0; 100;90;85];

T1_map = img;
T2_map = img;
PD_map = img;

for l = 0:length(T1)-1
    T1_map(img==l) = T1(l+1);
    T2_map(img==l) = T2(l+1);
    PD_map(img==l) = PD(l+1);
end

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

%% Simulate signals and add PD scaling
dict = qti_epg(T1,T2,0,seq,par);
dict.D_nn = bsxfun(@times,dict.D_nn,PD);
% plot(abs(dict.D_nn).')

%% Assign simulated signals to voxels
ref.X = zeros(Nreps,mtx_size,mtx_size);
for l = 0:length(T1)-1
    ref.X(:,img==l) = repmat(dict.D_nn(l+1,:).',[1 numel(img(img==l))]);
end
ref.X = permute(ref.X,[2 3 1]);

%% Assign qmap to ref and save
ref.qmap(:,:,1) = single(T1_map);
ref.qmap(:,:,2) = single(T2_map);
ref.pd = single(PD_map);
ref.msk = ref.pd > 0;
save(['data' filesep 'ref.mat'],'ref');

%% setup par
par.f.apply_tempsubspace = false;
par.f.compress = false;
par.ind.mtx_reco = 256;
par.recon.Jk = [2 2];
par.ind.Nslices = 1;
par.ind.Ncoils = 1;
par.f.match = true;
par.ind.NTimepoints = Nreps;
par = setup_qti_par(dict,ref,par);

%% Load acquisition trajectory
wf = load(['data' filesep 'wave.mat']);
k = wf.k;
ind = wf.ind;
tmp = zeros(1,size(k,2),2);
tmp(1,:,1) = real(k); tmp(1,:,2) = imag(k);
k = tmp;
par.ind.dsamp = sum(ind(1,:));
par.ind.nint = size(ind,1);

%% Get nufft operator
ck = k(1,:,1) + 1i*k(1,:,2); %get complex k-space
par.recon.Gn = NUFFT(ck,1,1,[0 0],[par.ind.mtx_reco par.ind.mtx_reco],2,par.recon.Jk);  %get nufft operator, no dcf: NUFFT(ck,1,...), dcf: NUFFT(ck,wf.dcf,...)             

%% Go to k-space with raw data
data_in = zeros(par.ind.mtx_reco,par.ind.mtx_reco,par.ind.Ncoils,par.ind.NTimepoints);
data_in(:,:,1,:) = ref.X;
raw = fwd_nufft(data_in,par);

%% reshape raw and save
raw = single(permute(raw,[4,3,1,2]));
raw_out = zeros(size(raw,1),size(ind,2));
raw_out(:,ind(1,:)>0) = raw;
raw = raw_out;
save(['data' filesep 'raw.mat'],'raw');

%% view examples
close all
show_ind = [1 50 100];
figure;
set(gcf,'units','normalized','outerposition',[0 0.4 0.6 0.6]);

for t=1:3
    subplot(2,3,t)
    imagesc(squeeze(abs(ref.X(:,:,show_ind(t)))))
    axis image
    axis off
    colormap gray
    title(['T = ' num2str(show_ind(t))]);
end
% 
subplot(2,3,4)
imagesc(T1_map*1e3)
axis image
axis off
ax = gca;
colormap(ax,hot);
cb = colorbar(ax,'south');
cb.Position(1) = cb.Position(1)+0.01; %move
cb.Position(3) = cb.Position(3)-0.02; %reduce height
cb.Position(4) = cb.Position(4)/4; %reduce width
cb.Color = 'w';
cb.AxisLocation = 'in';
caxis([0 4000])
title('T1 (ms)')
    
subplot(2,3,5)
imagesc(T2_map*1e3)
axis image
axis off
ax = gca;
colormap(ax,hot);
cb = colorbar(ax,'south');
cb.Position(1) = cb.Position(1)+0.01; %move
cb.Position(3) = cb.Position(3)-0.02; %reduce height
cb.Position(4) = cb.Position(4)/4; %reduce width
cb.Color = 'w';
cb.AxisLocation = 'in';
caxis([0 300])
title('T2 (ms)')


subplot(2,3,6)
imagesc(PD_map)
axis image
axis off
ax = gca;
colormap(ax,hot);
cb = colorbar(ax,'south');
cb.Position(1) = cb.Position(1)+0.01; %move
cb.Position(3) = cb.Position(3)-0.02; %reduce height
cb.Position(4) = cb.Position(4)/4; %reduce width
cb.Color = 'w';
cb.AxisLocation = 'in';
caxis([70 100])
title('PD (au)')


%% view raw data
figure;
set(gcf,'units','normalized','outerposition',[0.6 0.4 0.4 0.4]);
plot(abs(raw(:,3:10)))
title('k-space components')


