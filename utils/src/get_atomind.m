function atomind = get_atomind(mask,amask,w_dist,L,Nd,datadims,S)
%% atomind: get atom indexes that are in mask and search window
% =========================================================================
% Finds atoms inside search window for an atlas-based dictionary
% =========================================================================
% v1.1: 20.02.16
% =========================================================================
% P. Gomez
% =========================================================================

%Extend padded mask for all subjects
atomind = zeros(L,Nd); %L x Nd
if length(datadims)==3 %2D + Q data
    mask_sub  = repmat(mask,[1 1 S]);
elseif length(datadims)==4 %3D + Q data
    mask_sub  = repmat(mask,[1 1 1 S]);
end
mask_ind = find(mask_sub==1);
if length(datadims)==3
    ws1 = w_dist(1);
    ws2 = w_dist(2);
elseif length(datadims)==4
    ws1 = w_dist(1);
    ws2 = w_dist(2);
    ws3 = w_dist(3);
end

for n=1:Nd
   if length(datadims)==3 %2D + Q data
        [si, sj] = ind2sub(size(mask_sub),mask_ind(n)); %subindex of mask
        si = si-ws1:si+ws1;
        si = si(si>0); %select only positive indexes
        si = si(si<=datadims(1)); %select only indexes within image
        sj = sj-ws2:sj+ws2;
        sj = sj(sj>0); %select only positive indexes
        sj = sj(sj<=datadims(2)); %select only indexes within image
        ss = 1:S; %apply window on all subjects
        window_mask =zeros(size(mask_sub));
        window_mask(si,sj,ss) = 1;
   elseif length(datadims)==4 %3D + Q data
        [si,sj,sk] = ind2sub(size(mask_sub),mask_ind(n)); %subindex of mask
        si = si-ws1:si+ws1;
        si = si(si>0); %select only positive indexes
        si = si(si<=datadims(1)); %select only indexes within image
        sj = sj-ws2:sj+ws2;
        sj = sj(sj>0); %select only positive indexes
        sj = sj(sj<=datadims(2)); %select only indexes within image
        sk = sk-ws3:sk+ws3;
        sk = sk(sk>0); %select only positive indexes
        sk = sk(sk<=datadims(3)); %select only indexes within image
        ss = 1:S; %apply window on all subjects
        window_mask = zeros(size(mask_sub));
        window_mask(si,sj,sk,ss) = 1;
   end           
   %Find entries that are: in search window and in atlas mask
   window_mask = window_mask(amask>0);
   atomind(:,n) = window_mask(:);
end
            
end