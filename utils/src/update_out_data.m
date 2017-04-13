function  out = update_out_data(data,par)
%% Update out data
% =========================================================================
% Re-arranges data to output struct
% =========================================================================
% v1.1: 15.08.16
% v2.1: 30.03.17 - added 2x2 patches. Tag: pgd.2x2
% =========================================================================
% P. Gomez
% =========================================================================

%% Get variables
datadims = par.ind.datadims;
P = par.ind.P;
Q = par.ind.Q;
T = par.ind.NTimepoints;
ix = par.ind.ix;
iy = par.ind.iy;
if length(datadims)==4
    iz = par.ind.Nslices;
else
    iz = 1;
end
p_side = par.ind.p_side;
% --- 30.03.2017 pgd.2x2 --- %
if P==4
   par.ind.N = (mod(ix,2)+ix)*(mod(iy,2)+iy)/P*iz;
end
N = par.ind.N;

% --- 30.03.2017 pgd.2x2 --- %

if(par.f.Xout)
    Xout = zeros(N,P,T,'single');
    Xout(data.mask>0,:) = single(data.Xfit);
end

if(par.f.qout) && isfield(data,'qmapfit')
    qmap = zeros(N,P,Q,'single');
    qmap(data.mask>0,:) = single(data.qmapfit);
    qmap(isnan(qmap)) = 0; %remove possible NaN from infeasible search window
end

% --- 30.03.2017 pgd.2x2 --- %
if P == 4
   if length(datadims)==3
       Xcout = zeros([(mod(ix,2)+ix),(mod(iy,2)+iy),T]);
       if par.f.use_parallel
           parfor t = 1:T
               Xcout(:,:,t) = col2im(squeeze(Xout(:,:,t))',[p_side p_side],[(mod(ix,2)+ix) (mod(iy,2)+iy)],'distinct');
           end
       else
           for t = 1:T
               Xcout(:,:,t) = col2im(squeeze(Xout(:,:,t))',[p_side p_side],[(mod(ix,2)+ix) (mod(iy,2)+iy)],'distinct');
           end
       end
       Xout = Xcout;
       clear Xcout;
       if(par.f.qout) && isfield(data,'qmapfit')
           cqmap = zeros([(mod(ix,2)+ix),(mod(iy,2)+iy),Q]);
           if par.f.use_parallel
               parfor q = 1:Q
                   cqmap(:,:,q) = col2im(squeeze(qmap(:,:,q))',[p_side p_side],[(mod(ix,2)+ix) (mod(iy,2)+iy)],'distinct');
               end
           else
               for q = 1:Q
                   cqmap(:,:,q) = col2im(squeeze(qmap(:,:,q))',[p_side p_side],[(mod(ix,2)+ix) (mod(iy,2)+iy)],'distinct');
               end
           end
           qmap = cqmap;
           clear cqmap; 
       end
   elseif length(datadims)==4
       if par.f.kernel3D
            error('Only slice-wise processing with 2x2 patches');
       else
          Xout = reshape(Xout,[(mod(ix,2)+ix)*(mod(iy,2)+iy)/P,iz,P,T]);
          Xcout = zeros([(mod(ix,2)+ix),(mod(iy,2)+iy),iz,T]);
          if par.f.use_parallel
              parfor t = 1:T
                  for z=1:iz
                    Xcout(:,:,z,t) = col2im(squeeze(Xout(:,z,:,t))',[p_side p_side],[(mod(ix,2)+ix) (mod(iy,2)+iy)],'distinct');
                  end
              end
          else
              for t = 1:T
                  for z=1:iz
                    Xcout(:,:,z,t) = col2im(squeeze(Xout(:,z,:,t))',[p_side p_side],[(mod(ix,2)+ix) (mod(iy,2)+iy)],'distinct');
                  end
              end
          end
          Xout = Xcout;
          clear Xcout; 
          if(par.f.qout) && isfield(data,'qmapfit')
              qmap = reshape(qmap,[(mod(ix,2)+ix)*(mod(iy,2)+iy)/P,iz,P,Q]);
              cqmap = zeros([(mod(ix,2)+ix),(mod(iy,2)+iy),iz,Q]);
              if par.f.use_parallel
                  parfor q = 1:Q
                      for z=1:iz
                        cqmap(:,:,z,q) = col2im(squeeze(qmap(:,z,:,q))',[p_side p_side],[(mod(ix,2)+ix) (mod(iy,2)+iy)],'distinct');
                      end
                  end
              else
                  for q = 1:Q
                      for z=1:iz
                        cqmap(:,:,z,q) = col2im(squeeze(qmap(:,z,:,q))',[p_side p_side],[(mod(ix,2)+ix) (mod(iy,2)+iy)],'distinct');
                      end
                  end
              end
              qmap = cqmap;
              clear cqmap; 
           end
       end
   end
else
% --- 30.03.2017 pgd.2x2 --- %
    if length(datadims)==3
        switch par.recon.update
            case 'voxel'
                if(par.f.Xout)
                    Xout = reshape(Xout,[datadims(1:end-1),P,T]);             
                    Xout = squeeze(Xout(:,:,(P+1)/2,:));
                end
                if(par.f.qout) && isfield(data,'qmapfit')
                    qmap = reshape(qmap,[datadims(1:end-1),P,Q]);
                    qmap = squeeze(qmap(:,:,(P+1)/2,:));
                end
            case 'patch'
                if(par.f.Xout)
                    Xout = reshape(Xout,[N,P,T]);
                    Xout = col2imn(Xout,[ix iy],p_side); 
                end
                if(par.f.qout) && isfield(data,'qmapfit')
                    qmap = reshape(qmap,[N,P,Q]);
                    qmap = col2imn(qmap,[ix iy],p_side);
                end
        end
    elseif length(datadims)==4

        switch par.recon.update
            case 'voxel'
                if(par.f.Xout)
                    Xout = reshape(Xout,[datadims(1:end-1),P,T]);
                    Xout = squeeze(Xout(:,:,:,(P+1)/2,:));
                end
                if(par.f.qout) && isfield(data,'qmapfit')
                    qmap = reshape(qmap,[datadims(1:end-1),P,Q]);
                    qmap = squeeze(qmap(:,:,:,(P+1)/2,:));
                end
            case 'patch'
                % --- 30.05.2016 pgd.3Dkernel --- %
                if par.f.kernel3D 
                    if par.f.Xout
                        Xout = reshape(Xout,[N,P,T]);
                        Xout = col2vol(Xout,[ix iy iz],p_side);
                    end

                    if(par.f.qout) && isfield(data,'qmapfit')
                        qmap = reshape(qmap,[N,P,Q]);
                        qmap = col2vol(qmap,[ix iy iz],p_side);
                    end
                else
                    if par.f.Xout
                        if P == 1
                            Xout = squeeze(reshape(Xout,[datadims(1:end-1),P,T]));                        
                        else
                            Xout = reshape(Xout,[ix*iy,iz,P,T]);
                            Xcout = zeros([ix,iy,iz,T]);
                            if par.f.use_parallel
                                parfor z = 1:iz
                                    Xcout(:,:,z,:) = col2imn(squeeze(Xout(:,z,:,:)),[ix iy],p_side);
                                end
                            else
                                for z = 1:iz
                                    Xcout(:,:,z,:) = col2imn(squeeze(Xout(:,z,:,:)),[ix iy],p_side);
                                end
                            end
                            Xout = Xcout;
                            clear Xcout;
                        end
                    end
                    if(par.f.qout) && isfield(data,'qmapfit')
                        if P == 1
                            qmap = squeeze(reshape(qmap,[datadims(1:end-1),P,Q]));                        
                        else
                            qmap = reshape(qmap,[ix*iy,iz,P,Q]);
                            cqmap = zeros([ix iy iz Q]);
                            if par.f.use_parallel
                                parfor z = 1:iz
                                    cqmap(:,:,z,:) = col2imn(squeeze(qmap(:,z,:,:)),[ix iy],p_side);
                                end
                            else
                                for z = 1:iz
                                    cqmap(:,:,z,:) = col2imn(squeeze(qmap(:,z,:,:)),[ix iy],p_side);
                                end
                            end
                            qmap = cqmap;
                            clear cqmap;
                        end
                    end
                end
                % --- 30.05.2016 pgd.3Dkernel --- %
        end
    end
end

%% Out
if par.f.Xout
    out.Xfit = single(Xout);
    if isfield(data,'X')
        out.X = data.X;            
    end
end

if par.f.qout && isfield(data,'qmapfit')
    out.qmap = single(qmap);
end

if par.f.pdout && isfield(data,'pdfit')
    out.pd = single(zeros(N,1));
    out.pd(data.mask>0) = single(data.pdfit);
    % --- 30.03.2017 pgd.2x2 --- %
    if P == 4
      pd = reshape(repmat(out.pd,[1 P]),[(mod(ix,2)+ix)*(mod(iy,2)+iy)/P,iz,P]);
      cpd = zeros([(mod(ix,2)+ix),(mod(iy,2)+iy),iz]);
      if par.f.use_parallel
          parfor z=1:iz
              cpd(:,:,z) = col2im(squeeze(pd(:,z,:))',[p_side p_side],[(mod(ix,2)+ix) (mod(iy,2)+iy)],'distinct');
          end
      else
          for z=1:iz
              cpd(:,:,z) = col2im(squeeze(pd(:,z,:))',[p_side p_side],[(mod(ix,2)+ix) (mod(iy,2)+iy)],'distinct');
          end
      end
      out.pd = squeeze(cpd);
      
    else
    % --- 30.03.2017 pgd.2x2 --- %
        out.pd =  reshape(out.pd,par.ind.datadims(1:end-1));
    end
end

if par.f.mtout && isfield(data,'mtfit')
    % --- pgd.knn --- %
    out.mt = single(zeros(N,par.ind.knn));
    out.mt(data.mask>0,:) = single(data.mtfit);
    % --- 30.03.2017 pgd.2x2 --- %
    if P == 4
      knn = par.ind.knn;  
      mt = reshape(repmat(out.mt,[1 P]),[(mod(ix,2)+ix)*(mod(iy,2)+iy)/P,iz,P,knn]);
      cmt = zeros([(mod(ix,2)+ix),(mod(iy,2)+iy),iz,knn]);
      if par.f.use_parallel
          parfor z=1:iz
              for k=1:knn
                cmt(:,:,z,k) = col2im(squeeze(mt(:,z,:,k))',[p_side p_side],[(mod(ix,2)+ix) (mod(iy,2)+iy)],'distinct');
              end
          end
      else
          for z=1:iz
              for k=1:knn
                cmt(:,:,z,k) = col2im(squeeze(mt(:,z,:,k))',[p_side p_side],[(mod(ix,2)+ix) (mod(iy,2)+iy)],'distinct');
              end
          end
      end
      out.mt = squeeze(cmt);
      
    else
    % --- 30.03.2017 pgd.2x2 --- %
        out.mt =  reshape(out.mt,[par.ind.datadims(1:end-1),par.ind.knn]);
    % --- pgd.knn --- %
    end
end

if par.f.dmout && isfield(data,'dmfit')
     % --- pgd.knn --- %
    out.dm = single(zeros(N,par.ind.knn));
    out.dm(data.mask>0,:) = single(data.dmfit);
    % --- 30.03.2017 pgd.2x2 --- %
    if P == 4
      knn = par.ind.knn;  
      dm = reshape(repmat(out.dm,[1 P]),[(mod(ix,2)+ix)*(mod(iy,2)+iy)/P,iz,P,knn]);
      cdm = zeros([(mod(ix,2)+ix),(mod(iy,2)+iy),iz,knn]);
      if par.f.use_parallel
          parfor z=1:iz
              for k=1:knn
                cdm(:,:,z,k) = col2im(squeeze(dm(:,z,:,k))',[p_side p_side],[(mod(ix,2)+ix) (mod(iy,2)+iy)],'distinct');
              end
          end
      else
          for z=1:iz
              for k=1:knn
                cdm(:,:,z,k) = col2im(squeeze(dm(:,z,:,k))',[p_side p_side],[(mod(ix,2)+ix) (mod(iy,2)+iy)],'distinct');
              end
          end
      end
      out.dm = squeeze(cdm);
      
    else
    % --- 30.03.2017 pgd.2x2 --- %
        out.dm =  reshape(out.dm,[par.ind.datadims(1:end-1),par.ind.knn]);
    % --- pgd.knn --- %
    end
end

if par.f.Yout && isfield(data,'Y')  
   out.Y = data.Y; 
end


end