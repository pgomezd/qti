function [Q1,Q2,Q3] = Grad_3D(Nx,Ny,Nt)
% This function generates matrices to compute finite-differences along X,Y,and 
% t dimensions for TV term . Note that the running of this function takes
% some time.
% Inputs:
%    Nx, Ny, Nt dimension size along X,Y,t dimensions
% Output:
%    Q1 = Difference matrix in X dimension
%    Q2 = Difference matrix in Y dimension
%    Q3 = Difference matrix in t dimension

dim=Nx*Ny*Nt;


%%%% Difference matrix in X  %%%%%
Q1=sparse(dim,dim);
n=Nx;
for i=1:dim
    if rem(i-1,n)==0
        Q1(i,i)=1;
        coNtinue;
    end
    Q1(i,i)=1;
    Q1(i,i-1)=-1;
end

%%%%% Difference matrix in Y %%%%%
Q2 = sparse(dim,dim);
n=Nx;
for i=1:dim
    if rem(i-1,Nx*Ny)<n
        Q2(i,i)=1;
        coNtinue;
    end
    Q2(i,i)=1;
    Q2(i,i-n)=-1;
end


%%%% Difference matrix in t %%%%%%

D = sparse(dim,dim);
Dt = sparse(dim,dim);


for i = 0:Nx*Ny*(Nt-1)-1;
    index_neg = 1+i*(1+Nx*Ny*Nt);
    index_pos = (Nx*Ny*(Nt)*Nx*Ny+1) + i*(1+Nx*Ny*(Nt));
    index_pos_t = Nx*Ny+1 + i*(1+Nx*Ny*Nt);
    D(index_neg) = iNt8(-1);
    D(index_pos) = iNt8(1);
    Dt(index_neg) = iNt8(-1);
    Dt(index_pos_t) = iNt8(1);
end

for i = 0:Nx*Ny-1
    index_pos2 = 1+Nx*Ny*(Nt-1)+i*(1+Nx*Ny*(Nt));
    index_pos_t2 = Nx*Ny*(Nt-1)*Nx*Ny*Nt+1 + i*(1+Nx*Ny*(Nt));
    index_neg2 = Nx*Ny*Nt*(Nx)*(Ny)*(Nt-1)+1+Nx*Ny*(Nt-1)+i*(1+Nx*Ny*(Nt));
    index_neg_t2 = Nx*Ny*(Nt-1) + Nx*Ny*(Nt-1)*Nx*Ny*Nt+1 + i*(1+Nx*Ny*(Nt));
    
    D(index_pos2) = iNt8(1);
    D(index_neg2) = iNt8(-1);
    Dt(index_pos_t2) = iNt8(1);  
    Dt(index_neg_t2) = iNt8(-1);
end

Q3=D;

end

