function y = ctr_ftn(x,oper)
%ctr_ftn
%centered n-D forward/inverse Fourier transform
%
%Syntax
%y = ctr_ftn(x)
%y = ctr_ftn(x,oper)
%
%Description
%y = ctr_ftn(x) assumes the origin is in the center of x and y axes and
%perform forward Fourier transform for all dimensions of x.
%y = ctr_ftn(x,oper) uses the 1D array oper to specify the actual operation
%on x. For example, 'finf' represents forward FT in the first dimension,
%inverse FT in the second, no operation in the third, and forward FT again
%in the fourth dimension. If the length oper is smaller than the number of
%dimensions of x, then forward FT will be the default operations for the
%last a few dimensions.
%
%Remarks
%
%Examples
%
%Algorithm
%
%See Also
%
%References
%
%Written By 
%Dan Xu
%Department of Electrical and Computer Engineering
%University of Illinois at Urbana-Champaign
%October 11, 2005
%Email: danxu@uiuc.edu
%
%Version 1.0

%% Set defaults
xdim = ndims(x);

if nargin < 2 | ~isvector(oper)
    oper = repmat('f',[1,xdim]); %default operation is all forward FT's
end

if length(oper) > xdim
    oper = oper(1:xdim); %ignore the operations that exceed the number of dimensions
end

for ind = 1:xdim
    if ind > length(oper)
        new_oper(ind) = 'f'; %set unspecified digit to forward FT
    else
        if oper(ind) == 'f' | oper(ind) == 'i' | oper(ind) == 'n'
            new_oper(ind) = oper(ind); %inherit the recognizable digit
        else
            new_oper(ind) = 'f'; %set unrecognized digit to forward FT
        end
    end
end

oper = new_oper;

clear new_oper;


%% Forward and inverse Fourier transforms
for ind = 1:xdim
    switch oper(ind)
        case 'f'
            x = fftshift(fft(ifftshift(x,ind),[],ind),ind);
        case 'i'
            x = fftshift(ifft(ifftshift(x,ind),[],ind),ind);
        case 'n'
        otherwise
            error('A situation that is supposed not to happen!');
    end
end

y = x;
