function PSNR_out = PSNR(seq1,seq2)
% PSNR - AUXILIARY FUNCTION
%  Calculate the PSNR between two *NORMALISED* sequences.
%
%  Inputs:
%   [seq1, seq2] : 3D matrices representing two sequences of the same size. They
%                  need to be normalised.


%  Jose Caballero
%  Biomedical and Image Analysis Group
%  Department of Computing
%  Imperial College London, London SW7 2AZ, UK
%  jose.caballero06@imperial.ac.uk
%
%  October 2012

MSE = mean(abs(seq2(:)-seq1(:)).^2);
PSNR_out = 10*log10(1/MSE);

end