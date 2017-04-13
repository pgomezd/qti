%
%	Very simple simulation of Gradient-Spoiled 
%	sequence using EPG methods.
%

N = 200;	% Number of Sequence repetitions.
TR = .01;	% 10ms

F = [0;0;1];	% Equilibrium Magnetization.
		% [F+; F-; Z],  all longitudinal in Z0 state.


S1 = zeros(N,1);
S2 = zeros(N,1);

for k=1:N
  F = epg_rf(F,pi/6,0);		% RF Rotation.
  S1(k) = F(1,1);		% Record "Signal"
  F = epg_grelax(F,1,.1,TR,1,0,1,0);	% Relaxation, spoiling
  S2(k) = F(1,1);		% Record "Signal"
end;
 

plot([abs(S1) abs(S2)]);


 
