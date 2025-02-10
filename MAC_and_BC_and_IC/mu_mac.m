% function [Bu, GU, WU, S0, MSWMFU] = mu_mac(H, AU, Usize , cb)
%--------------------------------------------------------------------------
% Inputs: Hu, A , Uind
% Outputs: snrGDFEu, Gub, Wub, S0, MSWMFU
%
% H: noise-whitened channel matrix [HU ... H1]
% AU: Block Diag square root discrete modulators, blkdiag([AU ... A1])
% Usize: # of dimensions for each user U ... 1 in 1 x U row vector
% cb: = 1 if complex baseband or 2 if real baseband channel
%
% G: feedback matrix
% GU: unbiased feedback matrix
% W: feedfoward linear equalizer
% WU: unbiased feedforward linear equalizer
% S0: sub-channel channel gains
% MSWMFU: unbiased mean-squared whitened matched filter
% Bu - user u's bits/symbol 
%    the user should recompute SNR if there is a cyclic prefix
% 
%
%--------------------------------------------------------------------------
function [Bu, GU, WU, S0, MSWMFU] = mu_mac(H, A, Usize, cb)

%
[dum , U] = size(Usize);
Bu=zeros(1,U);

% Computing Ht: Ht = H*A
Ht = H*A;
%--------------------------------------------------------------------------

     
% Computing Rf, Rbinv, Gbar
Rf = Ht' * Ht;
Rbinv = Rf + eye(size(Rf));
Gbar = chol(Rbinv);

% Computing the matrices of interest
G = inv(diag(diag(Gbar)))*Gbar;
S0 = diag(diag(Gbar))*diag(diag(Gbar));
W = inv(S0)*inv(G');

GU = eye(size(G)) + S0*inv(S0-eye(size(G)))*(G-eye(size(G)));
WU = inv(S0-eye(size(G)))*inv(G');
MSWMFU = WU*Ht';
index=0;
for u=1:U
  for l=1:Usize(u)
  Bu(u) = Bu(u)+ (1/cb)*log2(S0(index+l,index+l));
  end
  index=index+Usize(u);
end
Bu=Bu;

