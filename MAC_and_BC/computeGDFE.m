% function [snrGDFEu, GU, WU, S0, MSWMFU] = computeGDFE(H, A, cb, Nx)
%--------------------------------------------------------------------------
% Inputs: H, A , Nx
% Outputs: snrGDFEu, Gub, Wub, S0, MSWMFU
%
% H: noise-whitened channel
% A: any square root of input autocorrelation matrix
%        This can be generalized non-square square-root
% cb: =1 if H is complex baseband; =2 if H is real baseband
% Nx: optional input of Nx when not equal to right-size of A
%
% G: feedback matrix
% GU: unbiased feedback matrix
% W: feedfoward linear equalizer
% WU: unbiased feedforward linear equalizer
% S0: sub-channel channel gains
% MSWMFU: unbiased mean-squared whitened matched filter
% snrGDFEu - unbiased SNR in dB; assumes size of R_sqrt input
%    the user should recompute SNR if there is a cyclic prefix
%
% Note: The R_sqrt need not be square non-singular, as long as
% R_sqrt*R_sqrt' = input autocorrelation matrix Rxx. 
%
%--------------------------------------------------------------------------
function [snrGDFEu, GU, WU, S0, MSWMFU] = computeGDFE(H, A, cb, Nx)
%-------------------------------------------------------------------------
% Computing Ht: Ht = H*A
Ht = H*A; % Ly-by-nx
%--------------------------------------------------------------------------

Ng = size(A, 2); % can be U*Lx if called by mu_bc
if (nargin<4)
  Nx = Ng;
end

     
% Computing Rf, Rbinv, Gbar, all Ng-by-Ng
Rf = Ht' * Ht; 
Rbinv = Rf + eye(Ng);
[G,S0] = ldl(Rbinv, 'upper');
invS0unb = pinv(S0-eye(Ng));

% Computing the matrices of interest
W = inv(G'*S0);

GU = eye(Ng) + S0*invS0unb*(G-eye(Ng));
WU = invS0unb/G';
MSWMFU = WU*Ht';
% Section to handle numerical is with large-dim determinants/geo-ave
if Ng < 10
   snrGDFEu = 10*log10(det(S0)^(1/(cb*Nx))-1);
else
   snrGDFEu = 10*log10(exp(sum(log(diag(S0)))/(cb*Nx)) - 1);
end

