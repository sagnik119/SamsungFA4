% function [Bu, Gunb, S0, MSWMFunb] = mu_bc(H, AU, Usize , cb)
%--------------------------------------------------------------------------
% Inputs: Hu, AU , Usize, cb
% Outputs: Bu, Gunb, Wunb, S0, MSWMFunb
%
% H: noise-whitened BC matrix [H1 ; ... ; HU] (with actual noise, not wcn)
%    sum-Ly x Lx x N, assume already normalize by sqrt(N)
% AU: Block-row square-root discrete modulators, [A1 ... AU]
%     Lx  x (U * Lx) x N
% Usize: # of (output, Lyu) dimensions for each user U ... 1 in 1 x U row vector
% cb: = 1 if complex baseband or 2 if real baseband channel
%
% GU: unbiased precoder matrices
% S0: sub-channel dimensional channel SNRs
% MSWMFunb: users' unbiased diagonal mean-squared whitened matched matrices
% Bu - users bits/symbol 1 x U
%    the user should recompute SNR if there is a cyclic prefix
% 
%
%--------------------------------------------------------------------------
function [Bu, Gunb, S0, MSWMFunb] = mu_bc(H, AU, Usize, cb)

%
[~ , U] = size(Usize);
[~,Lx,N]=size(H);
idx_end = cumsum(Usize);
idx_start = [1, idx_end(1:end-1)+1];

% Computing Ht: Ht = H*A

%--------------------------------------------------------------------------
S0=cell(U,N);
% b=cell(U,N);
for u=1:U
    for n=1:N
     S0{u,n} = zeros(Lx,Lx);
%     b{u,n}=0;
    end
end
MSWMFunb=cell(U,N);
Gunb=cell(U,N);
B = zeros(U,N);

% Use compute GDFE

for u=1:U
  for n=1:N
   [~,Gunbtemp,~,S0temp,MSWMFtemp]= computeGDFE( ...
    H(idx_start(u):idx_end(u),:,n), AU(:,:,n),cb);
   Gunb{u,n}=Gunbtemp(Lx*(u-1)+(1:Lx),:);
   S0{u,n}=S0temp(Lx*(u-1)+(1:Lx), Lx*(u-1)+(1:Lx));
   MSWMFunb{u,n}=MSWMFtemp(Lx*(u-1)+(1:Lx),:);
   B(u,n) = sum(log2(diag(S0{u,n})))/cb;
  end
end
Bu = sum(B,2);






