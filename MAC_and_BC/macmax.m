% function [Rxx, bsum , bsum_lin] = SWF(Eu, h, user_ind, Rnn, N , cb)
%
% the inputs are: 
% Eu  the U x 1 energy vector. If a single scalar, same energy for all users
%     Each user energy should be scaled by N/(N+nu)if there is cyclic prefix 
%     This energy is the trace of the corresponding user Rxx (u)
% h   The TIME-DOMAIN Ly x sum(Lx(u)) x N channel for all users
% user_ind  The start index for each user, in the same order as Eu
% Rnn The Ly x Ly noise autocorrelation
% N   The number of used tones (equally spaced over (0,1/T) at N/T.
% cb  cb = 1 for complex, cb=2 for real baseband
%
% the outputs are:
% Rxx A block-diagonal psd matrix with the input autocorrelation for each
%     user on each tone. Rxx has size (sum(Lx(u)) x sum(Lx(u)) x N .
%     sum trace(Rxx) over tones and spatial dimensions equal the Eu 
% bsum the maximum rate sum.
% bsum  bsum_lin - the maximum sum rate with a linear receiver
% 
%     b is an internal convergence (vector, rms) value, but not sum rate

function [Rxx, bsum, bsum_lin] = macmax(Esum, h, Lxu, N , cb)

H = fft(h, N, 3); % scaled so that N factor no longer needed in rate calculation
[Ly, Lx, ~] = size(H);
if numel(Lxu) > 1
    U = numel(Lxu);
    if sum(Lxu)~=Lx
        error('mismatch between sum of Lxu and Lx');
    end
elseif (Lx/Lxu)~= floor(Lx/Lxu)
    error('invalid Lxu');
else
    U = Lx/Lxu;
    Lxu = Lxu*ones(1,U);
end

idx_end = cumsum(Lxu);
idx_start = [1, idx_end(1:end-1)+1];
idx_exp_end = N*idx_end;
idx_exp_start = [1, idx_exp_end(1:end-1)+1];

Hcell=mat2cell(H, Ly, Lxu, ones(1,N));
Hexpand = zeros(N*Ly,N*Lx);
for u = 1:U
    Hexpand(:,idx_exp_start(u):idx_exp_end(u)) = blkdiag(Hcell{1,u,:});
end

Lxu_max = max(Lxu);
cvx_begin quiet
cvx_solver mosek
    variable Rxxun(Lxu_max,Lxu_max,N,U) hermitian semidefinite
    Rxxun_cell = reshape(num2cell(Rxxun, [1,2,3]), 1, U);
    for u = 1:U
        Rxxun_cell{u} = Rxxun_cell{u}(1:Lxu(u), 1:Lxu(u),:);
        Rxxun_cell{u} = reshape(num2cell(Rxxun_cell{u},[1,2]), 1, N);
    end
    Rxxun_cell = [Rxxun_cell{:}];
    Rxx = blkdiag(Rxxun_cell{:});
    maximize (1/cb*log2(exp(1))*log_det(eye(N*Ly) + Hexpand*Rxx*Hexpand'))
    subject to
        trace(Rxx) <= Esum;
cvx_end

% cvx_begin quiet
% cvx_solver mosek
%     variable Rxx1(Lx*N) nonnegative
%     variable Rxx2(Lx*N) nonnegative
%     maximize 0.5*(1/cb*log2(exp(1))*log_det(eye(N*Ly) + Hexpand*diag(Rxx1)*Hexpand')) + 0.5*(1/cb*log2(exp(1))*log_det(eye(N*Ly) + Hexpand*diag(Rxx2)*Hexpand'))
%     subject to
%         0.5*sum(Rxx1) + 0.5*sum(Rxx2) <= Esum;
%         Rxx2(10)==0;
% cvx_end

bsum=cvx_optval;

Rxxun_cell = reshape(num2cell(Rxxun, [1,2,3]), U, 1);
for u = 1:U
    Rxxun_cell{u} = Rxxun_cell{u}(1:Lxu(u), 1:Lxu(u),:);
    Rxxun_cell{u} = reshape(num2cell(Rxxun_cell{u},[1,2]), 1, N);
end
Rxxun_cell = vertcat(Rxxun_cell{:}); % U*N cell array

Rxx = zeros(Lx,Lx,N);
for n=1:N
    Rxx(:,:,n) = blkdiag(Rxxun_cell{:,n});
end

bs=zeros(1,U);
bsum_lin=0;
for u=1:U
    indices=idx_start(u):idx_end(u);
    for n=1:N
     bs(u)=bs(u)+(1/cb)*(log2(det(eye(Ly)+H(:,:,n)*Rxx(:, ...
       :,n)*H(:,:,n)')) - log2(det(eye(Ly)+H(:,:,n)*Rxx(:, ...
       :,n)*H(:,:,n)'- H(:,indices,n)*Rxx(indices, ...
       indices,n)*H(:,indices,n)')));
    end
    bsum_lin=bsum_lin+real(bs(u));
end

end

% % handle some special simple cases separately to speed it up
% if max(Lxu)==min(Lxu)
%     cvx_begin
%         variable Rxxu(Lxu(1),Lxu(1),N,U) hermitian semidefinite
%         tmp = num2cell(Rxxu, [1,2]);
%         maximize (1/cb*log2(exp(1))*log_det(eye(N*Ly) + Hexpand*blkdiag(tmp{:,:})*Hexpand'))
%         subject to
%             trace(sum(sum(Rxxu,3),4)) <= Esum;
%     cvx_end
% elseif U==1
%     cvx_begin
%         variable Rxx(Lx,Lx,N) hermitian semidefinite
%         tmp = num2cell(Rxx, [1 2]);
%         maximize (1/cb*log2(exp(1))*log_det(eye(N*Ly) + Hexpand*blkdiag(tmp{:})*Hexpand'))
%         subject to
%             trace(sum(Rxx, 3)) <= Esum;
%     cvx_end
% elseif U==2
%     cvx_begin
%         variable Rxx1(Lxu(1),Lxu(1),N) hermitian semidefinite
%         variable Rxx2(Lxu(2),Lxu(2),N) hermitian semidefinite
%         tmp1 = num2cell(Rxx1, [1 2]);
%         tmp2 = num2cell(Rxx2, [1 2]);
%         maximize (1/cb*log2(exp(1))*log_det(eye(N*Ly) + Hexpand*blkdiag(tmp1{:},tmp2{:})*Hexpand'))
%         subject to
%             trace(sum(Rxx1, 3)) + trace(sum(Rxx2, 3)) <= Esum;
%     cvx_end
% elseif U==3
%     cvx_begin
%         variable Rxx1(Lxu(1),Lxu(1),N) hermitian semidefinite
%         variable Rxx2(Lxu(2),Lxu(2),N) hermitian semidefinite
%         variable Rxx3(Lxu(3),Lxu(3),N) hermitian semidefinite
%         tmp1 = num2cell(Rxx1, [1 2]);
%         tmp2 = num2cell(Rxx2, [1 2]);
%         tmp3 = num2cell(Rxx3, [1 2]);
%         maximize (1/cb*log2(exp(1))*log_det(eye(N*Ly) + Hexpand*blkdiag(tmp1{:},tmp2{:}, tmp3{:})*Hexpand'))
%         subject to
%             trace(sum(Rxx1, 3)) + trace(sum(Rxx2, 3)) + trace(sum(Rxx3, 3)) <= Esum;
%     cvx_end
