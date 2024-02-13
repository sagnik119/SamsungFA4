function [Rxx, b_linear, b_gdfe, Rxxu] = linear_sol(H, Lxu, Eu, cb)
%linear_sol Computes the achievable rate of a linear system given channel
%and power constraint
%   Input arguments:
%       - H: Ly-by-Lx-by-N channel matrix
%       - Lxu: each user's number of transmit antennas. Scalar or length-U
%           vector. If scalar, all users have the same number of transmit
%           antennas
%       - Eu: each user's power contraint. Length-U vector
%       - cb: 1 if complex channel, 2 if real baseband channel
%       - fixband: 1 if total bandwidth is fixed, 0 if tonal bandwidth is
%           fixed
%   Outputs:
%       - Rxx: Lx-by-Lx-by-N block diagonal matrix, representing the
%       combined transmit autocorrelation matrix
%       - b_linear: sum rate/dimension for linear receiver, scalar
%       - b_gdfe: max sum rate/dimension for GDFE receiver, scalar
%       - Rxxu: Each user's Rxx. U-by-N cell array, Rxxu{u,n} is a
%           Lxu(u)-by-Lxu(u) matrix representing the transmit
%           autocorrelation matrix of user u on tone n

U = length(Eu);
Eu = reshape(Eu, [], 1);
[Ly, Lx, N] = size(H);
if length(Lxu) == 1
    Lxu = Lxu*ones(1, U);
elseif length(Lxu) ~= U
    error('Invalid length of Lxu');
else
    Lxu = reshape(Lxu, 1, U);
end
if sum(Lxu) ~= Lx
    error('Invalid value of Lxu');
end
idx_end = cumsum(Lxu);
idx_start = [1, idx_end(1:end-1) + 1];


% to speed the program up, separate the cases when U is small
switch U
    case 1
        cvx_begin quiet
        cvx_solver mosek
            variable Rxx(Lxu,Lxu,N) hermitian semidefinite
            expressions bn(N) E
            for n = 1:N
                bn(n) = log_det(eye(Ly) + N*H(:,:,n)*Rxx(:,:,n)*H(:,:,n)');
                E = E + trace(Rxx(:,:,n));
            end
            maximize sum(bn)
            subject to
                E <= Eu;
        cvx_end
        Rxxu = reshape(mat2cell(Rxx, Lxu, Lxu, ones(1,N)), 1, N);
    case 2 
        cvx_begin quiet
        cvx_solver mosek
            variable Rxx1(Lxu(1),Lxu(1),N) hermitian semidefinite
            variable Rxx2(Lxu(2),Lxu(2),N) hermitian semidefinite
            expressions bn(N) E(2)
            for n = 1:N
                bn(n) = log_det(eye(Ly) + N*H(:,:,n)*blkdiag(Rxx1(:,:,n), Rxx2(:,:,n))*H(:,:,n)');
                E(1) = E(1) + trace(Rxx1(:,:,n));
                E(2) = E(2) + trace(Rxx2(:,:,n));
            end
            maximize sum(bn)
            subject to
                E <= Eu;
        cvx_end
        Rxx = zeros(Lx, Lx, N);
        Rxxu = cell(2,N);
        for n = 1:N
            Rxx(:,:,n) = blkdiag(Rxx1(:,:,n), Rxx2(:,:,n));
            Rxxu(:,n) = [{Rxx1(:,:,n)}; {Rxx2(:,:,n)}];
        end
    case 3
        cvx_begin quiet
        cvx_solver mosek
            variable Rxx1(Lxu(1),Lxu(1),N) hermitian semidefinite
            variable Rxx2(Lxu(2),Lxu(2),N) hermitian semidefinite
            variable Rxx3(Lxu(3),Lxu(3),N) hermitian semidefinite
            expressions bn(N) E(3)
            for n = 1:N
                bn(n) = log_det(eye(Ly) + N*H(:,:,n)*blkdiag(Rxx1(:,:,n), Rxx2(:,:,n), Rxx3(:,:,n))*H(:,:,n)');
                E(1) = E(1) + trace(Rxx1(:,:,n));
                E(2) = E(2) + trace(Rxx2(:,:,n));
                E(3) = E(3) + trace(Rxx3(:,:,n));
            end
            maximize sum(bn)
            subject to
                E <= Eu;
        cvx_end
        Rxx = zeros(Lx, Lx, N);
        Rxxu = cell(3,N);
        for n = 1:N
            Rxx(:,:,n) = blkdiag(Rxx1(:,:,n), Rxx2(:,:,n), Rxx3(:,:,n));
            Rxxu(:,n) = [{Rxx1(:,:,n)}; {Rxx2(:,:,n)}; Rxx3(:,:,n)];
        end
    otherwise
        Lxu_max = max(Lxu);
        
        cvx_begin quiet
        cvx_solver mosek
            variable Rxxs(Lxu_max, Lxu_max, U, N) hermitian
            expressions tmp(Ly, Ly, N) bn(N) E(U)
            for n = 1:N
                for u = 1:U
                    tmp(:,:,n) = tmp(:,:,n) + N*H(:,idx_start(u):idx_end(u),n)...
                        *Rxxs(1:Lxu(u),1:Lxu(u),u,n)...
                        *H(:,idx_start(u):idx_end(u),n)';
                    E(u) = E(u) + trace(Rxxs(1:Lxu(u),1:Lxu(u),u,n));
                end
                bn(n) = log_det(eye(Ly) + tmp(:,:,n));
            end
            maximize sum(bn)
            subject to
                for u = 1:U
                    for n = 1:N
                        Rxxs(1:Lxu(u),1:Lxu(u),u,n) == semidefinite(Lxu(u));
                    end
                end
                E <= Eu;
        cvx_end
        
        Rxx = zeros(Lx, Lx, N);
        Rxxu = cell(U,N);
        
        for u = 1:U
            Rxxu(u,:) = reshape(mat2cell(Rxxs(1:Lxu(u),1:Lxu(u),u,:), Lxu(u), Lxu(u), 1, ones(1,N)), 1, N);
        end
        for n = 1:N
            Rxx(:,:,n) = blkdiag(Rxxu{:,n});
        end
end

bun_lin = zeros(U,N);
bn_gdfe = zeros(1,N);
for n = 1:N
    tmp_total = eye(Ly) + N*H(:,:,n)*Rxx(:,:,n)*H(:,:,n)';
    bn_gdfe(n) = 1/cb*log2(real(det(tmp_total)));
    for u = 1:U
        tmp_u = tmp_total - N*H(:,idx_start(u):idx_end(u),n)*Rxxu{u,n}*H(:,idx_start(u):idx_end(u),n)';
        bun_lin(u,n) = bn_gdfe(n) - 1/cb*log2(real(det(tmp_u)));
    end
end

b_linear = sum(bun_lin, 'all');
b_gdfe = sum(bn_gdfe);
end

