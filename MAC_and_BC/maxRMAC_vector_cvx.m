function [Rxxs, Eun, w, bun] = maxRMAC_vector_cvx(H, Lxu, Eu, theta)
%maxRMAC_vector_cvx Maximize weighted rate sum subject to power constraint
%   Input arguments:
%       - H: Ly-by-Lx-by-N channel matrix. H(:,:,n) denotes the channel for
%           the n-th tone.
%       - Lxu: number of transmit antennas of each user. It can be either a
%           scalar or a length-U vector. If it is a scalar, every user has
%           Lxu transmit antennas; otherwise user u has Lxu(u) transmit
%           antennas.
%       - Eu: Power constraint on each user, length-U vector.
%       - theta: weight on each user's rate, length-U vector.
%   Outputs:
%       - Rxxs: U-by-N cell array containing Rxx(u,n)'s if Lxu is a
%           length-U vector; or Lxu-by-Lxu-by-U-by-N tensor if Lxu is a
%           scalar.
%       - Eun:  U-by-N matrix showing users' transmit power on each tone
%       - w: U-by-1 Lagrangian multiplier w.r.t. power constraints
%       - bun: U-by-N matrix showing users' rate on each tone

UNIFORM_FLAG = 0;
[Ly, ~, N] = size(H);
if N == 1
    H = reshape(H,Ly,[],1);
end
Eu = reshape(Eu,[],1);
theta = reshape(theta,[],1);
[stheta, idx] = sort(theta, 'descend');
delta = -diff([stheta;0]);
U = length(Eu);
if length(Lxu) == 1
    Lxu = ones(1,U)*Lxu;
    UNIFORM_FLAG = 1;
end
Lxu_max = max(Lxu);
index_end = cumsum(Lxu);
index_start = [1,index_end(1:end-1)+1];
Lxu = Lxu(idx);
index_start = index_start(idx);
index_end = index_end(idx);
cvx_begin quiet
cvx_solver mosek
    variable rxxs(Lxu_max, Lxu_max, U, N) hermitian
    dual variable w
    expressions r(U,N) S(Ly, Ly) Eun(U,N)
    for n = 1:N
        S=eye(Ly);
        for u = 1:U
            S = S + H(:,index_start(u):index_end(u),n)...
                *rxxs(1:Lxu(u),1:Lxu(u),u,n)...
                *H(:,index_start(u):index_end(u),n)';
            r(u,n) = log_det(S);
            Eun(u,n) = trace(rxxs(1:Lxu(u),1:Lxu(u),u,n));
        end
    end
    maximize 0.5*sum(delta'*r)
    subject to
        w: sum(Eun,2)<=Eu;
        for u = 1:U
            for n=1:N
                rxxs(1:Lxu(u),1:Lxu(u),u,n)==hermitian_semidefinite(Lxu(u));
            end
        end
cvx_end

cumrate = zeros(U,N);
for n = 1:N
    S=eye(Ly);
    for u = 1:U
        S = S + H(:,index_start(u):index_end(u),n)...
            *rxxs(1:Lxu(u),1:Lxu(u),u,n)...
            *H(:,index_start(u):index_end(u),n)';
        cumrate(u,n) = 0.5*real(log2(det(S)));
    end
end
bun(idx,:) = diff([zeros(1,N); cumrate]);
w(idx)=w;

if UNIFORM_FLAG
    Rxxs(:,:,idx,1:N) = rxxs;
    for u = 1:U
        for n = 1:N
            Eun(idx(u),n) = trace(Rxxs(:,:,u,n));
        end
    end
else
    Rxxs = cell(U,N);
    for u = 1:U
        for n = 1:N
            Rxxs{idx(u),n} = rxxs(1:Lxu(u),1:Lxu(u),u,n);
            Eun(idx(u),n) = trace(Rxxs{idx(u),n});
        end
    end
end

end


