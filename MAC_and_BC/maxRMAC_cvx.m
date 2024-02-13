function [Eun, w, bun] = maxRMAC_cvx(H, Eu, theta)
%maxRMAC_cvx Maximize weighted rate sum subject to power constraint, each
%user has ONLY ONE transmit antenna
%   Input arguments:
%       - H: Ly-by-U-by-N channel matrix. H(:,:,n) denotes the channel for
%           the n-th tone.
%       - Eu: Power constraint on each user, length-U vector.
%       - theta: weight on each user's rate, length-U vector.
%   Outputs:
%       - Eun:  U-by-N matrix showing users' transmit power on each tone
%       - w: U-by-1 Lagrangian multiplier w.r.t. power constraints
%       - bun: U-by-N matrix showing users' rate on each tone

[Ly, U, N] = size(H);
Eu = reshape(Eu,[],1);
theta = reshape(theta,[],1);
[stheta, idx] = sort(theta, 'descend');
delta = -diff([stheta;0]);
sH = H(:,idx, :);
Gs = zeros(Ly, Ly, U, N);
parfor n=1:N
    for u=1:U
        Gs(:,:,u,n) = sH(:,u,n)*sH(:,u,n)';
    end
end
cvx_begin quiet
cvx_solver mosek
    variable Eun(U,N) nonnegative
    dual variable w
    expression r(U,N)
    S = cumsum(Gs.*repelem(reshape(Eun,1,1,U,N),Ly,Ly,1,1), 3);
    for u = 1:U
        for n = 1:N
            r(u,n) = log_det(S(:,:,u,n) + eye(Ly));
        end
    end
    maximize 0.5*sum(delta'*r)
    subject to
        w: sum(Eun,2)<=Eu;
cvx_end

%S = cumsum(pagemtimes(Gs,reshape(Eun,1,1,U,N)),3) + eye(Ly);
S = cumsum(Gs.*repelem(reshape(Eun,1,1,U,N),Ly,Ly,1,1), 3) + eye(Ly);
cumrate = zeros(U,N);
for u = 1:U
    for n = 1:N
        cumrate(u,n) = 0.5*real(log2(det(S(:,:,u,n))));
    end
end
bun(idx,:) = diff([zeros(1,N); cumrate]);
w(idx)=w;
end

