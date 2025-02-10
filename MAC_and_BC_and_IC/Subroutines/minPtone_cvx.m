function [f, bn, En] = minPtone_cvx(H, theta, w)
% minPtone_cvx minimizes f = - sum_{u=1}^U theta_u * b_{u,n} +
% sum_{u=1}^U w_u * E_{u,n} subject to b \in C_g(H(n),E). This program uses
% cvx package and mosek solver

% the inputs are:
% 1)  H, an Ly by Lx channel matrix. Ly is the number of receiver antennas, 
%     Lx=U is the total number of transmit antennas, where each receiver
%     only has one transmit antenna
% 2)  theta, a U by 1 vector containing the weights for the rates.
% 3)  w, a U by 1 vector containing the weights for each user's power.

% the outputs are:
% 1)  f, -1 * minimum value (or maximum value of the -1 * function).
% 2)  bn, a U by 1 vector containing the rates for each user
%     that optimizes the given function.
% 3)  En, a U by 1 vector containing the power of each user

[Ly, ~] = size(H);
U = length(theta);
[stheta, order] = sort(theta, 'descend');
D = eye(U)+diag(-ones(U-1,1),1);
max_bound = max(-diff([stheta;0]));
if max_bound > 1e8
    stheta = stheta/max_bound*1e4;
    w = w/max_bound*1e4;
end

cvx_begin quiet
cvx_solver mosek
    variable En(U) nonnegative
    expression cumRu(U)
    for v = 1:U
        cumRu(v) = log_det(H(:,order(1:v))*diag(En(order(1:v)))*H(:,order(1:v))'+eye(Ly));
    end
    cumRu = 0.5*cumRu;
    minimize (w'*En - cumRu'*pos(D*stheta))
cvx_end
for v = 1:U
    cumRu(v) = 0.5*log_det(H(:,order(1:v))*diag(En(order(1:v)))*H(:,order(1:v))'+eye(Ly));
end
bn(order) = diff([0;cumRu]);
En(order) = En;
f = -cvx_optval;

end

