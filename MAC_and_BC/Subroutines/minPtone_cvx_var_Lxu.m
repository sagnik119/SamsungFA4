function [f, bn, Rxxs, En] = minPtone_cvx_var_Lxu(H, Lxu, theta, w)
% minPtone_cvx_var_Lxu minimizes f = - sum_{u=1}^U theta_u * b_u +
% sum_{u=1}^U w_u * e_u 
% subject to b \in C_g(H,e). This program uses cvx package and mosek solver

% the inputs are:
% 1)  H, an Ly by Lx channel matrix. Ly is the number of receiver antennas, 
%     Lx is the total number of transmit antennas.
%     H(:,index_start(u):index_end(u)) is the channel for user u.
% 2)  Lxu, a scalar or a length-U vector containing the number of transmit
%     antennas of each user. If Lxu is a scalar, each user has Lxu antennas
% 3)  theta, a U by 1 vector containing the weights for the rates.
% 4)  w, a U by 1 vector containing the weights for each user's power.

% the outputs are:
% 1)  f, -1 * minimum value (or maximum value of the -1 * function).
% 2)  bn, a U by 1 vector containing the rates for each user
%     that optimizes the given function.
% 3)  Rxxs, a 1 by U cell array containing the Rxx's for each user that
%     optimizes the given function
% 4)  En, a U by 1 vector containing the power of each user

[Ly, ~] = size(H);
U = length(theta);
if length(Lxu) == 1
    Lxu = ones(1,U)*Lxu;
end
Lxu_max = max(Lxu);
index_end = cumsum(Lxu);
index_start = [1,index_end(1:end-1)+1];

[stheta, order] = sort(theta, 'descend');
D = eye(U)+diag(-ones(U-1,1),1);
max_bound = max(-diff([stheta;0]));
if max_bound > 1e8
    stheta = stheta/max_bound*1e4;
    w = w/max_bound*1e4;
end
cvx_begin sdp quiet
    cvx_solver mosek
    variable Rxx(Lxu_max,Lxu_max, U) hermitian semidefinite
    expressions Ru(U) En(U) tmp(Ly, Ly)
    tmp = eye(Ly);
    for v = 1:U
        Hv = H(:,index_start(order(v)):index_end(order(v)));
        rxxv = Rxx(1:Lxu(order(v)),1:Lxu(order(v)),order(v));
        tmp = tmp + Hv*rxxv*Hv';
        Ru(v) = 0.5*log_det(tmp);
        En(order(v)) = trace(rxxv);
    end
    minimize (w'*En - Ru'*pos(D*stheta))
cvx_end

Rxxs = cell(U,1);
tmp = eye(Ly);
for v = 1:U
    Hv = H(:,index_start(order(v)):index_end(order(v)));
    rxxv = Rxx(1:Lxu(order(v)),1:Lxu(order(v)),order(v));
    tmp = tmp + Hv*rxxv*Hv';
    Ru(v) = 0.5*real(log(det(tmp)));
    Rxxs{order(v)} = rxxv;
    En(order(v)) = trace(rxxv);
end
bn = zeros(U,1);
bn(order) = diff([0;Ru]);
f = -cvx_optval;

end

