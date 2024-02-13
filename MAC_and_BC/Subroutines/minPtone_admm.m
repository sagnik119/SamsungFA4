function [f, bn, Rxxs, En] = minPtone_admm(H, Lxu, theta, w)
% minPtone minimizes f = - sum_{u=1}^U theta_u * b_u + sum_{u=1}^U w_u *
% e_u + sum_{u=1}^U log(e_u)
% subject to b \in C_g(G,e)

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
% 3)  En, a U by 1 vector containing the powers Of each user

[Ly, ~] = size(H);
U = length(theta);
if length(Lxu) == 1
    Lxu = ones(1,U)*Lxu;
end

[stheta, order] = sort(theta, 'descend');
sdeltas = -diff([stheta;0]);
Hu = mat2cell(H, Ly, Lxu);
sHu = Hu(order)';
Q = cell(U,1);
R = cell(U,1);
for u = 1:U
    [Q{u}, R{u}] = qr(pinv(sHu{u}'),0);
end
sw = w(order);

err = 1;
err_threshold = 1e-6;
iter = 1;
max_iter = 1e6;
rho = 1/norm(H)^2/U;

coeff = min((U:-1:1)',repmat(U:-1:1,U,1)) + eye(U);
%inv_coeff = inv(coeff);

% initialize Rxx, T, A, Y, Z
sRxxs = cell(U,1);
Ts = zeros(Ly, Ly, U);
Ys = zeros(Ly, Ly, U);
Zs = zeros(Ly, Ly, U);
dZs = zeros(Ly, Ly, U);
tmp_T = zeros(Ly, Ly, U);
for u = 1:U
    sRxxs{u} = eye(Lxu(order(u)));
    %sRxxs{u} = zeros(Lxu(order(u)));
    Ts(:,:,u) = sHu{u}*sRxxs{u}*sHu{u}';
end
As = cumsum(Ts, 3) + repmat(eye(Ly), 1, 1, U);

while err > err_threshold && iter < max_iter
    % Rxx(u,n) and A_un update
    tmp_Rxx = Ts - Ys;
    tmp_A = As - dZs - Zs;
    for u = 1:U
        [P,D] = eig(Q{u}'*tmp_Rxx(:,:,u)*Q{u} - sw(u)/rho*R{u}*R{u}');
        D = max(real(D),0);
        sRxxs{u} = R{u}'*P*D*P'*R{u};
        [V, Lambdas] = eig(tmp_A(:,:,u));
        Lambdas = real(diag(Lambdas));
        As(:,:,u) = V*diag((Lambdas+sqrt(Lambdas.^2+4*sdeltas(u)/rho))/2)*V';
    end
    % T_u update
    for u = 1:U
        tmp_T(:,:,u) = sHu{u}*sRxxs{u}*sHu{u}'; 
    end
    tmp_rhs = cumsum(As - repmat(eye(Ly), 1, 1, U) + Zs, 3, 'reverse') + Ys + tmp_T;
    rhs_flat = reshape(permute(tmp_rhs,[2,1,3]),Ly,Ly*U)';
    %Ts = kron(inv_coeff, eye(Ly))*rhs_flat;
    Ts = kron(coeff, eye(Ly))\rhs_flat;
    Ts = permute(reshape(Ts',Ly, Ly, U),[2,1,3]);
    % Y_u update
    dYs = tmp_T - Ts;
    Ys = Ys + dYs;
    % Z_u update
    dZs = As - cumsum(Ts, 3) - repmat(eye(Ly), 1, 1, U);
    Zs = Zs + dZs;
    %dZs{:}
    % calculate error
    err = max(max(abs(dYs),[],'all'), max(abs(dZs),[],'all'));
    iter = iter + 1;
end

for u = 1:U
    Ts(:,:,u) = sHu{u}*sRxxs{u}*sHu{u}';
end
As = cumsum(Ts, 3) + repmat(eye(Ly), 1, 1, U);

Rxxs(order) = sRxxs;
Rxxs = Rxxs';
En = zeros(U,1);
cum_bn = zeros(U,1);
for u = 1:U
    En(u) = real(trace(Rxxs{u}));
    cum_bn(u) = real(log(det(As(:,:,u))));
end
bn(order) = 0.5*diff([0;cum_bn]);
bn = bn';
f = -theta'*bn + w'*En; 

end

