function [Rxxs, E, theta, bun] = minPMAC_var_Lxu(H, Lxu, bu_min, w, cb)
% minPMAC_var_Lxu solves minimize sum_u(w_u*tr(Rxx(u))), subject to
% bu>=bu_min
% This program calls startEllipse_var_Lxu and minPtone_cvx_var_Lxu. CVX
% package with mosek solver is needed.
%
% The inputs are: 
% 1)    H, an Ly by Lx by N channel matrix. Ly is the number of receiver
%       antennas, Lx is the total number of transmit antennas and N is the
%       total number of tones.
%       H(:,:,n) is the channel matrix for all users on tone n.
%       To generate H from a time-domain channel h, use fft(h, N, 3)
% 2)    Lxu, a scalar or a length-U vector containing the number of
%       transmit antennas of each user. If Lxu is a scalar, each user has
%       Lxu antennas.
% 3)    bu_min, a U by 1 vector containing the target rates for each user.
% 4)    w, a U by 1 vector containing the weights for each user's power.
% 5)    cb, a scalar indicator of whether a real-baseband channel (cb=2) or
%       a complex-baseband channel (cb=1) is used. The default value is
%       cb=1.

% The outputs are:
% 1)    Rxxs, cell array containing the Rxx's for each user on each tone.
%       Rxx{n}{u} is the Rxx for user u on tone n.
% 2)    E, a U by N matrix containing the transmit power of each user on
%       each tone.
% 3)    theta, the optimal U by 1 dual variable vector containing optimal
%       weights of rates. theta also determines the decoding order.
% 4)    bun, a U by N matrix containing the rates of each user on each tone

tstart=tic;
if nargin < 5
    cb = 1;
end
bu_min = bu_min/(3-cb); % rate calculations are w.r.t. bit/real-dimension

err = 1e-9;                              % error tolerance                             
[~, ~, N] = size(H);
U = length(w);
if length(Lxu) == 1
    Lxu = ones(1,U)*Lxu;
end
w = reshape(w, U, 1);
bu_min = reshape(bu_min, U, 1);

% initialize
count = 0;
bun = zeros(U,N);
E = zeros(U,N);
Rxxs = cell(1,N);
[A, theta] = startEllipse_var_Lxu(H, Lxu, bu_min, w);

bu_min = bu_min*log(2);     % conversion from bits to nuts 

while 1
    % solve for each tone
    parfor tone = 1:N
        [~, bun(:,tone), Rxxs{tone}, E(:,tone)] = minPtone_cvx_var_Lxu(H(:,:,tone), Lxu, theta, w); 
    end
    % update ellipsoid for theta
    g = sum(bun,2) - bu_min;                   
    % stopping criteria
    tmp_err = sqrt(g'*A*g)
    if tmp_err <= err           
        break
    end
    gt = g/tmp_err;
    tmp = A*gt;
    theta = theta - 1/(U + 1)*tmp;
    A = U^2/(U^2 - 1)*(A - 2/(U + 1)*(tmp*tmp'));
    ind = find(theta < zeros(U,1));
    
    while ~isempty(ind)                  % This part is to make sure that theta is feasible,
        g = zeros(U,1);                  % it was not covered in the lecture notes and you may skip this part
        g(ind(1)) = -1;
        tmp = A * g / sqrt(g' * A * g);
        theta = theta - 1 / (U + 1) * tmp;
        A = U^2 / (U^2 - 1) * (A - 2 / (U + 1) * (tmp * tmp'));
        ind = find(theta < zeros(U,1));
    end   
    count = count+1
end

if max(Lxu) == 1
    Rxxs = E;
end
bun = bun /log(2)*(3-cb);          % conversion from nut-rates/real-dimension to bit rates
toc(tstart)
end
    