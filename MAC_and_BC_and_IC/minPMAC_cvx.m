function [Eun, theta, bun] = minPMAC_cvx(H, bu_min, w, cb)
% minPMAC_cvx solves minimize sum_u(w_u*tr(Rxx(u))), subject to
% bu>=bu_min
% This function uses CVX package

% the inputs are: 
% 1)  H, an Ly by U by N channel matrix. Ly is the number of receiver antennas, 
%     U is the total number of users and N is the total number of tones.
%     H(:,:,n) is the channel matrix for all users on tone n
%     and H(:,u,n) is the channel for user u on tone n. In this code we assume each user 
%     only has single transmit antenna, thus H(:,u,n) is a column vector.
%     To get H from a time-domain channel h, use fft(h, N, 3)
% 2)  bu, a U by 1 vector containing the target rates for all the users.
% 3)  w, a U by 1 vector containing the weights for each user's power.
% 4)  cb, a scalar indicator of whether a real-baseband channel (cb=2) or
%     a complex-baseband channel (cb=1) is used. The default value is cb=1.

% the outputs are:
% 1)  Eu, a U by N matrix containing the powers for all users and over all tones
%     that minimizes the weighted-sum power. E(u,:) is the power allocation for
%     user u over all tones.
% 2)  theta, the optimal U by 1 dual variable vector containing optimal weights
%     of rates. theta determines the decoding order.
% 3)  bu, a U by N matrix containing the rates of all users on all tones after convergence.

tstart = tic;
if nargin < 4
    cb = 1;
end
bu_min = bu_min / (3-cb);
err_threshold = 1e-2;

[Ly, U, N] = size(H);
w = reshape(w,U,1);
bu_min = reshape(bu_min,U,1);

% initialize
err = 1;
[A, theta] = startEllipse(H, bu_min, w);         % starting ellipsoid 
bu_min = bu_min*log(2);     % conversion from bits to nuts 

count = 1;
bun = zeros(U,N);
Eun = zeros(U,N);
while (err > err_threshold)
    % Rxx step
    parfor n = 1:N
        [~,bun(:,n), Eun(:,n)] = minPtone_cvx(H(:,:,n), theta, w)
    end
    % theta step
    g = sum(bun,2) - bu_min;                   
    % stopping criteria
    err = sqrt(g'*A*g)
    if err <= err_threshold           
        break
    end
    gt = g/err;
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
    count = count + 1
end
bun = bun /log(2)*(3-cb);          % conversion from nut-rates/real-dimension to bit rates
toc(tstart)
end
    
  