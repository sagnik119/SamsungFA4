function [E, b, theta, w, f] = admMAC(G, bu, Eu)
% This is the ellipsoid method part and the main function to call.

% the inputs are: 
% 1)  G, an Ly by U by N channel matrix. Ly is the number of receiver antennas, 
%     U is the total number of users and N is the total number of tones.
%     G(:,:,n) is the channel matrix for all users on tone n
%     and G(:,u,n) is the channel for user u on tone n. In this code we assume each user 
%     only has single transmit antenna, thus G(:,u,n) is a column vector. 
% 2)  Eu, a U by 1 vector containing power constraints for all the users.
% 3)  bu, a U by 1 vector containing the target rates for all the users.

% the outputs are:
% 1)  E, a U by N matrix containing the powers for all users and over all tones
%     that support the rate vector given in bu. E(u,:) is the power allocation for
%     user u over all tones. An all zero E indicates that the given rates
%     in bu are not achievable with given powers in Eu.
% 2)  b, a U by N matrix containing the rates of all users on all tones after convergence.
%     Again, all zero rates indicates that given bu is not achievable with
%     Eu.
% 3)  theta, the optimal U by 1 dual variable vector containing optimal weights
%     of rates. theta determines the decoding order.
% 4)  w, the optimal U by 1 dual variable vector containing optimal weights
%     of powers.
% 5) f, the dual optimal value.

bu = bu  * log(2);             % conversion from bits to nuts. 
err = 1e-9;                    % error tolerance
w_0 = 1;                       % arbitrary starting ellipsoid.
count = 0;
[Ly, U, N] = size(G);

w = w_0 * ones(U,1);
theta = w_0 * ones(U,1);
A = eye(2 * U) * (2 * U * w_0^2);
g = w_0 * ones(2 * U,1);

while 1  
                               % Calculate the dual function and its
                               % sub-gradient
    [f, b, E] = Dual_adm(G, theta, w, bu, Eu);
 
    if f < -err                % f < 0, terminate, the rates are not achievable
        E = zeros(U,N);
        b = zeros(U,N);
        return
    end                        % Ellipsoid method starts here,
                               % sub-gradient:
    g = [Eu - sum(E,2);sum(b,2) - bu];
    
    if sqrt(g' * A * g) <= err || min(g) >= 0 % stopping critera
        break
    end
                               % Update the ellipsoid
    gt = g / sqrt(g' * A * g);
    temp = [w;theta] - 1 / (2 * U + 1) * A * gt;
    w = temp(1:U);
    theta = temp(U + 1:2 * U);
    A = (2 * U)^2 / ((2 * U)^2 - 1) * (A - 2 / (2 * U + 1) * A * gt * gt' * A);
    
    ind = find([w;theta] < zeros(2 * U,1));
    
    while ~isempty(ind)        % This part is to make sure that [w;theta] is feasible
        g = zeros(2 * U,1);    % it was not covered in lecture notes, you can ignore this part.
        g(ind(1)) = -1;
        gt = g / sqrt(g' * A * g);
        temp = [w;theta] - 1 / (2 * U + 1) * A * gt;
        w = temp(1:U);
        theta = temp(U + 1:2 * U);
        A = (2 * U)^2 / ((2 * U)^2 - 1) * (A - 2 / (2 * U + 1) * A * gt * gt' * A);
        ind = find([w;theta] < zeros(2 * U,1));
    end                        
    count = count + 1;          % Just a counter
end
count
b = b /log(2);                 % conversion from nuts to bits
   
    
    
        

    