% This function calculates the rate pair (b_1,b_2) such that b_2 = alpha b_1 and 
% it lies on the boundary point of the capacity region of the MAC in problem 13.11.
function b = capregion(alpha)

G = zeros(2,2,2);
G(:,:,1) = [-0.5 0.1; -1.7 0.2];    % Whitened equivalent channel
G(:,:,2) = [-1.2 1.1;1.1 -0.1];
Eu = [1;1];                         % Energy constraints

b_l = 0;
b_u = 10;

err = 1e-6;                         % error tolerance

while b_u - b_l > err
    
    b = (b_l + b_u) / 2;
    [E, r, theta, w, f] = admMAC(G, [b; alpha * b], Eu);
    if f >= 0
                                    % (b, alpha b) is feasible with given
                                    % energy constraints
        b_l = b;                    % increase b_l
    else
        b_u = b;                    % otherwise decrease b_u
    end
end
