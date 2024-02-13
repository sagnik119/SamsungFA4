function [bmac] = MACrates(H,Rxx,order)
%MACrates Calculate MAC individual user rates. Assume each (sub)user has
%one transmit antenna
%   Inputs:
%       H: U by Ly channel matrix
%       Rxx: U by U diagonal matrix
%       order: length U vector. The Rx decodes from order(1) to order(end).
%   Output:
%       bmac: bit rate of each user

[~,Lx] = size(H);
H = H(order(end:-1:1), order(end:-1:1))*sqrtm(Rxx);
Rbinv = H'*H+eye(Lx);
G = chol(Rbinv);
bmac = log2(diag(G).^2);
bmac = bmac(order(end:-1:1));
end

