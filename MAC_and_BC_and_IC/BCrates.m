function [bbc] = BCrates(H, Rxx, order)
%BCrates Calculate BC individual user rates. Assume each (sub)user has
%one receive antenna
%   Inputs:
%       H: Lx by U channel matrix
%       Rxx: Lx by Lx by U matrix. Rxx(:,:,u) denotes the autocorrelation
%       matrix of user u.
%       order: length U vector. The Tx encodes from order(end) to order(1)
%   Output:
%       bbc: bit rate of each user

[~,U] = size(H);
bbc = zeros(1,U);
H = H(order,:);
Rxx = Rxx(:,:,order);
Rxxcum = cat(3, zeros(U), cumsum(Rxx, 3));
for u = 1:U
    bbc(u) = log2(1 + H(u,:)*Rxxcum(:,:,u+1)*H(u,:)')...
        - log2(1 + H(u,:)*Rxxcum(:,:,u)*H(u,:)');
bbc(order) = real(bbc);
end

