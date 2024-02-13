% function Rxxb = mac2bc(Rxxm, Hvecu) 
%
% mac2bc converts Gaussian MAC autocorrelation matrices to the 
% corresponding autocorrelation matrices of its dual Gaussian BC.
% These autocorrelation matrices achieve 
% the same user rates for a GDFE MAC and a GDFE noiseless precoder. The
% BC encoding order reverses the MAC decoding order in the input order pi.
%
% U is the number of users; Lx (per-user) and Ly are the number 
% of MAC transmit and receive antennas respectively. These are extracted/
% inferred from the dimensions of inputs Hvecu and Rxxm.  The MAC user order is
% reversed.
%
% Inputs
%
% Rxxm is an Lx by Lx by U matrix containing the MAC input Rxx matrices. 
%   where Rxxm(:,:,u) is the covariance matrix for MAC user u. 
% Hvecu is an Ly by Lx by U matrix that contains all MAC user-channel matrices.
%   Hvecu(:,:,u) is user u's channel matrix, the transpose of Hvecutrans 
%   in bc2mac prog.
%
% Output
%
% Rxxb is an Ly x Ly x U array of dual-BC autocorrelation matrices.
%   Rxxb(:,:,u) is user u's aucorrelation matrix.
%
% When there are unequal numbers of dimensions per MAC-user input, so 
% Lx --> Lxu, the input autocorrelation matrices Rxxm should all use the 
% maximum size Lxmax and thus place zeros in any rows/columns necessary on 
% some users to get to Lxmax.  In this situation, the correspond H_u should
% be extended with zeroed columns in the corresponding zeroed-input Rxxm(u)
% dimensions. 

function Rxxb = mac2bc(Rxxm, Hvecu)

H = Hvecu(:,:,end:-1:1);
Rxxm = Rxxm(:,:,end:-1:1);


[Ly, ~, U] = size(H);
Lx = size(Rxxm, 1);

Rxxbsum = zeros(Ly,Ly);
Rxxb = zeros(Ly,Ly,U);

% Rxxbsum is the transmit covariance matrix of the BC

B = zeros(Ly,Ly,U);
A = zeros(Lx,Lx,U);

% A and B matrices are the noise-plus-xtalk for BC and MAC respectively.


B(:,:,U) = eye(Ly);

for u = U:-1:2
   B(:,:,u-1) = B(:,:,u) + H(:,:,u)*Rxxm(:,:,u)*H(:,:,u)';
end

A(:,:,1) = eye(Lx);

for u = 1:U
   temp_A = inv(sqrtm(A(:,:,u)));
   temp_B = inv(sqrtm(B(:,:,u)));
   
   [F,~,~] = svd(temp_B * H(:,:,u) * temp_A, 'econ');
   Rxxb(:,:,u) = temp_B * F * M' * sqrtm(A(:,:,u)) * Rxxm(:,:,u) * sqrtm(A(:,:,u)) * M * F' * temp_B;
   
   Rxxbsum = Rxxbsum + Rxxb(:,:,u);
   
   if u~=U
      A(:,:,u+1) = eye(Lx) + H(:,:,u+1)' * Rxxbsum * H(:,:,u+1);
   end
end

% Rxxb(:,:,pi) = Rxxb; 
