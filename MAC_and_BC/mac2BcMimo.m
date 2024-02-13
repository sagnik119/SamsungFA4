% This function converts the covariace matrices in a Gaussian MAC to corresponding
% covariance matrices in its Dual BC. These covariance matrices achieve the
% same set of rates in the dual BC by dirty paper coding scheme. The
% encoding order in the BC is the reverse of the decoding order in the MAC.
% The total number of users is denoted by K. The decoding order in the MAC
% is given by K by 1 vector pi. pi(k) is the user that is decoded kth in
% the successive decoding. H is a t by r by K matrix containing all the
% channel matrices of the MAC. H(:,:,k) is the channel matrix for user k in
% the MAC. t and r are the number of transmit and receive antennas in the
% BC. S is a r by r by K matrix containing the covariance matrices of the
% MAC. S(:,:,k) is the covariance matrix for user k in the MAC. G is the
% function output which will contain the covariance matrices of the BC.
% G(:,:,k) is the covariance matrix for user k. This function is for just
% an MIMO-BC and does not include parallel MIMO-BCs.

function G = mac2BcMimo(S, H, pi)


H = H(:,:,pi);
S = S(:,:,pi);


[t, r, K] = size(H);

Gtot = zeros(t,t);

% Gtot is the transmit covariance matrix of the BC

B = zeros(t,t,K);
A = zeros(r,r,K);

% A and B matrices are the Rtildenoise for BC and MAC respectively


B(:,:,K) = eye(t);

for k = K:-1:2
   B(:,:,k-1) = B(:,:,k) + H(:,:,k)*S(:,:,k)*H(:,:,k)';
end

A(:,:,1) = eye(r);

for k = 1:K
   temp_A = inv(sqrtm(A(:,:,k)));
   temp_B = inv(sqrtm(B(:,:,k)));
   
   [F L M] = svd(temp_B * H(:,:,k) * temp_A, 'econ');
   G(:,:,k) = temp_B * F * M' * sqrtm(A(:,:,k)) * S(:,:,k) * sqrtm(A(:,:,k)) * M * F' * temp_B;
   
   Gtot = Gtot + G(:,:,k);
   
   if k~=K
      A(:,:,k+1) = eye(r) + H(:,:,k+1)' * Gtot * H(:,:,k+1);
   end
end

G(:,:,pi) = G;
