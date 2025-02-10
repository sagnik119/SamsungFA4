% This function converts the covariace matrices in a Gaussian MAC to corresponding
% covariance matrices in its Dual BC. These covariance matrices achieve the
% same set of rates in the dual BC by dirty paper coding scheme. The
% encoding order in the BC is the reverse of the decoding order in the MAC.
% The total number of users is denoted by K. The decoding order in the MAC
% is given by K by 1 vector ind. ind(k) is the user that is decoded kth in
% the successive decoding. H is a n by K by N matrix containing all the
% channel matrices of the MAC. H(:,k,i) is the channel matrix for user k on tone i in
% the MAC. n is the number transmit antennas in the BC (receiver antennas
% in the MAC). S is a n by n by K by N matrix containing the covariance matrices of the
% each user of BC on each tone. Gtot is a n by n by K matrix containing the
% total covariance matrices of the BC on each tone. r is a K by N matrix 
% containing the individual user rates. This function is 
% written for just one receive antenna per user for the BC and is for 
% parallel BCs and MACs (multi tone).

% the inputs are:
% 1)  H, an n by K by N channel matrix. H(:,k, n) is the channel 
%     vector for user k on tone n. This code assumes each user has just one
%     transmit antenna. n is the number of receive antennas and n=K in DSL.
% 2)  P, a K by N vector of powers for each user on each tone.
% 3)  ind, a K by 1 vector that specified the decoding order in the MAC. 
%     ind(k) is the user that is decoded kth in the successive decoding.


function [S, Gtot, pBC,r] = mac2BcMultiTone(H, P, ind)

% The code is written for decoding order [1 2 3 ... K]

[n, K, N] = size(H);


H = H(:,ind,:);
P = P(ind,:);

Gtot = zeros(n,n,N);
S = zeros(n,n,K,N);

% after re-ordering the matrices, the decoding order is 1,2,3...,K

% Gtot is the transmit covariance matrix of the BC on each tone
pBC = zeros(K,N);
% MAC rates to be matched
r = zeros(K,N);
for i = 1:N
    for k = 1:K
        r(k,i) = log2(det(eye(n) + H(:,k:K,i) * diag(P(k:K, i)) * H(:,k:K,i)') ...
            / det(eye(n) + H(:,k+1:K,i) * diag(P(k+1:K, i)) * H(:,k+1:K,i)'));
    end
    g = 2.^r(:,i) - 1;
    w = zeros(n,K);
    L = zeros(K,K);
    for k = 1:K
        w(:,k) = inv(eye(n) + H(:,k+1:K,i) * diag(P(k+1:K, i)) * H(:,k+1:K,i)') * H(:,k,i);
        w(:,k) = w(:,k) / norm(w(:,k));
        L(k,k) = abs(w(:,k)' * H(:,k,i))^2/g(k);
        for j = 1:k-1
            L(k,j) = -abs(w(:,j)' * H(:,k,i))^2;
        end
    end
    pBC(:,i) = L \ ones(K,1);
    for k = 1:K
        S(:,:,k,i) = w(:,k) * w(:,k)' * pBC(k,i);
        Gtot(:,:,i) = Gtot(:,:,i) + S(:,:,k,i);
    end
end

pBC(ind, :) = real(pBC);
S(:,:,ind,:) = S;
r(ind,:) = real(r);