function [Rxxbaru] = MAC2BCmimo_new(Rxxu, Hu, pi)
%MAC2BCmimo This function converts the covariace matrices in a Gaussian
%MAC to corresponding covariance matrices in its Dual BC. These covariance
%marices achieve the same set of rates in the dual BC by dirty paper coding
%scheme. The encoding order in the BC is the reverse of the decoding order
%in the MAC.
%   Inputs:
%       Rxxu: Tmac=Rbc by Tmac by U matrix containing the covariance
%       matrices of the MAC.
%       Hu: Rmac=Tbc by Tmac by U matrix containing all the channel
%       matrices of the MAC.
%       pi: U by 1 vector denoting the decoding order in the MAC.     
%   Outputs:
%       Rxxbaru: Tbc by Tbc by U matrix, the covariance matrix of each user
%       in the dual BC

Hu = Hu(:,:,pi);
Rxxu = Rxxu(:,:,pi);
[Tbc, Rbc, U] = size(Hu);
Rmac = Tbc;

% transmit covariance matrices in BC
Rxxbaru = zeros(Tbc, Tbc, U);
% sum of transmit covariance matrices in BC
Rxxbar = zeros(Tbc,Tbc); % Gtot
% noise+interference in MAC
Rnn = zeros(Rmac, Rmac, U); % B
% noise+interference in BC
Rnnbar = eye(Rbc); % A

Rnn(:,:,U) = eye(Rmac);
for u = U:-1:2
    Rnn(:,:,u-1) = Rnn(:,:,u) + Hu(:,:,u)*Rxxu(:,:,u)*Hu(:,:,u)';
end

for u = 1:U
    sRnn = sqrtm(Rnn(:,:,u));
    sRnnbar = sqrtm(Rnnbar);
    Htmp = sRnn\Hu(:,:,u)/sRnnbar;
    [F,~,M] = svd(Htmp, 'econ');
    tmp = sRnn\F*M';
    Rxxbaru(:,:,u) = tmp*(sRnnbar*Rxxu(:,:,u)*sRnnbar)*tmp';
    Rxxbar = Rxxbar + Rxxbaru(:,:,u);
    if u ~= U
        Rnnbar = eye(Rbc) + Hu(:,:,u+1)'*Rxxbar*Hu(:,:,u+1);
    end
end

Rxxbaru(end:-1:1,end:-1:1,pi) = Rxxbaru;

end

