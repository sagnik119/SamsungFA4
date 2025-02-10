function [Rnn, sumRatebar, S1, S2, S3, S4, S5, Z] = cvx_wcnoise(Rxx, H, Lyu)
%cvx_wcnoise This function computes the worst-case noise for a given input
%autocorrelation Rxx and channel matrix.
%   Arguments:
%       - Rxx: input autocorrelation, size(Lx, Lx)
%       - H: channel response, size (Ly, Lx)
%       - Lyu: number of antennas at each user, scalar/vector of length U
%   Outputs:
%       - Rnn: worst-case noise autocorrelation, with white local noise
%       - sumRatebar: maximum sum rate/real-dimension

[Ly,Lx] = size(H);
if length(Lyu) == 1
    U = Ly/Lyu;
    Lyus = ones(1,U)*Lyu;
else
    U = length(Lyu);
    Lyus = Lyu;
end
cum_Lyu = cumsum(Lyus);
cum_Lyu = [0, cum_Lyu(1:end-1)];
us = 1:U;

[Q,D,~]=svd(Rxx);
rD = rank(D);
sqD = sqrt(diag(D));
sqD(rD+1:end) = 0;
Htilde = H*Q*diag(sqD);
cvx_begin sdp quiet
cvx_solver mosek
cvx_precision high
variable Rnn(Ly, Ly) hermitian 
variable Z(Lx, Lx) hermitian
dual variables S1{U} S2{U} S3 S4 S5
maximize log_det(eye(Lx) - Z)
subject to
for u = us
    S1{u}: eye(Lyus(u))==real(Rnn(cum_Lyu(u) + (1:Lyus(u)), cum_Lyu(u) + (1:Lyus(u))));
    S2{u}: 0==imag(Rnn(cum_Lyu(u) + (1:Lyus(u)), cum_Lyu(u) + (1:Lyus(u))));
end
S3: Rnn == hermitian_semidefinite(Ly);
S4: [Rnn + Htilde*Htilde', Htilde; Htilde', Z] == hermitian_semidefinite(Ly+Lx);
S5: Z == hermitian_semidefinite(Lx);
cvx_end

sumRatebar = -0.5*log2(exp(1))*cvx_optval;
