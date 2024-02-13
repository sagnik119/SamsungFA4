% function [Rxx, bsum , bsum_lin] = SWF(Eu, H, Lxu, Rnn, cb)
%
% Simultaneous water-filling MAC max rate sum (linear and nonlinear GDFE)
% The input is space-time domain h, and the user can specify a temporal
% block symbol size N (essentially an FFT size). 
%
% Inputs: 
% Eu  U x 1 energy/SAMPLE vector. Single scalar equal energy all users
%     any (N/N+nu) scaling should occur BEFORE input to this program.    
% H   The FREQUENCY-DOMAIN Ly x sum(Lx(u)) x N MIMO channel for all users.
%     N is determined from size(H) where N = # used tones 
% Lxu 1xU vector of each user's number of antennas
% Rnn The Ly x Ly x N noise-autocorrelation tensor (last index is per tone)
% cb  cb = 1 for complex, cb=2 for real baseband
%       cb=2 corresponds to a frequency range at an sampling rate 1/T' of
%       [0, 1/2T'] while with cb=1, it is [0, 1/T'].  The Rnn entered for
%       these two situations may differ, depending on how H is computed.
%
% Outputs:
% Rxx A block-diagonal psd matrix with the input autocorrelation for each
%     user on each tone. Rxx has size (sum(Lx(u)) x sum(Lx(u)) x N .
%     sum trace(Rxx) over tones and spatial dimensions equal the Eu 
% bsum the maximum rate sum.
% bsum  bsum_lin - the maximum sum rate with a linear receiver
%     b is an internal convergence sum rate value, not output
%
% This program significantly modifies one originally supplied by student
% Chris Baca

function [Rxx, bsum, bsum_lin, bu_a_lin, bun_lin] = linear_only(Eu, H, Lxu , Rnn, cb)
U = numel(Lxu);
[Ly, Lx_sum, N] = size(H);

if numel(Eu) == 1
    Eu = Eu*ones(U);
end

i = 0;
for n=1:N
Rxx(:,:,n) = diag(zeros(Lx_sum, N));
H(:,:,n)=inv(sqrtm(Rnn(:,:,n)))*H(:,:,n);
end
b = zeros(U, 1);
bsum=sum(b);

cum_index=1;
user_ind(1)=1;
for u=2:U
    cum_index = cum_index + Lxu(u);
    user_ind(u) = cum_index;
end

for u = 1:U
   
    st = user_ind(u);
    
    if u < U
        en = user_ind(u+1)-1;
    else
        en = Lx_sum;
    end
    
    for n = 1:N
        Rxx(st:en, st:en,n) = Eu(u)*eye(en-st+1);
    end
end

%b = zeros(U);
while 1
    bsum_prev = bsum;
    b_prev = b;
    
    for u = 1:U
        %new user
        st = user_ind(u);

            if u < U
                en = user_ind(u+1)-1;
            else
                en = Lx_sum;
            end
        M_u = zeros(en-st+1,en-st+1, N);
        g_u = [];
       Rnn_un=zeros(Ly,Ly,N);
        for n = 1:N
            %new tone
            Rnn_un(:,:,n) = eye(Ly);
            for v = 1:U
                if v ~= u
                    st = user_ind(v);
                    if v < U
                        en = user_ind(v+1)-1;
                    else
                        en = Lx_sum;
                    end

                    Hvn = H(:, st:en, n);
                    Rxxvn = Rxx(st:en, st:en, n);
                    Rnn_un(:,:,n) = Rnn_un(:,:,n) + Hvn*Rxxvn*Hvn';
                end
            end
            
            st = user_ind(u);

            if u < U
                en = user_ind(u+1)-1;
            else
                en = Lx_sum;
            end
            
            
            Hun = H(:, st:en, n);         
            Hun_til = sqrtm(inv(Rnn_un(:,:,n)))*Hun;

            [F, D, M_n] = svd(Hun_til);
            s = svd(Hun_til);
            
            M_u(:,:,n) = M_n;
            g_u(:,n) = s.^2;
           
        end
        
        g_flat = reshape(g_u, [1, numel(g_u)]);
%        [g_sort, ind] = sort(g_flat, 'descend');

        [B(u,:),E(u,:),L_star] = waterfill_gn(g_flat,Eu(u),0, cb);
        
        L = numel(g_flat);
%        Etot = N*Eu(u);
%        j = L;
        
%        e = zeros(1,L);
%       size(E)
%        e(1:L_star)=E;
        
        e = reshape(E(u,:), [en-st+1, N]);
        b=sum(B(u,:));
        
        for n = 1:N
            Rxx(st:en, st:en, n) = M_u(:,:,n)*diag(e(:,n))*M_u(:,:,n)';
        end
    end
        bsum=0;
        for n=1:N
        bsum=bsum+(1/cb)*log2(det(eye(Ly)+H(:,:,n)*Rxx(:,:,n)*H(:,:,n)'));
        end
        bsum=real(bsum);

    i = i+1;
    
    if abs(bsum-bsum_prev) <= 1e-6
  % if norm(b-b_prev) <= 0
        break
    end
    if i>1000
        i
        break
    end
end
%b = b(:,1);
%Hcell = mat2cell(H,Ly,Lx_sum,ones(1,N));
%Hexpand = [blkdiag(Hcell{1,1,:})]
%Rcell = mat2cell(Rxx,Ly,Lx_sum,ones(1,N));
%Rcell = mat2cell(Rxx,Lx_sum,Lx_sum,ones(1,N));
%Rxxexpand = [blkdiag(Rcell{1,1,:})]
%bsum = log2(det(eye(Lx_sum*N)+Hexpand*Rxxexpand*Hexpand'));
%bsum = log2(det(eye(Ly*N)+Hexpand*Rxxexpand*Hexpand'));
%bsum=0;
%for n=1:N
%    bsum=bsum+(1/cb)*log2(det(eye(Ly)+H(:,:,n)*Rxx(:,:,n)*H(:,:,n)'));
%end
%bsum=real(bsum);
bs=zeros(1,U);
bsum_lin=0;
bu_lin = zeros(U, N);
for u=1:U
    indices=user_ind(u):user_ind(u)+Lxu(u)-1;
    for n=1:N
     bu_lin(u, n) = (1/cb)*(log2(det(eye(Ly)+H(:,:,n)*Rxx(:, ...
       :,n)*H(:,:,n)')) - log2(det(eye(Ly)+H(:,:,n)*Rxx(:, ...
       :,n)*H(:,:,n)'- H(:,indices,n)*Rxx(indices, ...
       indices,n)*H(:,indices,n)')));

     bs(u)=bs(u)+(1/cb)*(log2(det(eye(Ly)+H(:,:,n)*Rxx(:, ...
       :,n)*H(:,:,n)')) - log2(det(eye(Ly)+H(:,:,n)*Rxx(:, ...
       :,n)*H(:,:,n)'- H(:,indices,n)*Rxx(indices, ...
       indices,n)*H(:,indices,n)')));
    end
    bsum_lin=bsum_lin+real(bs(u));
end
bu_a_lin = real(bs);
bun_lin = real(bu_lin);