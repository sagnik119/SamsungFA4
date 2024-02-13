% This function initializes the ellipsoid method for power minimization
% problem for MAC, with variable number of transmit antennas at each user.
% Called by minPMAC_new program

% function [A, g] = startEllipse_var_Lxu(H, Lxu, bu_min, w)
% H containes the channel matrices as is defined in minPMAC.m, Lxu is
% either a scalar or a length-U vector indicating number of transmit
% antennas of each user. bu_min is the target U by 1 rate vector, w is the
% U by 1 power weights. 
% For all fixed margin WF steps, decoding order is 1,2,...U.

% A is the matrix describing the staring ellipsoid and g is its center

function [A, g] = startEllipse_var_Lxu(H, Lxu, bu_min, w)


[Ly, ~, N] = size(H);
U = length(w);
if length(Lxu) == 1
    Lxu = ones(1,U)*Lxu;
end
index_end = cumsum(Lxu);
index_start = [1,index_end(1:end-1)+1];

order = 1:U;                       % decoding order, arbitrary
tmax = zeros(U,1);

for u = 1:U
    Rxxs = cell(U,N);
    for u1 = 1:U
        for tone = 1:N
            Rxxs{u1,tone} = zeros(Lxu(u1));
        end
    end
    b = bu_min;
    b(u) = b(u) + 1;
    en = zeros(U, N);
    for u1 = order
        M = cell(1,N);
        gs = zeros(Lxu(u1),N);
        for tone = 1:N
            Rnoise = eye(Ly) +...
                H(:,[1:index_start(u1)-1,index_end(u1)+1:end],tone)...
                *blkdiag(Rxxs{[1:u1-1,u1+1:U], tone})...
                *H(:,[1:index_start(u1)-1,index_end(u1)+1:end],tone)';
            [~,h_tmp,M{tone}] = svd(Rnoise^(-1/2)*H(:,index_start(u1):index_end(u1),tone),0);
            gs(1:length(diag(h_tmp)),tone) = diag(h_tmp).^2;
        end
        [~, e_tmp] = fmwaterfill_gn(reshape(gs,1,[]), b(u1)/N, 0);
        en_tmp = reshape(e_tmp, Lxu(u1), N);
        en(u1,:) = sum(en_tmp, 1);
        for tone = 1:N
            Rxxs{u1,tone} = M{tone}*diag(en_tmp(:,tone))*M{tone}';
        end
    end
    
    tmax(u) = w'* sum(en,2);
end

g = tmax / 2;
A = g'*g*eye(U);        % an sphere centerced at tmax/2






    
    