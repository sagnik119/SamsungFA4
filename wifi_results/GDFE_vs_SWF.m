function [E_GDFE, b_GDFE, E_lin, b_lin] = GDFE_vs_SWF(H, noise_var, power)
    
    cb = 1;
    [Ly, Lx, fft_length] = size(H);
    w   = ones(1, Lx);
    Lxu = ones(1, Lx);
    Rnn = repmat(eye(Ly), [1, 1, fft_length]);
    
    
    [Rxx, bsum, bsum_lin, bu_a_lin, bun_lin] = linear_only( ...
        fft_length/(fft_length+10) * ones(1, Lx) * power, ...
        H/sqrt(noise_var), ...
        Lxu , Rnn, cb);
    tic 
%     [Eun, bun, ~, ~, ~] = minPMAC_reformulated(H/sqrt(noise_var), bu_a_lin, w, cb, ones(2^Lx, fft_length));  
    [Eun, ~, bun] = maxRMAC_cvx(H/sqrt(noise_var), ... 
                                fft_length/(fft_length+10) * ones(1, Lx) * power * fft_length,...
                                ones(1, Lx), cb);
    toc 

     E_GDFE = sum(Eun, "all");
     b_GDFE = sum(bun, "all");
     E_lin  = sum(Rxx, "all") ;
     b_lin  = sum(bu_a_lin, "all");
end