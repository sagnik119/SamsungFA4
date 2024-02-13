close all; clear all

fd_dump = '_dump';

n_pkt = 10;
num_tx = 2;
num_rx = 2;
N_fft = 2048;
num_symb = 28;

% noise variance
sigma2 = 1e-4;

% QPSK
mod_order = 4;
qamtable = [ ...
    0.7071    0.7071; ...
    0.7071   -0.7071; ...
    -0.7071    0.7071; ...
    -0.7071   -0.7071];

num_re = N_fft*num_symb;
rnd_idx = randi([1,mod_order],num_re,1);
qam_symbols = qamtable(rnd_idx,1) + 1i*qamtable(rnd_idx,2);
x_tx = qam_symbols(:,ones(1,2));

for nn=1:n_pkt
    % channels in freq., are loaded
    fn_dump = sprintf('%s/ch_fd_dump_pkt%03d.mat',fd_dump,nn);
    dump_mat = load(fn_dump,'ce_fd_dump');

    % channel dimension
    % size(dump_mat.ce_fd_dump) 
    % N_fft x num_symbol x num_rx x num_tx

    ce_fd_dump = dump_mat.ce_fd_dump;
    chan = reshape(ce_fd_dump, [], num_rx, num_tx);

    % y = Hx
    for ii = 1:num_rx
        y_rx(:,ii) = sum(squeeze(chan(:,ii,:)).*x_tx,2);
    end

    % y = Hx + noise
    noise = rand(N_fft*num_symb,1) + 1i*rand(N_fft*num_symb,1);
    y_rx = y_rx + sqrt(sigma2)*noise;

    % demod
    ce_rx = squeeze(sum(chan,3));
    qam_symb_rx = sum(y_rx.*conj(ce_rx),2)./sum(abs(ce_rx).^2,2);

    % EVM computation
    evm_rx = 10*log10(mean(abs(qam_symb_rx(:)-qam_symbols(:)).^2));
    evm_rx_avg = mean(evm_rx);

    if (1)
        figure(1);
        sc_plot(qam_symb_rx(:),'x','r');
        hold on;
        sc_plot(qam_symbols(:),'s','b');
        hold off;
        grid on;

        xlim([-1.3 1.3]);
        ylim([-1.3 1.3]);
    end
end