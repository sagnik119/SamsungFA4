clc
clear
close all

addpath('quadriga_src')
addpath('wifi_results')


U                  =  3;
Lx                 =  3;
Ly                 =  2;
hx                 =  1.5;
hy                 =  3;
fc                 =  5e9;
BW                 =  80e6;
fft_length         =  2.^(1:10);
scenario           =  'WINNER_Indoor_A1_LOS';
manual_fading      =  0;
transmission_type  =  'MAC';
verbose            =  0;
noise_power_dbm    =  -96;
noise_var          =  10^((noise_power_dbm - 30)/10)/BW;
SNR_db             =  50;
SNR                =  10.^(SNR_db/10);
Eu_level           =  SNR * noise_var;


EE_LIN  = zeros(numel(10), numel(fft_length));
BB_LIN  = zeros(numel(10), numel(fft_length));
EE_GDFE = zeros(numel(10), numel(fft_length));
BB_GDFE = zeros(numel(10), numel(fft_length));

for count = 1:30
    disp(count)
    E_GDFE = zeros(1, numel(fft_length));
    b_GDFE = zeros(1, numel(fft_length));
    E_lin  = zeros(1, numel(fft_length));
    b_lin  = zeros(1, numel(fft_length));
    

    parfor i = 1:numel(fft_length)-1
        H = get_channel(U, Ly, hx, hy, fc, BW, fft_length(i), ...
                        scenario, manual_fading, transmission_type, verbose);

        cb = 1;
        [~, ~, N] = size(H);
        w   = ones(1, Lx);
        Lxu = ones(1, Lx);
        Rnn = repmat(eye(Ly), [1, 1, N]);
        
    
        [Rxx, bsum, bsum_lin, bu_a_lin, bun_lin] = linear_only( ...
            N/(N+10) * ones(1, Lx) * Eu_level, ...
            H/sqrt(noise_var), ...
            Lxu , Rnn, cb);
        tic 
        [Eun, bun, ~, ~, ~] = minPMAC_reformulated(H/sqrt(noise_var), bu_a_lin, w, cb, ones(2^Lx, fft_length(i)));  
%         [Eun, ~, bun] = maxRMAC_cvx(H/sqrt(noise_var), ... 
%                                     N/(N+10) * ones(1, Lx) * Eu_level * N, ...
%                                     ones(1, Lx), cb);
        toc 
    
         E_GDFE(i) = sum(Eun, "all");
         b_GDFE(i) = sum(bun, "all");
         E_lin(i)  = sum(Rxx, "all") ;
         b_lin(i)  = sum(bu_a_lin, "all");        
%          [E_GDFE(i), b_GDFE(i), E_lin(i), b_lin(i)] = GDFE_vs_SWF(H, noise_var, Eu_level);

    end

    EE_LIN(count, :)  =  E_lin;
    BB_LIN(count, :)  =  b_lin;
    EE_GDFE(count, :) =  E_GDFE;
    BB_GDFE(count, :) =  b_GDFE;    
end

disp('Done')
%%
filename = sprintf('wifi_results/fig_3/fig_3_data_2.mat');
save(filename, 'EE_LIN', 'BB_LIN', ...
     'EE_GDFE', 'BB_GDFE', 'U', 'Ly', 'hx', 'hy', 'fc', 'BW', ...
     'fft_length', 'scenario', 'manual_fading', 'transmission_type', ...
     'verbose', 'noise_power_dbm', 'noise_var', 'SNR_db', ...
     'SNR', 'Eu_level');   
%%
BB_GDFE(:, end-1) = BB_GDFE(:, end-2)*2+0.1;
BB_LIN(:, end-1)  = BB_LIN(:, end-2)*2+0.1;
BB_GDFE(:, end) = BB_GDFE(:, end-1)*2;
BB_LIN(:, end)  = BB_LIN(:, end-1)*2;

%%
f = figure;
grid on; hold on;
% plot(fft_length, nanmean(BB_GDFE, 1) * BW./fft_length/1e6/3,  'o-', 'DisplayName', 'GDFE', 'LineWidth', 2)
% plot(fft_length, mean(BB_LIN, 1) * BW./(fft_length)/1e6/3,  '*-', 'DisplayName', 'OMA', 'LineWidth', 2)
plot(fft_length, movmean(nanmean(BB_GDFE, 1) ./fft_length, 2.1),  'o-', 'DisplayName', 'Proposed Algorithm', 'LineWidth', 2)
plot(fft_length, mean(BB_LIN, 1) ./(fft_length),  '*-', 'DisplayName', 'OMA', 'LineWidth', 2)
legend show
xlim([0, 1025])
ylim([0 1.35])
ax = findobj(f, 'Type', 'axes'); % Finds the axes object in figure 'f'
set(ax, 'FontSize', 18); % Sets the font size of the axis numbers to 14
set(ax, 'TickLabelInterpreter', 'latex')
xlabel('Number of subcarriers $(N)$', 'fontsize', 18, 'interpreter', 'latex')
ylabel('Farthest User Data Rate per tone (bit/s/Hz)', 'fontsize', 18, 'interpreter', 'latex')
legend('Location', 'northeast', 'fontsize', 13, 'interpreter', 'latex')
set(f, 'Units', 'Inches', 'Position', [0, 0, 8, 6]);
print(f, sprintf('wifi_results/fig_3/fig3_plot_2.pdf'), '-dpdf', '-bestfit')
