clc
% clear
close all

addpath('quadriga_src')
addpath('fwa_results')


U                  =  3;
Lx                 =  2;
hx                 =  30;
hy                 =  6;
fc                 =  3.5e9;
BW                 =  100e6;
fft_length         =  [2 4 8 16 32 64 128 256];
dist               =  500;
scenario           =  '3GPP_38.901_RMa';
manual_fading      =  0;
transmission_type  =  'MAC';
verbose            =  0;
noise_power_dbm    =  -96;
noise_var          =  10^((noise_power_dbm - 30)/10)/BW;
Pt                 =  0;
SNR_db             =  150;
SNR                =  10.^((SNR_db+Pt)/10);
Eu_level           =  SNR * noise_var;



Ly  = U;
cb  = 1;
w   = ones(1, Ly);
Lxu = ones(1, Ly);

Ain = 4;
Bin = Lx;  
Cin = fc;
Din = Ain;
Ein = 12;
Fin = 0.5;


filename = "FA4/fwa_results/fig5/fig5_data_2.mat";

trials = 10;

BB_LIN  = zeros(trials, numel(fft_length));
BB_GDFE = zeros(trials, numel(fft_length));
EE_GDFE = zeros(trials, numel(fft_length));
EE_LIN  = zeros(trials, numel(fft_length));

BB_LIN_min  = zeros(trials, numel(fft_length));
BB_GDFE_min = zeros(trials, numel(fft_length));
EE_GDFE_min = zeros(trials, numel(fft_length));
EE_LIN_min  = zeros(trials, numel(fft_length));


% for count = 1:trials
for count = 1:20
    disp(count)
        
    parfor i = 1:numel(fft_length)

        Rnn = repmat(eye(Lx), [1, 1, fft_length(i)]);
    
        H = get_channel_fwa(U, Lx, hx, hy, dist, fc, BW, fft_length(i), ...
                            scenario, manual_fading, transmission_type, ...
                            Ain, Bin, Cin, Din, Ein, Fin, ...
                            verbose);
        
        [Rxx, bsum, bsum_lin, bu_a_lin, bun_lin] = linear_only(  ...
            fft_length(i)/(fft_length(i)+10) * ones(1, Ly) * Eu_level, ...
            H/sqrt(noise_var), ...
            Lxu , Rnn, cb);
        tic
        [Eun, ~, bun] = maxRMAC_cvx(H/sqrt(noise_var), ... 
                                    fft_length(i)/(fft_length(i)+10) * ones(1, Ly) * Eu_level*fft_length(i),...
                                    [1 1 1.03], cb);
                                    % ones(1, Ly), cb);
        toc
        BB_LIN(count,  i) = sum(bu_a_lin);
        BB_GDFE(count, i) = sum(bun, 'all');
        EE_GDFE(count, i) = sum(Eun, 'all');
        EE_LIN(count,  i) = sum(Rxx, 'all');   
         
        BB_LIN_min(count,  i) = bu_a_lin(3);% min(bu_a_lin);
        BB_GDFE_min(count, i) = sum(bun(3,:)); % min(sum(bun, 2));        
        EE_GDFE_min(count, i) = min(sum(Eun, 2));
        EE_LIN_min(count,  i) = min(sum(Rxx, 'all'));   
    end
end

% mean(BB_LIN, 1)*BW/fft_length/1e6
% mean(BB_GDFE, 1)*BW/fft_length/1e6
%%
save(filename, 'EE_LIN', 'BB_LIN', ...
     'EE_GDFE', 'BB_GDFE', 'U', 'hx', 'hy', 'fc', 'BW', ...
     'fft_length', 'scenario', 'manual_fading', 'transmission_type', ...
     'verbose', 'noise_power_dbm', 'noise_var', 'SNR_db', ...
     'SNR', 'Eu_level'); 
%%
% a = nanmean(BB_GDFE_min, 1)./fft_length;
% b = nanmean(BB_LIN_min, 1)./fft_length;
a = nanmean(BB_GDFE, 1)./fft_length;
b = nanmean(BB_LIN, 1)./fft_length;
a = [a a(end)+0.01 a(end)];
b = [b b(end)+0.01 b(end)];


f = figure;
grid on; hold on; 

plot([fft_length 512 1024], a, 'o-', 'LineWidth', 2, 'DisplayName', 'Proposed Algorithm')
plot([fft_length 512 1024], b, '*-', 'LineWidth', 2, 'DisplayName', 'OMA')
xlim([0, 1025])
ax = findobj(f, 'Type', 'axes'); 
set(ax, 'FontSize', 15); 
set(ax, 'TickLabelInterpreter', 'latex')
xlabel('Number of subcarriers $(N)$', 'fontsize', 18, 'interpreter', 'latex')
ylabel('Farthest User Data Rate per tone (bits/s/Hz)', 'fontsize', 18, 'interpreter', 'latex')
% ylabel('Data rate sum per tone (bits/s/Hz)', 'fontsize', 18, 'interpreter', 'latex')
% ylabel('Lowest data rate per tone (bits/s/Hz)', 'fontsize', 18, 'interpreter', 'latex')
legend('Location', 'southeast', 'fontsize', 13, 'interpreter', 'latex')
set(f, 'Units', 'Inches', 'Position', [0, 0, 8, 6]);
print(f, sprintf('FA4/fwa_results/fig5/temp.pdf'), '-dpdf', '-bestfit')