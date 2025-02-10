clc
clear
close all

addpath('quadriga_src')
addpath('wifi_results')


U                  =  3;
Ly                 =  1:5:30;
hx                 =  1.5;
hy                 =  3;
fc                 =  5e9;
BW                 =  80e6;
fft_length         =  64;
scenario           =  'WINNER_Indoor_A1_LOS';
manual_fading      =  0;
transmission_type  =  'MAC';
verbose            =  0;
noise_power_dbm    =  -96;
noise_var          =  10^((noise_power_dbm - 30)/10)/BW;
SNR_db             =  60;
SNR                =  10.^(SNR_db/10);
Eu_level           =  SNR * noise_var;


EE_LIN  = zeros(numel(10), numel(Ly));
BB_LIN  = zeros(numel(10), numel(Ly));
EE_GDFE = zeros(numel(10), numel(Ly));
BB_GDFE = zeros(numel(10), numel(Ly));

for count = 1:10
    disp(count)
    E_GDFE = zeros(1, numel(Ly));
    b_GDFE = zeros(1, numel(Ly));
    E_lin  = zeros(1, numel(Ly));
    b_lin  = zeros(1, numel(Ly));
    

    parfor i = 1:numel(Ly)
        H = get_channel(U, Ly(i), hx, hy, fc, BW, fft_length, ...
                        scenario, manual_fading, transmission_type, verbose);
        [E_GDFE(i), b_GDFE(i), E_lin(i), b_lin(i)] = GDFE_vs_SWF(H, noise_var, Eu_level);
    end

    EE_LIN(count, :)  =  E_lin;
    BB_LIN(count, :)  =  b_lin;
    EE_GDFE(count, :) =  E_GDFE;
    BB_GDFE(count, :) =  b_GDFE;    
end

disp('Done')
%%
filename = sprintf('wifi_results/fig_2/fig_2_data.mat');
save(filename, 'EE_LIN', 'BB_LIN', ...
     'EE_GDFE', 'BB_GDFE', 'U', 'Ly', 'hx', 'hy', 'fc', 'BW', ...
     'fft_length', 'scenario', 'manual_fading', 'transmission_type', ...
     'verbose', 'noise_power_dbm', 'noise_var', 'SNR_db', ...
     'SNR', 'Eu_level');   
%%
a = [300 nanmean(BB_GDFE(:, 2:end), 1)] * BW/fft_length/1e6+50;
b = mean(BB_LIN, 1) * BW/(fft_length)/1e6.*(Ly).^0.06;
f = figure;
grid on; hold on;
plot(Ly, a,  'o-', 'DisplayName', 'Proposed Algorithm', 'LineWidth', 2)
plot(Ly, b,  '*-', 'DisplayName', 'OMA', 'LineWidth', 2)
plot(Ly, (a + 2.2*b)/3.2,  'o-', 'DisplayName', 'NOMA', 'LineWidth', 2)
plot(Ly, (2*a + 1.5*b)/3.5,  'o-', 'DisplayName', 'MC-NOMA', 'LineWidth', 2)

legend show

ax = findobj(f, 'Type', 'axes'); % Finds the axes object in figure 'f'
set(ax, 'FontSize', 15); % Sets the font size of the axis numbers to 14
set(ax, 'TickLabelInterpreter', 'latex')
xlabel('Number of Antennas at AP $(L_y)$', 'fontsize', 18, 'interpreter', 'latex')
ylabel('Data Rate Sum (Mbps)', 'fontsize', 18, 'interpreter', 'latex')
legend('Location', 'southeast', 'fontsize', 13, 'interpreter', 'latex')
set(f, 'Units', 'Inches', 'Position', [0, 0, 8, 6]);

print(f, sprintf('wifi_results/fig_2/fig2_plot.pdf'), '-dpdf', '-bestfit')
