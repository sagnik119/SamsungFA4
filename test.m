% To use get_channel and get_fading to generate channels
% Use generated channel to run minPMAC program

addpath('channel_files', 'FA4', 'MAC_and_BC', 'single_user', 'tmp', 'utils', ...
    'FA4\Wi-Fi_ChannelModels_B_D', 'MAC_and_BC\Subroutines');

%% Generating channels
[h1, t_dly, t_dly_10ns] = get_channel('B', 1, 2, 1);
[h2, t_dly, t_dly_10ns] = get_channel('B', 1, 2, 2);
[h3, t_dly, t_dly_10ns] = get_channel('B', 1, 2, 3);

[fad_m1, fad_std1] = get_fading('B', 3, 2.4e9);
fad_db1 = fad_m1 + randn()*fad_std1;
[fad_m2, fad_std2] = get_fading('B', 5, 2.4e9);
fad_db2 = fad_m2 + randn()*fad_std2;
[fad_m3, fad_std3] = get_fading('B', 8, 2.4e9);
fad_db3 = fad_m3 + randn()*fad_std3;

h1 = h1*10^(-fad_db1/20);
h2 = h2*10^(-fad_db2/20);
h3 = h3*10^(-fad_db3/20);
h = cat(2, h1, h2, h3);
H = fft(h, 16, 3);

%% Applying minPMAC on generated channel
[Eun, theta, bun, FEAS_FLAG, bu_a, info] = minPMAC(H*30000, 16*[2,3,4], [1,1,1], 1);
%[Eun, theta, bun] = minPMAC_cvx(H*30000, 32*[2.5,2.5,2.5], [1,1,1], 1)