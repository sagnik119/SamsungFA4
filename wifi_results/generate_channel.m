clc
clear
close all

addpath('quadriga_src')
addpath('wifi_results')

U = 3;
manual_fading = 0;


fc = 5e9;
BW = 80e6;
fft_length = 64;

noise_power_dbm = -96;
noise_var = 10^((noise_power_dbm - 30)/10)/BW;

SNR_db = 10:10:80;
SNR = 10.^(SNR_db/10);

Eu_level = SNR * noise_var;

num_AP_antennas = 2;

HH = zeros(num_AP_antennas, U, fft_length, 1000);

for count = 1:1000
    disp(count)
    
    s = qd_simulation_parameters;
    s.use_absolute_delays = 1;
    s.center_frequency = fc;
    s.show_progress_bars = 0;

    l = qd_layout(s);
    l.no_tx = 1;
    l.no_rx = U;
    

    
    h = zeros(l.no_rx, num_AP_antennas, fft_length);
    
    theta = [0, 5, 10, 15, 20, 25] * pi/180;
    dist  = [3, 3, 3, 3, 3, 3];
        
    
    for antenna = 1:num_AP_antennas
        l.tx_position(:,1) = [0  0  3]';
        
        for i = 1:l.no_rx
            l.rx_position(:, i) = dist(i) * [cos(theta(i)), sin(theta(i)), 1.5];
        end
        
        
            
        l.set_scenario('WINNER_Indoor_A1_LOS'); % Example: NLOS scenario
        p = l.init_builder;
        
        
        if manual_fading == 1
            [m, n] = size(p);
            for i = 1:m
                for j = 1:n
                    disp('Removing Quadriga Fading')
                    p(i,j).plpar = [];                                           % Disable path-loss
                end
            end
        end
        
        p.gen_parameters;                                       % Generate small-scale-fading parameters
        c = p.get_channels;                                     % Generate channel coefficients
        
        return
        
        
        [fad_m, fad_std] = get_fading("B", dist, fc);
        
        
        
        for i = 1:numel(c)
           split = strsplit(c(i).name, '_');
           r_idx = split(3);
           r_idx = r_idx{1};
           r_idx = str2num(r_idx(3:end));
           
           t_idx = split(2);
           t_idx = t_idx{1};
           t_idx = str2num(t_idx(3:end));
        
           fad_db = fad_m(r_idx) + randn()*fad_std;
        
           
           if manual_fading == 1
               disp('Adding Manual Fading')
               h(r_idx, antenna, :) = c(i).fr(BW, fft_length) * 10^(-fad_db/20); % * 10^2.5;% * 10^2.5;% * 1e2;% * 1e2;
           else
               h(r_idx, antenna, :) = c(i).fr(BW, fft_length);% * 10^(-fad_db/20); % * 10^2.5;% * 10^2.5;% * 1e2;% * 1e2;
           end
        end
    end
    H = permute(h, [2 1 3]);
    HH(:,:,:, count) = H;

end


filename = sprintf('wifi_results/H_wifi_%d_users.mat', U);
save(filename, 'HH', 'fc', 'BW', 'fft_length', 'noise_var', 'SNR_db', 'SNR', 'Eu_level')

    
    