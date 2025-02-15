% Conducts comparison of linear receiver based algorithm with minPMAC
%
% get_channel variables:
% model_letter          : Channel model letter, 'A' to 'F'
% N_tx                  : Number of transmit antennas
% N_rx                  : Number of receive antennas
% index                 : index of the channel, 1 to 5000
% h                     : The channel matrix, N_rx by N_tx by N_taps
% tap_delay             : the tap delays in s
% tap_delay_10ns        : the tap delays in 10ns
%
% get_fading variables:
% dist                  : distance between Tx and Rx in meter
% fc                    : carrier frequency
% fad_mean_db           : mean of the fading (in dB)
% fad_std_db            : standard deviation of the log-normal fading (in dB)
%
% minPMAC variables:
% H                     : channel matrix of dimensions (Ly * U * N)
% noise_var             : noise variance
% bu_min                : (U * 1) matrix containing minimum required data rates for each user
% w                     : (U * 1) matrix containing the weights for each user's power
% cb                    : 1 if complex baseband, or 2 if real baseband
% Eun                   : U by N energy distribution that minimizes the weighted-sum energy. 
%                           E(u,n) is user u's energy allocation on tone n.
% theta                 :the optimal U by 1 dual variable vector containing optimal weights
%                           of rates. Theta determines the decoding order. Largest theta is
%                           decoded last, and smallest first.
%
% bun                   :U by N bit distributions for all users.
%
% FEAS_FLAG             : indicator of achievability. 
%                       FEAS_FLAG=1 if the target is achieved by a single ordering; 
%                       FEAS_FLAG=2 if the target is achieved by time-sharing
%
% bu_a                  : U-by-1 vector showing achieved sum rate of each user. 
%
% info                  : various length output depending on FEAS_FLAG
%                           --if FEAS_FLAG=1: 1 x 4 cell array containing
%                           {Rxxs, Eun, bun, theta} corresponds to the single vertex
%                           there are no equal-theta user sets in this case
%                           --if FEAS_FLAG=2: 1-x 6 cell array, with each row representing
%                           a time-shared vertex {Rxxs, Eun, bun, theta, frac}
%                           there are numclus equal-theta user sets in this case.
%
%                           info's row entries in detail (one row for each vertex shared
%                                - Eun:  U-by-N matrix showing users' transmit energy on each tone.
%                                        If infeasible, output 0.
%                                - bun: U-by-N matrix showing users' rate on each tone. If
%                                        infeasible, output 0.
%                                - theta: U-by-1 Lagrangian multiplier w.r.t. target rates
%                                - order: produces the order from left(best) to right for vertex
%                                - frac: fraction of dimensions for each vertex in time share (FF =2 ONLY)
%                                - cluster: index (to which cluster the user belongs; 0 means no cluster)

% Setup
clc 

addpath('channel_files', 'FA4', 'MAC_and_BC', 'single_user', 'tmp', 'utils', ...
    'FA4/Wi-Fi_ChannelModels_B_D', 'MAC_and_BC/Subroutines', 'FA4/ML', ...
    'FA4/ML/outputs');

try
    load FA4/ML/outputs/Data_SIMO.mat
catch
    clear all;
    Data_SIMO = [];
    save FA4/ML/outputs/Data_SIMO.mat
end

% Number of data samples to generate
num_samples = 1;

% fix random generator seed for reproducible results
% rng(43);

%% Define ranges or sets for input parameters
channel_types = 'B';        
users_range = 3; % No. of users
N_tx = 1; % no. of antennas per user
N_rx_range = 2; % No. of antennas at AP
index_range = 1:5000; % channel sample index
dist_range = [1,5]; % have to generate random float value between 1m and 5m
fc = 2.4e9; % carrier frequency
fft_length_range = 2.^[6]; % range of FFT lengths
noise_var = 1.0/9800000;
transmit_power_dB = linspace(-10, 51, 50);
transmit_power = 10.^(transmit_power_dB/10);

bu_min_single_user_single_tone_range = 1:5; % to be multiplied by FFT length to get total energy per user
theta_single_user_range = 1;
cb = 1; % complex baseband channel

% fixing samples to monitor average data statistics
index = [1, 2, 3] +100;
% get_fading parameters
dist = [3, 3, 3];


%% Preallocate storage for data samples
data_samples = struct();

% Initialize the wait bar
wait_bar = waitbar(0, 'Initializing...');

idx = 1;
while idx <= num_samples
    % try
        % Update the wait bar with the current progress
        waitbar(idx / num_samples, wait_bar, sprintf('Processing: %d of %d', idx, num_samples));
    
        fprintf ("Starting data aggregation for index %d\n", idx);
        % Randomly select values from the specified ranges or sets
        num_users = users_range;
        % get_channel parameters
        channel_type = randsample(channel_types, 1);
        N_rx = N_rx_range;
        % fft parameters
        fft_length = fft_length_range(1);
        theta = theta_single_user_range*ones(1, num_users);
        % minPMAC parameters
    
    
    
    
        % bu_min = [4.479 4.479 4.479];%5 * ones(1, num_users); % bit/Hz
    
    
    
    
    
        % Print parameters used for this iteration
        fprintf ("Parameters: num_users: %d , N_rx: %d, fft_length: %d\n", ...
            num_users, N_rx, fft_length );
    
        fprintf ("Generating channel\n");
        % Generating channels with random values
        try
            [h_user, t_dly, t_dly_10ns] = get_channel(channel_type, N_tx, ...
                    N_rx, index(1)); % to get number of delay taps
        catch exception
            fprintf ("Channel file missing: skipping to next iteration\n")
            continue; % skip this iteration since channel file missing
        end
        h = zeros(N_rx, num_users*N_tx, size(t_dly, 2)); 
        for user=1:num_users
            [h_user, t_dly, t_dly_10ns] = get_channel(channel_type, N_tx, ...
                N_rx, index(user));
            
            % [temp, ~,~] = size(h_user);
            % for tempp = 1:temp
            %     h_user(tempp,1,:) = channel_resample(h_user(tempp,1,:), 100e6, 4, 5, 9, 0);
            % end
            % h_user = channel_resample(h_user, 100e6, 4, 5, 8, 0);
            [fad_m, fad_std] = get_fading(channel_type, dist(user), fc);
            fad_db = fad_m + randn()*fad_std;
            h_user = h_user*10^(-fad_db/20);
            h(:, (user-1)*N_tx+1:user*N_tx, :) = h_user;
        end
    
        % FFT
        H = fft(h, fft_length, 3);
    
        % Gathering linear receiver data rates
        Lxu = ones(1, num_users);
        Rnn = repmat(eye(N_rx), [1, 1, fft_length]);
        Eun = eye(num_users);
        [Rxx, bsum, bsum_lin, bu_a_lin, bun_lin] = linear_only(fft_length/(fft_length + 8)* ...
            sum(Eun, 2)/fft_length, H/sqrt(noise_var), ...
            Lxu , Rnn, cb);


    
        % fprintf ("Starting minPMAC\n");
        % Applying minPMAC on generated channel with random values
        Eu = ones (1, num_users);
        [Eun, w, bun] = maxRMAC_cvx(H/sqrt(noise_var), Eu, theta);
    
    
        % fprintf ("minPMAC finished\n");
    
        % fprintf ("Starting maxRMAC\n");
        % Applying maxRMAC on generated channel with random values
        % [Eun, w, bun] = maxRMAC_cvx(H/noise_var,...
        % sum(Eun, 1), ones(1, num_users));
        % fprintf ("maxRMAC finished\n");
    
        % Store inputs and outputs
        data_samples(idx).N_tx = N_tx;
        data_samples(idx).fc = fc;
        data_samples(idx).num_taps = size(t_dly, 2);
        data_samples(idx).num_users = num_users;
        data_samples(idx).channel_type = channel_type;
        data_samples(idx).N_rx = N_rx;
        data_samples(idx).sample_index = index;
        data_samples(idx).dist = dist;
        data_samples(idx).fft_length = fft_length;
        % data_samples(idx).bu_min = bu_min;
        data_samples(idx).w = w;
        data_samples(idx).h = h;
        data_samples(idx).H = H;
        data_samples(idx).Eun = Eun ;
        % data_samples(idx).theta = theta;
        data_samples(idx).bun = bun  ;
        % data_samples(idx).FEAS_FLAG = FEAS_FLAG;
        % data_samples(idx).bu_a = bu_a ;
        % data_samples(idx).info = info;
        data_samples(idx).bu_a_lin = bu_a_lin ;
        data_samples(idx).b_sum_lin = bsum_lin ;
        data_samples(idx).b_sum = bsum ;
        data_samples(idx).bun_lin = bun_lin ;
        data_samples(idx).Rxx = Rxx;
    
        % save collected data cumulatively after every sample
        temp_filename = strcat('FA4/comparison_models/outputs/data_samples_compare_linear_to_minPMAC' ...
            , int2str(idx), '.mat');
        save(temp_filename, 'data_samples');
    
        Data_SIMO = [Data_SIMO, data_samples(idx)];
        save('FA4/comparison_models/outputs/Data_SIMO.mat', 'Data_SIMO')

        try
            load FA4/comparison_models/outputs/Data_SIMO.mat
        catch
            clear all;
            Data = [];
            save FA4/comparison_models/outputs/Data_SIMO.mat
        end
    
        idx = idx+1;
    % catch
    %     disp('TTTTTTTRRRRRRRRRYYYYYYYY   AAAAAAGGGGGGAAAAAAIIIIINNNN')
    %     continue
    % end
end

% Close the wait bar when the loop is done
close(wait_bar);

% Save the data samples to a .mat file
save('FA4/comparison_models/outputs/data_samples.mat', 'data_samples');