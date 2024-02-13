% Conducts data aggregation for ML-based solutions to margin-adaptive 
%   multi-user data-transmission
%   For every instance, collected variables are:
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
addpath('channel_files', 'FA4', 'MAC_and_BC', 'single_user', 'tmp', 'utils', ...
    'FA4/Wi-Fi_ChannelModels_B_D', 'MAC_and_BC/Subroutines', 'FA4/ML', ...
    'FA4/ML/outputs');

% Number of data samples to generate
num_samples = 1;

%% Define ranges or sets for input parameters
channel_types = ['B']; 
users_range = 3; % No. of users
N_tx = 2; % no. of antennas per user
% N_rx_range = 4; % No. of antennas at AP
N_rx = 2; % No. of antennas at AP
index_range = 1:5000; % channel sample index
dist_range = [1,5]; % have to generate random float value between 1m and 5m
fc = 2.4e9; % carrier frequency
fft_length_range = 2.^[4]; % range of FFT lengths
noise_var = 1.0/30000; % noise variance
bu_min_single_user_single_tone_range = 1:5; % to be multiplied by FFT length to get total energy per user
w_single_user_range = 1; % have to generate random float value between 1 and 10
cb = 1; % complex baseband channel


%% Preallocate storage for data samples
data_samples = struct();

% Initialize the wait bar
wait_bar = waitbar(0, 'Initializing...');

idx = 1;
while idx <= num_samples
    % Update the wait bar with the current progress
    waitbar(idx / num_samples, wait_bar, sprintf('Processing: %d of %d', idx, num_samples));

    fprintf ("Starting data aggregation for index %d\n", idx);
    % Randomly select values from the specified ranges or sets
    % num_users = randsample(users_range, 1);
    num_users = users_range;
    % number of antennas at all users
    Lxu = ones(1, num_users)*N_tx;
    % get_channel parameters
    channel_type = randsample(channel_types, 1);
    % N_rx = randsample(N_rx_range, 1);
    index = datasample(index_range, num_users);
    % get_fading parameters
    % dist = dist_range(1) + (dist_range(2) - dist_range(1)).*rand(1, num_users);
    dist = [3, 3, 3];
    % fft parameters
    % fft_length = randsample(fft_length_range, 1);
    fft_length = fft_length_range(1);
    % minPMAC parameters
    % bu_min = transpose(datasample(bu_min_single_user_single_tone_range, num_users));
    bu_min = [10;10;10];
    % w = w_single_user_range(1) + (w_single_user_range(2) - ...
    %     w_single_user_range(1)).*rand(num_users, 1);
    w = w_single_user_range(1).*rand(num_users, 1);

    % Print parameters used for this iteration
    fprintf ("Parameters: num_users: %d , N_rx: %d, fft_length: %d\n", ...
        num_users, N_rx, fft_length );
    disp ("bu_min");
    disp (bu_min);

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
        [fad_m, fad_std] = get_fading(channel_type, dist(user), fc);
        fad_db = fad_m + randn()*fad_std;
        h_user = h_user*10^(-fad_db/20);
        h(:, (user-1)*N_tx+1:user*N_tx, :) = h_user;
    end

    % FFT
    H = fft(h, fft_length, 3);
    
    fprintf ("Starting minPMACMIMO\n");
    % Applying minPMAC on generated channel with random values
    [FEAS_FLAG, bu_a, info] = minPMACMIMO(H/sqrt(noise_var), ...
        Lxu, fft_length*bu_min, w, cb);
    fprintf ("minPMACMIMO finished\n");

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
    data_samples(idx).bu_min = bu_min;
    data_samples(idx).w = w;
    data_samples(idx).h = h;
    data_samples(idx).H = H;
    % data_samples(idx).Eun = Eun;
    % data_samples(idx).theta = theta;
    % data_samples(idx).bun = bun;
    data_samples(idx).FEAS_FLAG = FEAS_FLAG;
    data_samples(idx).bu_a = bu_a;
    data_samples(idx).info = info;

    % save collected data cumulatively after every sample
    temp_filename = strcat('FA4\ML\outputs\data_samples', int2str(idx), '.mat');
    save(temp_filename, 'data_samples');

    idx = idx+1;
   
end

% Close the wait bar when the loop is done
close(wait_bar);

% Save the data samples to a .mat file
save('FA4\ML\outputs\data_samples.mat', 'data_samples');
