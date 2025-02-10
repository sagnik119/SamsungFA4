% Parameters
Ly                 = 2;
hx                 = 1.5;
hy                 = 3;
fc                 = 5e9;
BW                 = 80e6;
fft_length         = 64;
scenario           = 'WINNER_Indoor_A1_LOS';
manual_fading      = 0;
transmission_type  = 'MAC';
verbose            = 0;
dist               = 3;

noise_var          = 3e-21;
cb                 = 1;
Rnn                = repmat(eye(Ly), [1, 1, fft_length]);

transmit_SNR       = 10.^((0:10:80)/10);
Eun_level          = transmit_SNR * noise_var;
num_monte_carlo    = 100;
num_snr            = length(transmit_SNR);

U_values = 2:6;
max_U    = max(U_values);

% Pre-allocate result array: [U_scenarios, MonteCarlo, max_U, SNRs]
bu_a_lins = zeros(length(U_values), num_monte_carlo, max_U, num_snr);

inv_noise_sqrt = 1 / sqrt(noise_var);

% Suppose we know the total number of iterations:
total_iterations = length(U_values)*num_snr*num_monte_carlo;

% Create a waitbar
hWait = waitbar(0, 'Processing...');

% Create a DataQueue
q = parallel.pool.DataQueue;

% Define a callback function that runs each time we receive data from the workers
afterEach(q, @(x)updateProgress(x, total_iterations, hWait));

% Nested function to update progress
function updateProgress(~, total_iterations, hWait)
    persistent count; 
    if isempty(count)
        count = 0;
    end
    count = count + 1;
    fraction = count / total_iterations;
    waitbar(fraction, hWait, sprintf('Progress: %.2f%%', fraction*100));
end

for U_idx = 1:length(U_values)
    U = U_values(U_idx);
    Lxu        = ones(1, U);
    Eun        = eye(U);
    fft_factor = (fft_length/(fft_length + 8)) * sum(Eun, 2);
    
    for s = 1:num_snr
        fixed_gain = fft_factor * Eun_level(s);
        
        % Temporary cell array to hold results for this SNR
        temp_results_s = cell(num_monte_carlo, 1);
        
        % Parallelize over Monte Carlo iterations
        parfor i = 1:num_monte_carlo
            H = get_channel(U, Ly, 1, hx, hy, dist, fc, BW, fft_length, ...
                            scenario, manual_fading, transmission_type, verbose);
            H_scaled = H * inv_noise_sqrt;
            
            [~, ~, ~, bu_a_lin, ~] = linear_only(fixed_gain, H_scaled, Lxu, Rnn, cb);
            % bu_a_lin is of length U
            
            temp_results_s{i} = bu_a_lin; 

            % Send a progress update signal
            send(q, 1);
        end
        
        % Store the results for this SNR after parfor completes
        for i = 1:num_monte_carlo
            bu_a_lins(U_idx, i, 1:U, s) = temp_results_s{i};
        end
    end
end

close(hWait);

% After completion, bu_a_lins contains data organized as:
%   bu_a_lins(U_idx, i, :, s)
% corresponding to each U scenario, Monte Carlo run, each user (up to U), and each SNR.