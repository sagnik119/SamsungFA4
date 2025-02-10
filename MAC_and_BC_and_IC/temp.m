% Test Script for minPMAC_reformulated and minPMAC_WMMSE

% Clear workspace and command window
clear; clc;

% Seed the random number generator for reproducibility
rng(1);

% Define system parameters
U = 4;          % Number of users
N = 1;          % Number of time slots (can be increased as needed)
Ly = 2;         % Number of receive antennas at the base station
Lx = 1;         % Number of transmit antennas per user (assuming single antenna)

% Generate random channel matrices H (Ly x U x N)
H = zeros(Ly, U, N);
for n = 1:N
    H(:, :, n) = sqrt(1/2) * (randn(Ly, U) + 1i * randn(Ly, U));
end

% Define minimum rate requirements bu_min (1 x U)
bu_min = [1, 1, 1, 1];  % Adjust as needed (in bits per channel use)

% Define weights w (1 x U)
w = ones(1, U);  % Equal weights for all users

% Define codebook scaling factor cb
cb = 1;  % Adjust as needed

% Run minPMAC_reformulated
disp('Running minPMAC_reformulated...');
tic;
[FEAS_FLAG_reformulated, bu_a_reformulated, info_reformulated] = minPMAC_reformulated(H, bu_min, w, cb);
time_reformulated = toc;
disp(['Time taken by minPMAC_reformulated: ', num2str(time_reformulated), ' seconds']);

% Run minPMAC_WMMSE
disp('Running minPMAC_WMMSE...');
tic;
[FEAS_FLAG_WMMSE, bu_a_WMMSE, info_WMMSE] = minPMAC_WMMSE(H, bu_min, w, cb);
time_WMMSE = toc;
disp(['Time taken by minPMAC_WMMSE: ', num2str(time_WMMSE), ' seconds']);

% Compare FEAS_FLAG
if FEAS_FLAG_reformulated == FEAS_FLAG_WMMSE
    disp('Both methods have the same feasibility flag.');
else
    disp('Methods have different feasibility flags.');
end

% Compare achieved rates bu_a
tolerance = 1e-3;  % Tolerance for comparing achieved rates
rate_difference = norm(bu_a_reformulated - bu_a_WMMSE, Inf);
if rate_difference < tolerance
    disp('Achieved rates are approximately the same for both methods.');
else
    disp('Achieved rates are different between the methods.');
    disp(['Maximum rate difference: ', num2str(rate_difference)]);
end

% Print achieved rates for both methods
disp('Achieved rates from minPMAC_reformulated:');
disp(bu_a_reformulated);
disp('Achieved rates from minPMAC_WMMSE:');
disp(bu_a_WMMSE);