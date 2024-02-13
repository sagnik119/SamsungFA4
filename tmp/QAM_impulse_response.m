% Parameters
N = 1000; % number of bits
M = 16; % 16-QAM
k = log2(M); % bits per symbol
rolloff = 0.25; % rolloff factor for the raised cosine filter
span = 6; % filter span
sps = 4; % samples per symbol

% Generate random bits
data = randi([0 1], N, 1);

% Map bits to 16-QAM symbols
symbols = bi2de(reshape(data, k, N/k).', 'left-msb');
qam_symbols = qammod(symbols, M, 'UnitAveragePower', true);

% Pulse shaping filter (Raised Cosine)
h = rcosdesign(rolloff, span, sps);

% Upsample and filter
up_sampled = upsample(qam_symbols, sps);
tx_signal = filter(h, 1, up_sampled);

% Plot impulse response of the raised cosine filter
figure;
stem(h, 'filled');
title('Impulse Response of Raised Cosine Filter');
xlabel('Sample');
ylabel('Amplitude');

% (Optional) Plot the 16-QAM constellation
% scatterplot(qam_symbols);                              