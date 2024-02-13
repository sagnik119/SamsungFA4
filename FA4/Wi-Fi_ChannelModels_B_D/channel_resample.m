function resampled_signal = channel_resample(signal, fs_signal, L, M, n_taps, plt)
    fs_new = fs_signal*L/M;
    idx_signal = 0:numel(signal)-1;
    idx_resampled_signal = 0:n_taps-1;
    resampled_signal = zeros(1, numel(idx_resampled_signal));
    for m = 1:numel(idx_resampled_signal)
       temp = 0;
       for n = 1:numel(idx_signal)
          temp = temp + signal(n)* ...
              sinc(fs_signal*(idx_resampled_signal(m)/fs_new - idx_signal(n)/fs_signal));
       end
       resampled_signal(m) = temp;%*fs_signal/fs_new;%* temp; 
    end
    % resampled_signal = resampled_signal * sqrt(L/M);

    resampled_signal(idx_resampled_signal==0) = sum(resampled_signal(idx_resampled_signal <= 0));
    resampled_signal = resampled_signal(idx_resampled_signal>=0);
    resampled_signal(abs(resampled_signal) < 1e-8) = 0;
    idx_resampled_signal = idx_resampled_signal(idx_resampled_signal >=0);
    
    if plt
        figure;
        plot(linspace(-fs_signal/2, fs_signal/2, 4096), ...
             abs(fftshift(fft(signal, 4096))), 'displayname', '100MHz')
        hold on
        plot(linspace(-fs_new/2, fs_new/2, 4096), ...
             abs(fftshift(fft(resampled_signal, 4096))), 'displayname', [num2str(fs_new/1e6) 'MHz (resampled)'])

        grid minor
        title('Impulse Response (Frequency Domain)')
        xlabel('Frequency (Hz)', 'interpreter', 'latex')
        legend show 

        figure;
        stem(idx_signal/fs_signal, abs(signal), 'filled', 'displayname', '100MHz')

        hold on

        stem(idx_resampled_signal(1:end)/fs_new, abs(resampled_signal(1:end)), 'displayname', '80MHz (resampled)')
        grid on
        title('Impulse Response (Time Domain)')
        xlabel('time(s)', 'interpreter', 'latex')
        legend show    
    
    
    end
end

