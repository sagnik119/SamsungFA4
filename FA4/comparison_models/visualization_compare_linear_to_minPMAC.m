% To compare the performance between minPMAC GDFE solution and linear
% (current WiFi architecture) solution
load('C:\GitHub\Samsung_FA4_Stanford_2023\FA4\comparison_models\outputs\Data_SIMO','data_samples');

curr_ind = 0;
% Loop through each element of data_samples
for i = 1:length(data_samples)
    % Check if the number of users is 3
    if data_samples(i).num_users == 3
        % Get the tone size for this sample
        [~, tone_size] = size(data_samples(i).bun);
        tone_sizes = 1:tone_size;

        % Extract the bun values for each user and store them
        user1_bun = data_samples(i).bun(1, 1:tone_size);
        user2_bun = data_samples(i).bun(2, 1:tone_size);
        user3_bun = data_samples(i).bun(3, 1:tone_size);
        user1_bun_lin = data_samples(i).bun_lin(1, 1:tone_size);
        user2_bun_lin = data_samples(i).bun_lin(2, 1:tone_size);
        user3_bun_lin = data_samples(i).bun_lin(3, 1:tone_size);
    end

    % if curr_ind < 1
        % Plot the max possible bun values for each user
        figure; % Create a new figure
        % subplot (1, 9, [1 2 3 4]);
        plot(tone_sizes, user1_bun, 'o-', 'DisplayName', 'User 1'); hold on;
        plot(tone_sizes, user2_bun, 'x-', 'DisplayName', 'User 2')
        plot(tone_sizes, user3_bun, 's-', 'DisplayName', 'User 3');


        % Add labels and legend
        xlabel('Tone');
        ylabel('bits/Hz');
        legend('show');
        title('Max Data Rates for Each User Across Different Tones');

    %     % Plot the max possible linear bun values for each user
    %     subplot (1, 9, [5 6 7 8]);
    %     plot(tone_sizes, user1_bun_lin, 'o-', 'DisplayName', 'User 1'); hold on;
    %     plot(tone_sizes, user2_bun_lin, 'x-', 'DisplayName', 'User 2')
    %     plot(tone_sizes, user3_bun_lin, 's-', 'DisplayName', 'User 3');
    % 
    % 
    %     % Add labels and legend
    %     xlabel('Tone');
    %     ylabel('bits/Hz');
    %     legend('show');
    %     title('Data Rates using Linear Receiver for Each User Across Different Tones');
    % 
    %     subplot (1, 9, 9);
    %     axis off;
    %     str = sprintf(['N_{tx}: %d\nN_{rx}: %d\ndist(m): %0.2f, %0.2f, %0.2f\n' ...
    %         'fft length: %d\nbu_{min}(Mbps): %d, %d, %d\n'], ...
    %         data_samples(i).N_tx, ...
    %         data_samples(i).N_rx, ...
    %         data_samples(i).dist(1, 1), data_samples(i).dist(1, 2), data_samples(i).dist(1, 3), ...
    %         data_samples(i).fft_length, ...
    %         data_samples(i).bu_min(1), ...
    %         data_samples(i).bu_min(2), ...
    %         data_samples(i).bu_min(2));
    %     annotation('textbox', [0.8, 0.1, 0.28, 0.8], 'String', str, 'EdgeColor', ...
    %         'none', 'HorizontalAlignment', 'center');
    % end 

    % curr_ind = curr_ind+1;

end

% Statistics across all collected samples
% counts = zeros(1,9);
% avg_rates_max = zeros(1,9);
% avg_rates_minPMAC = zeros(1,9);
% avg_rates_lin = zeros(1,9);
% for i = 1:length(data_samples)
%     curr_fft_length = data_samples(i).fft_length;
%     counter = int32(log2(curr_fft_length) -2);
%     avg_rates_max(i) = data_samples(i).b_sum;
%     % avg_rates_minPMAC(i) = sum(data_samples(i).bu_a);
%     avg_rates_minPMAC(i) = sum(data_samples(i).bun, 'all');
%     avg_rates_lin(i) = data_samples(i).b_sum_lin;
% end
% 
% figure; % Create a new figure
% plot(2.^(3:11), avg_rates_max, 'o-', 'DisplayName', 'Max Sum-rate'); hold on;
% plot(2.^(3:11), avg_rates_minPMAC  , '-', 'DisplayName', 'MinPMAC Sum-rate'); hold on;
% plot(2.^(3:11), avg_rates_lin, 'x-', 'DisplayName', 'Linear Sum-rate')
% 
% % Add labels and legend
% xlabel('FFT length');
% ylabel('Sum Rate per tone');
% legend('show');
% title('Max versus linear data rates for users across different FFT sizes');


