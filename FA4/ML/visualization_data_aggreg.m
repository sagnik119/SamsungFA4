
curr_ind = 0;
% Loop through each element of data_samples
for i = 1:length(data_samples)
    % Check if the number of users is 3
    if data_samples(i).num_users == 3
        % Get the tone size for this sample
        [~, tone_size] = size(data_samples(i).info.bun{1});
        tone_sizes = 1:tone_size;
        
        % Extract the bun values for each user and store them
        user1_bun = data_samples(i).info.bun{1}(1, 1:tone_size);
        user2_bun = data_samples(i).info.bun{1}(2, 1:tone_size);
        user3_bun = data_samples(i).info.bun{1}(3, 1:tone_size);
    end
    
    if curr_ind < 5
        % Plot the bun values for each user
        figure; % Create a new figure
        plot(tone_sizes, user1_bun, 'o-', 'DisplayName', 'User 1'); hold on;
        plot(tone_sizes, user2_bun, 'x-', 'DisplayName', 'User 2')
        plot(tone_sizes, user3_bun, 's-', 'DisplayName', 'User 3');
        
        
        % Add labels and legend
        xlabel('Tone');
        ylabel('bits/Hz');
        legend('show');
        title('Data Rates for Each User Across Different Tones');
    end 

    curr_ind = curr_ind+1;

end

