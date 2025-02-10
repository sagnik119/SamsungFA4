% Assuming the following variables are in the workspace:
% bu_a_lins: [length(U_values), num_monte_carlo, max_U, num_snr]
% U_values
% num_monte_carlo
% transmit_SNR
% max_U

SNR_dB = 10*log10(transmit_SNR)-50;

% Define line properties
highlight_line_width = 2;
other_line_width = 0.5;

% Light gray for non-highlighted lines
other_color = [0.7 0.7 0.7];

% Generate a set of distinct colors for each U value
% lines(n) returns an n-by-3 matrix of distinct colors
colors = lines(length(U_values));

num_U = length(U_values);

for U_idx_highlight = 1:num_U
    figure; hold on; grid on;
    
    % Plot all U_values
    for U_idx = 1:num_U
        U = U_values(U_idx);

        % Compute mean_min_R for this U
        min_across_users = squeeze(min(bu_a_lins(U_idx, :, 1:U, :), [], 3));
        mean_min_R = mean(min_across_users, 1);

        % If U_idx <= U_idx_highlight, this line is highlighted
        if U_idx <= U_idx_highlight
            plot(SNR_dB, mean_min_R, '-o', ...
                'LineWidth', highlight_line_width, ...
                'Color', colors(U_idx, :));
        else
            % Non-highlighted line (lighter, thinner)
            plot(SNR_dB, mean_min_R, '-o', ...
                'LineWidth', other_line_width, ...
                'Color', other_color);
        end
    end

    xlabel('SNR (dB)');
    ylabel('Average Min(R)');
    title_str = sprintf('Average Min(R) vs SNR');
    title(title_str);

    % Legend entries: Show each U with its respective color in the legend
    legend_entries = arrayfun(@(x) sprintf('U = %d', x), U_values, 'UniformOutput', false);
    
    % Create a legend with the correct colors:
    % We must ensure legend colors match the plotted lines.
    % The lines are plotted in order, so the legend picks them up in the same order.
    legend(legend_entries, 'Location', 'best');
    hold off;
end