function [H] = get_channel_fwa(U, Lx, hx, hy, dist, fc, BW, fft_length,    ...
                                  scenario, manual_fading, transmission_type, ...
                                  Ain, Bin, Cin, Din, Ein, Fin, ...
                                  verbose)
    
    
    s = qd_simulation_parameters;
    s.use_absolute_delays = 1;
    s.center_frequency    = fc;
    s.show_progress_bars  = verbose;

    l = qd_layout(s);

    l.no_tx = 1;
    l.no_rx = U;

%     theta = [0, 5, 10, 15, 20, 25 30 35 40 45] * pi/180;
%     dist  = [dist dist dist dist dist dist dist dist dist dist];
%     
    theta = (0:5:5*(U-1))*pi/180;
    % dist = dist*ones(1, U);
    dist = [500, 600, 700];
    for i = 1:l.no_rx
        l.rx_position(:, i) = dist(i) * [cos(theta(i)), sin(theta(i)), hy];
    end

%     l.rx_position(:,1)  = [dist,  0,  hy] + 5  * [randn() randn() 0];
%     l.rx_position(:,2)  = [dist,  10, hy] + 5  * [randn() randn() 0];
%     l.rx_position(:,3)  = [dist, -10, hy] + 5  * [randn() randn() 0];
    

    l.tx_position(:,1) = [0  0 hx]';
    l.tx_array = qd_arrayant('3gpp-3d', Ain, Bin, Cin, Din, Ein, Fin);

    for i = 1:U
        l.rx_array(i) = qd_arrayant('omni');
        l.rx_array(i).center_frequency = fc;
        normalize_gain(l.rx_array(i), 1, 5);
    end

    for i = 1:Bin
       normalize_gain(l.tx_array, i, 13);
    end
            
    l.set_scenario(scenario); 
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
        
        
        
    % [fad_m, fad_std] = get_fading("FWA", dist, fc);
        
    h = zeros(l.no_rx, c(1).no_txant, fft_length);
        
        
    for i = 1:numel(c)
       split = strsplit(c(i).name, '_');
       r_idx = split(3);
       r_idx = r_idx{1};
       r_idx = str2num(r_idx(3:end));
       
       t_idx = split(2);
       t_idx = t_idx{1};
       t_idx = str2num(t_idx(3:end));
    
       % fad_db = fad_m + randn()*fad_std;
    
       
       if manual_fading == 1
           disp('Adding Manual Fading')
           h(r_idx, :, :) = c(i).fr(BW, fft_length);% * 10^(-fad_db/20); 
       else
           h(r_idx, :, :) = c(i).fr(BW, fft_length);
       end
    end

    
    if strcmp(transmission_type, 'MAC') == 1
        H = permute(h, [2 1 3]);    
    else
        H = h;
    end
end