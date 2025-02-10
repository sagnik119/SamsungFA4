function [H] = get_channel(U, Ly, Lxu, hx, hy, dist, fc, BW, fft_length, scenario, manual_fading, transmission_type, verbose)
    
    addpath('quadriga_src/');

    s = qd_simulation_parameters;
    s.use_absolute_delays = 1;
    s.center_frequency    = fc;
    s.show_progress_bars  = verbose;

    l = qd_layout(s);
    l.no_tx = 1;
    l.no_rx = U;
    

    
    h = zeros(l.no_rx, Ly, fft_length);
    
%     theta = [0, 5, 10, 15, 20, 25 30 35 40 45] * pi/180;
%     dist  = [3, 3, 3, 3, 3, 3, 3, 3, 3, 3];
    theta = (0:5:5*(U-1))*pi/180;
    dist = dist*ones(1, U);

    for antenna = 1:Ly
        l.tx_position(:,1) = [0  0  hy]';
        
        for i = 1:l.no_rx
            l.rx_position(:, i) = dist(i) * [cos(theta(i)), sin(theta(i)), hx];
            if Lxu == 1
                l.rx_array(i) = qd_arrayant('omni');
            else
                l.rx_array(i) = qd_arrayant(['ula', num2str(Lxu)]);
            end
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
        
        
        
        
        [fad_m, fad_std] = get_fading("B", dist, fc);
        
        
        
        for i = 1:numel(c)
           split = strsplit(c(i).name, '_');
           r_idx = split(3);
           r_idx = r_idx{1};
           r_idx = str2num(r_idx(3:end));
           
           t_idx = split(2);
           t_idx = t_idx{1};
           t_idx = str2num(t_idx(3:end));
        
           fad_db = fad_m(r_idx) + randn()*fad_std;
        
           
           if manual_fading == 1
               disp('Adding Manual Fading')
               h(Lxu*r_idx-(Lxu-1):Lxu*r_idx, antenna, :) = c(i).fr(BW, fft_length) * 10^(-fad_db/20); 
           else
               h(Lxu*r_idx-(Lxu-1):Lxu*r_idx, antenna, :) = c(i).fr(BW, fft_length);
           end
        end
    end

    if strcmp(transmission_type, 'MAC') == 1
        H = permute(h, [2 1 3]);    
    end

end
