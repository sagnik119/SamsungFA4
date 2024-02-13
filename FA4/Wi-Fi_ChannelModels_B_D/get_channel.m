% ---------------------------------------------------------------------------------------------------------------------
% COMMERCIAL IN CONFIDENCE
% Copyright (c) 2010, CSR Ltd. All rights reserved.
% ---------------------------------------------------------------------------------------------------------------------
%
% function [h, tap_delay] = get_channel(model_letter, N_tx, N_rx, index)
%
%
% Inputs:
% ---------------------------------------------------------------------------------------------------------------------
%
% model_letter: Channel model letter, 'A' to 'F'
% N_tx        : Number of transmit antennas
% N_rx        : Number of receive antennas
% index       : index of the channel, 1 to 5000
%
% Outputs:
% ---------------------------------------------------------------------------------------------------------------------
%
% h           : The channel matrix, N_rx by N_tx by N_taps
% tap_delay   : the tap delays in s
% tap_delay_10ns   : the tap delays in 10ns
%

function [h, tap_delay tap_delay_10ns] = get_channel(model_letter, N_tx, N_rx, index, T_s)
persistent current_model tap_delays H

if ~exist('T_s', 'var')
    T_s = 50;
end
muap = length(index);

% remove trailing tilde
tilde_posn = strfind(model_letter, '~');
if tilde_posn
     model_letter = model_letter(1:tilde_posn-1);
end
if length(model_letter) == 2 % LLS channel
    model = sprintf('IEEE_CHANNEL_%sTX4RX4_L%g', model_letter(1), T_s);
    if ~strcmp(current_model, model)
        load (model);
        fileName = ['IEEE_CHANNEL_' model_letter(1) 'TX4RX4_40000.bin'];
        tmp = fopen(fileName);
        H = fread(tmp,'float');
        fclose(tmp);                
        current_model = model;
    end        
    tap_delay = tap_delays;
    if tilde_posn
        % return the even spaced tap delays used previously for backwards compatibility 
        tap_delay_10ns = 1:18;
    else
        tap_delay_10ns = round(tap_delays/10e-9);
    end

    h = zeros(N_rx,N_tx,length(tap_delay));
    for i_rx = 1:N_rx
        for i_STS = 1:N_tx/muap
            for i_ap = 1:muap
                chIdx=mod(index(i_ap)-1,5000)+1;
                pos = 2*(4*40000*(i_STS-1)+40000*(i_rx-1) + (chIdx-1))*length(tap_delay)+1;
                for i = 1:length(tap_delay)
                    h(i_rx,(i_ap-1)*N_tx/muap+i_STS,i) = H(pos+2*(i-1)) + 1i*H(pos+2*(i-1)+1);
                end
            end
        end
    end
else
    if length(model_letter) == 1 
        if (N_tx > 8) || (N_rx > 8)
            model = sprintf('IEEE_CHANNEL_%sTX16RX4_L%g', model_letter, T_s);
        elseif (N_tx > 4) || (N_rx > 4)
            model = sprintf('IEEE_CHANNEL_%sTX8RX8_L%g', model_letter, T_s);
        else
            model = sprintf('IEEE_CHANNEL_%sTX4RX4_L%g', model_letter, T_s);
        end
        if ~strcmp(current_model, model)
            load (model);
            disp(sprintf('loaded %s', model))
            H = H_nrx_ntx;
            current_model = model;
        end
    elseif (length(model_letter) == 3) && (strcmp(model_letter(1:2), 'tv')) && (~isempty(strfind('ABCDEF', model_letter(3))))
        T_s = 2600;
        model = sprintf('IEEE_CHANNEL_%sNLoSTX4RX4_L%g_CohT1000.mat', model_letter(3), T_s);
        if ~strcmp(current_model, model)
            load (model);
            disp(sprintf('loaded %s', model))
            H = H_nrx_ntx;
            current_model = model;
        end        
    elseif length(model_letter) == 3 && (model_letter(2) == 'L' || model_letter(2) == 'N')
        %% For NS3 testing
        channelFlag = model_letter(1);
        losFlag = model_letter(2);
        if strcmp(losFlag,'L')
            LosString = 'LoS';
        else
            LosString = 'NLoS';
        end
        ChannelString = [channelFlag LosString];
        if N_tx <= 4 && channelFlag == 'D'
            model = sprintf(['IEEE_CHANNEL_5GHz_' ChannelString '_TX4RX4_L13.dat']);    
        else
            model = sprintf(['IEEE_CHANNEL_5GHz_' ChannelString '_TX16RX4_L13.dat']);    
        end
        if ~strcmp(current_model, model)
            H = load (model);
            current_model = model;
            tap_delays = get_delays(model_letter)*10e-9;
        end
    end

    tap_delay = tap_delays;
    if tilde_posn
        % return the even spaced tap delays used previously for backwards compatibility 
        tap_delay_10ns = 1:18;
    else
        tap_delay_10ns = round(tap_delays/10e-9);
    end
    %h = squeeze(H(1:N_rx,1:N_tx,:,index));
    if length(model_letter) == 1 
        sizeH = size(H);
        h = zeros(N_rx,N_tx,size(H,3));
        for i_ap = 1:muap
            chIdx=mod(index(i_ap)-1,sizeH(end))+1;
            if sizeH(2) < N_tx/muap
                tmp = H((i_ap-1)*N_tx/muap+1:i_ap*N_tx/muap,1:N_rx,:,chIdx);
                h(:,(i_ap-1)*N_tx/muap+1:i_ap*N_tx/muap,:) = permute(tmp,[2 1 3]);
            else
                h(:,(i_ap-1)*N_tx/muap+1:i_ap*N_tx/muap,:) = H(1:N_rx,(i_ap-1)*N_tx/muap+1:i_ap*N_tx/muap,:,chIdx);
            end
        end

    elseif (length(model_letter) == 3) && (strcmp(model_letter(1:2), 'tv')) && (~isempty(strfind('ABCDEF', model_letter(3))))
        [~, ~, ~, len_h_tap] = size(H); 
        coh_length = floor(len_h_tap/ 1000);
        len_5ms = ceil(5e-03 / (1/2600));
        h = zeros(N_rx,N_tx,size(H,3));
        for i_ap = 1:muap
            index_range = [((index(i_ap)-1) * coh_length + 1) : ((index(i_ap)-1) * coh_length + len_5ms)];
            h(:,(i_ap-1)*N_tx/muap+1:i_ap*N_tx/muap,:) = H(1:N_rx,(i_ap-1)*N_tx/muap+1:i_ap*N_tx/muap,:,index_range);
        end
    elseif length(model_letter) == 3 && (model_letter(2) == 'L' || model_letter(2) == 'N')
        h = zeros(N_rx,N_tx,length(tap_delay_10ns));
        if N_tx<=4
            refNumTx = 4;
        else
            refNumTx = 16;
        end
        for rxIdx = 1:N_rx
            for txIdx = 1:N_tx/muap
                for i_ap = 1:muap
                    Chval = H(refNumTx*(rxIdx-1)+((i_ap-1)*N_tx/muap+txIdx-1)+index(i_ap),:);
                    h(rxIdx,(i_ap-1)*N_tx/muap+txIdx,:) = Chval(1:2:end)+1i*Chval(2:2:end);
                end
            end
        end
    end
                
end



