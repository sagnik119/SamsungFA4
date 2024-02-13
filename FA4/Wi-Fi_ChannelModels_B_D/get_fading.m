function [fad_mean_db, fad_std_db] = get_fading(model_letter, dist, fc)
%get_fading get log-normal distribution of path loss + shadow fading loss
%(in dB). To generate fading samples, use fad_mean_db + fad_std_db*randn().
%   Inputs:
%   - model_letter: channel model letter, 'B' or 'D'
%   - dist: distance between Tx and Rx in meter
%   - fc: carrier frequency
%   Outputs:
%   - fad_mean_db: mean of the fading (in dB)
%   - fad_std_db: standard deviation of the log-normal fading (in dB)


LFSd = @(d, f) 20*log10(d) + 20*log10(f) - 147.5;
if model_letter == 'B'
    if dist <= 5
        fad_mean_db = LFSd(dist, fc);
        fad_std_db = 3;
    else
        fad_mean_db = LFSd(5, fc) + 35*log10(dist/5);
        fad_std_db = 4;
    end
elseif model_letter == 'D'
    if dist <= 10
        fad_mean_db = LFSd(dist, fc);
        fad_std_db = 3;
    else
        fad_mean_db = LFSd(10, fc) + 35*log10(dist/10);
        fad_std_db = 5;
    end
else
    disp("Not Implemented");
end
end