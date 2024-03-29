% function [En, bn] = waterfill(total_en, gn, gap, cb)% 
% Standard Lagrange Multiplier waterfill
% Water accepts the channel gains as input, and finds water-fill solution
% for the given (input) total energy, with codes of gap "gap."
% The gap is specified as linear (so not in dB). 
% cb = 1 for complex channel; cb = 2 for real baseband channel
% The program uses any gain set and does not compute an SNR, nor does it
% know the original of the gains (so no guard-band penalty is presumed)
%
% Inputs 
%       total_en is the total energy for all dimensions
%       gn is a vector containing all the gains (energy xfer)
%       gap is thelinear gap of the code, so gap =1 means 0 dB
%       cb indicates whether the channel is complex baseband (cb=1) or real
%       baseband (cb = 2)
% Output
%       En is the set of energies for the original dimensions
%       bn is the set of rates for the original dimensions

function [En, bn] = waterfill(total_en, gn, gap, cb)
    Ntot = numel(gn);
    Ex_bar = total_en/Ntot;
    [gn_sorted, Index]=sort(gn, 'descend');  % sort gain, and get Index
    [~,InvIndex] = sort(Index);
    % Ntot = numel(gn);
    num_zero_gn = length(find(gn_sorted == 0)); %number of zero gain subchannels
    Nstar= Ntot - num_zero_gn;
    % Number of used channels,
    % start from Ntot - (number of zero gain subchannels)
    
    while(1)
        K=1/Nstar*(Ntot*Ex_bar+gap*sum(1./gn_sorted(1:Nstar)));
        En_min=K-gap/gn_sorted(Nstar);	% En_min occurs in the worst channel
        if (En_min<0)
            Nstar=Nstar-1;  % If negative En, continue with less channels
        else
            break;       % If all En positive, done.
        end
    end
    
    En=K-gap./gn_sorted(1:Nstar); 		% Calculate En
    bn = 1/cb * log2(1 + gn_sorted(1:Nstar).*En);
    En = [En; zeros(Ntot-Nstar, 1)];        % Pad zeros to revert to unsorted order
    bn = [bn; zeros(Ntot-Nstar, 1)];
    En = En(InvIndex);
    bn = bn(InvIndex);
end