clc
clear
close all

Nmax = 16;
Ly = 2;
h=cat(3,[1 0 .8 ; 0 1 1],[.9 -.3 0 ; .5 -1 -1],[0 .2 0 ; .4 -.63 0],[0 0 0 ; 0 .648 0])*10;
bsum=zeros(1, 2);
bsumlin=zeros(1, 2);
% bsum_minpmac = zeros (1, 2);
bu_min = [11, 11, 11];
w = [1, 1, 1];
cb = 1;
for index=2
i=8*index; % (don8t need to plot a point for every number of tones)
H = fft(h, i, 3);
Rnn=zeros(Ly,Ly,i);
for n=1:i
    pass
Rnn(:,:,n) = eye(2);
end
[Rxx, bsum(index), bsumlin(index)] = SWF(i/(i+3)*[1 1 1], H, [1 1 1], Rnn(:,:,:), 1);
% bsum(index)=bsum(index)/(i+3);
% bsumlin(index)=bsumlin(index)/(i+3);
[Eun, theta, bun, FEAS_FLAG, bu_a, info] = minPMAC(H, i*bu_min', w', cb);
% bsum_minpmac(index) = sum(bu_a, 'all');
end
% plot([8, 16], bsum,[8, 16],bsumlin, bsum,[8, 16],bsum_minpmac);