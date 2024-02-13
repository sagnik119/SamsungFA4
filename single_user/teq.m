function [P,ryy,rxy,rle,v,d,w,b,mmse,SNRmfb_in_dB,c]=teq(p,L,mu,delta,sigma,Ex_bar,filtertype)
% programmed by Ryoulhee Kwak
%
% filtertype 1=FIR , 2= IIR(just for one pole filter)
% p=channel impulse response 
% in case of FIR, p= 1+0.5D+2.3D^2->p=[1 0.5 2.3] 
% in case of IIR,  p=1/(2.4-0.5D)-> p=[2.4 -0.5]
% sigma =>noise power in linear scale
% delta =>delay
% L => number of taps
% c=> w*p
% v=>eigen  vectors
% d=>eigent values
%mu => for b

if filtertype==1% FIR 

    norm_p=norm(p);
    size_p=size(p, 2);
    P=toeplitz([p(1), zeros(1,L-1)]',[p,zeros(1,L-1)]);
    
    ryy=Ex_bar*(P*P')+sigma*eye(L,L);
    rxy=[zeros(mu+1,delta) eye(mu+1) zeros(mu+1,L+size_p-2-mu-delta)]*P'*Ex_bar;
    rle=eye(mu+1)*Ex_bar-rxy/ryy*rxy';
    [v, d]=eig(rle);
    d_temp=diag(d)';
    [m,n ]= min(d_temp);
    b=norm_p*v(:,n)';
    w=b*rxy/ryy;
    conv(w,p);
    mmse=b*rle*b';
%by using easy way to get error energy
    c=conv(w,p)
    alpa=c(delta+1)./b(1)
    biased_error_energy=norm_p^2*(m-Ex_bar*(1-alpa)^2)
    unbiased_error_energy=norm_p^2*(m-Ex_bar*(1-alpa)^2)/alpa^2

    SNRmfb=norm_p^2*Ex_bar/unbiased_error_energy;
    SNRmfb_in_dB=10*log10(SNRmfb)

else %IIR
    %in case of IIR filter P matrix is infinite dimesion but there is a trick to get exact rxy, ryy
    norm_p=sqrt(1/(p(1)^2*(1-(p(2)/p(1))^2 )));
    size_p=size(p,2);
    ptemp=[(-p(2)/p(1)).^(0:1:L-1)]/p(1); %-1 factor!
    P=toeplitz([ptemp(1) zeros(1,L-1)]',ptemp);
    Ptemp=toeplitz(ptemp',ptemp);

    ryy=Ex_bar*Ptemp*norm_p^2+sigma*eye(L,L);
    rxy=[zeros(mu+1,delta) eye(mu+1) zeros(mu+1,L-1-mu-delta)]*P'*Ex_bar;
    rle=eye(mu+1)*Ex_bar-rxy/ryy*rxy';
    [v, d]=eig(rle);
    d_temp=diag(d)';
    [m,n ]= min(d_temp);
    b=norm_p*v(:,n)';
    w=b*rxy/ryy;
    c=conv(w,p);
    sum(conv(w,p).^2)-norm(b)^2*(w(1)/b(1))^2;
    mmse=b*rle*b';
 %by using easy way to get error energy
    alpa=c(delta+1)./b(1)
    biased_error_energy=norm_p^2*(m-Ex_bar*(1-alpa)^2)
    unbiased_error_energy=norm_p^2*(m-Ex_bar*(1-alpa)^2)/alpa^2

    SNRmfb=norm_p^2*Ex_bar/unbiased_error_energy;
    SNRmfb_in_dB=10*log10(SNRmfb);
end