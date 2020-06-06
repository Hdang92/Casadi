clear all 
close all


addpath(genpath('/home/haidang/Matlab/LA_toolbox'))
addpath('/home/haidang/Matlab/casadi-linux-matlabR2014b-v3.4.5')
load('cycles_eyes_closed.mat')


for k=1:4
    resp=data(k).BS;
    stim=data(k).FS;
    [FD(k),TD(k)] = getFRF(stim,resp,100);
    loglog(FD(k).f,FD(k).Remnants)
    hold on
end

    xlabel('frequency (Hz)')
    ylabel('amplitude (deg^2)')


% function [x,y] = createnoise(butter1,noise1,cut1,cut2,q)
%     
%     for k=1:50
%         samp_freq = 100; % [Hz] Abtastrate der AD-Wandlerbox;
%         fco = [cut1,cut2]; % [Hz] Cut off Frequenzen fuer den bandpass filter
%         duration = 60.5; % [s] Stimulationsdauer pro noise amplitude
%         noise = wgn((duration)*samp_freq,1,noise1); % white noise Signal
%         fnyq = samp_freq/q;
%         [b,a] = butter(3,fco/fnyq,'bandpass'); % bandpass filter
%         noise = filter(b,a,noise);
% 
%         n_out(k,:) = noise*3;
%     end
%         [~,y,x] = getSpec(n_out,100);
%         [~,yoo,f] = getSpec(n_out,100);
% end
% 
% f = fit(createnoise
% [x,y] = creatnoise(coeff(1),coeff(2),coeff(3),coeff(4));
% loglog(x,mean(y))
% xlim([0.01 2])
hold on
for k=1:50
     samp_freq = 100; % [Hz] Abtastrate der AD-Wandlerbox;
        fco = [0.003, 0.9]; % [Hz] Cut off Frequenzen fuer den bandpass filter
        duration = 60.50; % [s] Stimulationsdauer pro noise amplitude
        noise = wgn((duration)*samp_freq,1,0); % white noise Signal
        fnyq = samp_freq/2;
        [b,a] = butter(3,fco/fnyq,'bandpass'); % bandpass filter
        noise = filter(b,a,noise);
        n_out(k,:) = noise*3;
end
[~,yoo,f] = getSpec(n_out,100);
yoo= yoo(1:62);
f=f(1:62);
loglog(f,mean(yoo))
xlim([0.01 2])