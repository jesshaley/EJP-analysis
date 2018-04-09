function [ind_burststart,ind_burstend,riseTime,decayTime,numEJPs] = findEJPbursts2(Vm,time,freq,ind_startEJP,ind_peakEJP,displayTimes)
% Finds and returns burst start and end intervals based on thresholding
%   INPUTS:
%       Vm: membrane potential trace (mV)
%       time: array of time points
%       threshold: mV above baseline used to find start and end
%       displayTimes: plot trace and start and end times
%   OUTPUTS:
%       ind_burststart: indices of burst starts
%       ind_burstend: indices of burst ends

sampling_freq = 1/time(2); % sampling frequency (Hz)

%% Find index of first and last EJP in putative bursts

Vm_diff = diff(Vm(ind_startEJP)');
ind_diff = diff(ind_startEJP);
not_decayed = find(ind_diff < sampling_freq/(freq*1.9) & Vm_diff > 0.2 | ...
        (ind_diff < sampling_freq/(freq*2.3) & ind_diff < sampling_freq/5))+1;
firstEJP = [1:length(ind_startEJP)]; firstEJP(not_decayed) = [];
lastEJP = [firstEJP(2:end)-1,length(ind_startEJP)];

%% Threshold is 67% rise and decay of first and last EJP

Vm_rise63 = (1-(1/exp(1))) * (Vm(ind_peakEJP(firstEJP))-Vm(ind_startEJP(firstEJP))); % 63% change in amplitude of Vm
Vm_decay63 = (1-(1/exp(1))) * (Vm(ind_peakEJP(lastEJP))-Vm(ind_startEJP(firstEJP)));

Vm_rise = Vm_rise63 + Vm(ind_startEJP(firstEJP));
Vm_decay = Vm(ind_peakEJP(lastEJP)) - Vm_decay63;

ind_burststart = NaN(size(firstEJP)); ind_burstend = ind_burststart;

for i = 1:length(firstEJP)
    data_rise = Vm(ind_startEJP(firstEJP(i)):ind_peakEJP(firstEJP(i))); % Vm during rise phase of EJP only
    if i == length(firstEJP)
        data_decay = Vm(ind_peakEJP(lastEJP(i)):length(Vm)); % Vm during decay phase of EJP only
    else
        data_decay = Vm(ind_peakEJP(lastEJP(i)):ind_startEJP(firstEJP(i+1))); % Vm during decay phase of EJP only
    end
    try
        ind_burststart(i) = ind_startEJP(firstEJP(i)) + min(find(data_rise > Vm_rise(i))); % index of time for 63% change
        ind_burstend(i) = ind_peakEJP(lastEJP(i)) + min(find(data_decay < Vm_decay(i))); % index of time for 63% change
    end
end

%% Remove first or last EJP if end before start and vice versa

while ind_burststart(1) > ind_burstend(1)
    ind_burstend(1) = [];
end
while ind_burststart(end) > ind_burstend(end)
    ind_burststart(end) = [];
end

%% Calculate Rise and Decay Constants

riseTime = ind_burststart-ind_startEJP(firstEJP);
decayTime = ind_burstend-ind_peakEJP(lastEJP);

%% Calculate EJPs per burst

numEJPs = diff([firstEJP,length(ind_startEJP)+1]);

%% Plot Burst start and end times

if strcmp(displayTimes,'on')
    figure
    set(gcf,'Position',[50 300 1200 300])
    hold on
    plot(time,Vm,'k','LineWidth',2)
    %scatter(time(ind_trough),-amp_trough,'b')
    scatter(time(ind_burststart),Vm(ind_burststart),'g')
    scatter(time(ind_burstend),Vm(ind_burstend),'r')
    %scatter(time(ind_burststart(remove_burststart)),Vm(ind_burststart(remove_burststart)),'c')
    %scatter(time(ind_burstend(remove_burstend3)),Vm(ind_burstend(remove_burstend3)),'c')
    xlim([0 max(time)])
    xlabel('Time')
    ylabel('V_m (mV)')
end

end

