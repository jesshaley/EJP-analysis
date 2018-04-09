function [riseTime] = riseEJP(Vm,time,ind_startEJP,ind_peakEJP,displayRise)
% Finds EJP rise constant
%   INPUTS:
%       Vm: membrane potential trace (mV)
%       time: array of time points
%       ind_startEJP: indices of EJP starts
%       ind_peakEJP: indices of EJP peaks
%       displayRise: choose 'on' or 'off' to plot rise times
%   OUTPUTS:
%       riseTime: time from start to 63% peak (s)


time_startEJP = time(ind_startEJP);
Vm_startEJP = Vm(ind_startEJP);
Vm_peakEJP = Vm(ind_peakEJP);

Vm_rise63 = (1-(1/exp(1))) * (Vm_peakEJP-Vm_startEJP); % 63% change in amplitude of Vm
Vm_rise = Vm_rise63 + Vm_startEJP; % Vm at 63% of peak

for i = 1:length(time_startEJP)
    data_rise = Vm(ind_startEJP(i):ind_peakEJP(i)); % Vm during rise phase of EJP only
    ind_rise(i) = ind_startEJP(i) + min(find(data_rise > Vm_rise(i))); % index of time for 63% change
end
time_rise = time(ind_rise); % time relative to file start for 63% change

riseTime = time_rise - time_startEJP; % rise times for 63% change

if strcmp(displayRise,'on')
    figure
    subplot(1,2,1)
    scatter(time_startEJP,riseTime,'k')
    xlim([0 time(end)])
    xlabel('EJP Start Time (s)')
    ylabel('EJP Rise Time (63% Amp) (s)')
    
    subplot(1,2,2)
    scatter(riseTime,Vm_rise63,'k')
    xlabel('EJP Rise Time (63% Amp) (s)')
    ylabel('EJP 63% Amp (mV)')
end

end

