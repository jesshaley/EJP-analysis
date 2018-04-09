function [dVm_dt_smooth,ind_startEJP,ind_peakEJP,Vm_filt] = findEJPpeaks(Vm,time)
% Finds EJP start and peak indices
%   INPUTS:
%       Vm: membrane potential trace (mV)
%       time: array of time points
%   OUTPUTS:
%       dVm_dt: first derivative of Vm
%       ind_startEJP: indices of EJP starts
%       ind_peakEJP: indices of EJP peaks

sampling_freq = 1/time(2); % sampling frequency (Hz)

%% Calculate derivate of Vm

% Low pass filter
fPass = 30; %filter frequency
wFilt = (2 * 1/sampling_freq) * fPass; [B,A] = butter(2, wFilt, 'low');
Vm_filt = filtfilt(B, A, Vm);

dVm_dt_filt = diff(Vm_filt)./diff(time);

dVm_dt = diff(Vm)./diff(time);
dVm_dt_smooth = smooth(dVm_dt,15); % smooth with moving window of 15

%% Check for small EJPs

% check if EJPs are small (amp < 5 mV); if yes, use Vm if no, use dVm/dt
if max(Vm)-min(Vm) < 10
    small = 'yes';
else
    small = 'no';
end

%% Find peaks using Vm

if strcmp(small,'yes')
    [Vm_peakEJP,ind_peakEJP] = findpeaks(Vm,'MinPeakProminence',0.5,'MinPeakDistance',sampling_freq/10);
    for i = 1:length(ind_peakEJP)
        [Vm_startEJP(i),ind_startEJP(i)] = min(Vm([-600:0]+ind_peakEJP(i)));
        ind_startEJP(i) = ind_startEJP(i) - 600 + ind_peakEJP(i);
        %Vm_EJP(i) = Vm_peakEJP(i)-Vm_startEJP(i);
    end
end

%% Find Peaks of dVm/dt

if ~strcmp(small,'yes')
    [dVm_dt_pkamp,dVm_dt_pkind] = findpeaks(dVm_dt_filt,...
        'MinPeakProminence',75,'MinPeakHeight',40,'MinPeakDistance',sampling_freq/1000);
    
    % Find EJP start and peak indices
    
    dVm_dt_0 = find(abs(dVm_dt_filt) < 15); % find indices where dVm/dt ~= 0
    ind_startEJP = NaN(size(dVm_dt_pkind))'; ind_peakEJP = ind_startEJP;
    
    % find indices where dVm/dt ~= 0 just before and right after peaks of dVm/dt
    for i = 1:length(dVm_dt_pkind)
        try
            ind_startEJP(i) = dVm_dt_0(max(find(dVm_dt_0 < dVm_dt_pkind(i)))); % start index of EJP
            ind_peakEJP(i) = dVm_dt_0(min(find(dVm_dt_0 > dVm_dt_pkind(i)))); % start index of EJP
            
            % find closest point to 0 within 10 ms jitter
            [min_0,startoffset(i)] = min(abs(dVm_dt_smooth(ind_startEJP(i)+[-sampling_freq/100:sampling_freq/100])));
            [max_0,peakoffset(i)] = min(abs(dVm_dt_smooth(ind_peakEJP(i)+[-sampling_freq/100:sampling_freq/100])));
        end
    end
    
    % reset start and peak indices to reflect exactly where dVm_dt inflects
    ind_startEJP = ind_startEJP - sampling_freq/100 + startoffset;
    ind_peakEJP = ind_peakEJP - sampling_freq/100 + peakoffset;
end

% figure
% subplot(2,1,1)
% hold on
% plot(Vm)
% scatter(ind_startEJP,Vm(ind_startEJP),'g')
% scatter(ind_peakEJP,Vm(ind_peakEJP),'r')
% 
% subplot(2,1,2)
% hold on
% plot(dVm_dt_smooth)
% scatter(ind_startEJP,dVm_dt_smooth(ind_startEJP),'g')
% scatter(ind_peakEJP,dVm_dt_smooth(ind_peakEJP),'r')

%% Remove duplicates and NaN

% find and remove duplicates due to deviations in rise of EJP
remove = [find(isnan(ind_startEJP)),find(diff(ind_startEJP)<150)+1];
ind_startEJP(remove) = [];
ind_peakEJP(remove) = [];

end

