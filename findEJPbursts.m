function [ind_burststart,ind_burstend] = findEJPbursts(Vm,time,threshold,displayTimes)
% Finds and returns burst start and end intervals based on thresholding
%   INPUTS:
%       Vm: membrane potential trace (mV)
%       time: array of time points
%       threshold: mV above baseline used to find start and end
%       displayTimes: plot trace and start and end times
%   OUTPUTS:
%       ind_burststart: indices of burst starts
%       ind_burstend: indices of burst ends

sampling_freq = 1/time(2);

%% Find troughs of membrane potential (i.e. points in between EJP bursts)

% Low pass filter
fPass = 10; %filter frequency
wFilt = (2 * 1/sampling_freq) * fPass; [B,A] = butter(2, wFilt, 'low');
Vm_filt = filtfilt(B, A, Vm);

[amp_trough,ind_trough] = findpeaks(-Vm_filt,...
    'MinPeakDistance',sampling_freq/3,'MinPeakProminence',3);%threshold/1.5,'MinPeakHeight',mean(-Vm)-5,);

% Remove anything > 3 mV above minumum for each second

too_high = [];
for i = 1:length(ind_trough)
    if ind_trough(i) < sampling_freq
        indices = [1:sampling_freq];
    elseif ind_trough(i) > length(time) - sampling_freq
        indices = length(time) + [-sampling_freq:0];
    else
        indices = ind_trough(i) + [-sampling_freq/2:sampling_freq/2];
    end
    min_Vm = min(Vm(indices));
    if -amp_trough(i) > min_Vm + 3;
        too_high = [too_high,i];
    end
end
ind_trough(too_high) = [];
amp_trough(too_high) = [];

amp_diff = diff(-amp_trough);
ind_diff = diff(ind_trough);
not_decayed = find(ind_diff < sampling_freq/1.9 & amp_diff > 0.2)+1;
ind_trough(not_decayed) = [];
amp_trough(not_decayed) = [];

figure
hold on
plot(Vm,'b')
plot(Vm_filt,'k','LineWidth',2)
scatter(ind_trough,-amp_trough,'g')
%scatter(ind_trough(too_high),-amp_trough(too_high),'r')
%scatter(ind_trough(not_decayed),-amp_trough(not_decayed),'b')

% for i = 1:length(ind_trough)
%     if -amp_trough(i) < -amp_trough(i+1)
%     end
% end
 
% trough_bottomhalf = sort(ind_trough);
% trough_bottomhalf = trough_bottomhalf(round(length(ind_trough)/2):end);
% min_dist = mean(diff(trough_bottomhalf));
% 
% [amp_trough,ind_trough] = findpeaks(-Vm,...
%     'MinPeakProminence',5,'MinPeakHeight',mean(-Vm)-5,'MinPeakDistance',min_dist/4);

%% Find indices of putative burst starts and ends above a threshold

ind_burststart = NaN(size(ind_trough)); % save indices of burst starts
ind_burstend = NaN(size(ind_trough)); % save indices of burst ends

for i = 1:length(ind_trough)
    % define threshold
    thresh = amp_trough(i)-threshold; % X mV above trough
    
    % find index of times around trough when Vm first passes threshold
    try
        ind_burststart(i) = min(find(-Vm(ind_trough(i):end)<thresh))+ind_trough(i);
    end
    try
        ind_burstend(i) = max(find(-Vm(1:ind_trough(i))<thresh));
    end
    
    % find amplitude from baseline to threshold
%     if ind_trough(i) < sampling_freq
%         start = 1;
%         finish = ind_trough(i) + 2*sampling_freq;
%     elseif ind_trough(i) > length(time) - sampling_freq
%         start = ind_trough(i) - 2*sampling_freq;
%         finish = length(time);
%     else
%         start = ind_trough(i) - sampling_freq;
%         finish = ind_trough(i) + sampling_freq;
%     end
%     min_Vm = min(Vm(start:finish));
%     amp_burststart(i) = Vm(ind_burstend(i)) - min_Vm;
%     amp_burstend(i) = Vm(ind_burststart(i)) - min_Vm;
end

% figure
% hold on
% scatter(time(ind_burstend),amp_burststart,'g')
% scatter(time(ind_burststart),amp_burstend,'r')

% remove_burststart = find(amp_burststart > median(amp_burststart));
% remove_burstend = find(amp_burstend > median(amp_burstend));
% 
% ind_burststart(remove_burststart) = [];
% ind_burstend(remove_burstend) = [];

%% Remove non-unique/repeats and NaN of burst starts and ends

remove_burststart2 = find(diff(ind_burststart)==0);
remove_burstend2 = find(diff(ind_burstend)==0);
ind_burststart(remove_burststart2) = [];
ind_burstend(remove_burstend2) = [];

ind_burststart(isnan(ind_burststart)) = [];
ind_burstend(isnan(ind_burstend)) = [];

%% Remove burst starts where inter-start interval is too short

diff_burststart = diff(ind_burststart);
if std(diff_burststart) > 10000
    remove_burststart = find(diff_burststart < quantile(diff_burststart,0.75)/1.5)+1;
    ind_burststart(remove_burststart(1:end-1)) = [];
end

%% Remove burst ends where duration of bursts is too long (i.e. movement artifact)

durations = NaN(size(ind_burstend));

for i = 2:length(ind_burstend)
    start = max(find(ind_burststart < ind_burstend(i))); % burst start just prior to burst end
    durations(i) = ind_burstend(i) - ind_burststart(start);
end

if nanstd(durations) > 1000
    remove_burstend = find(durations > nanmedian(durations)*1.6);
    ind_burstend(remove_burstend) = [];
end

%% Remove burst ends where inter-end interval is too short

% diff_burstend = diff(ind_burstend);
% remove_burstend2 = find(diff_burstend < quantile(diff_burstend,0.75)/1.8);
% ind_burstend(remove_burstend2(1:end-1)) = [];

%% Remove burst multiple burst ends for same burst start

% remove_burstend3 = [];
% 
% for i = 2:length(ind_burststart)
%     ends = find(ind_burstend > ind_burststart(i-1) & ind_burstend < ind_burststart(i));
%     if length(ends) > 1
%         remove_burstend3 = [remove_burstend3,ends(1:end-1)'];
%     end
% end
% ind_burstend(remove_burstend3) = [];

remove_burststart3 = [];
for i = 2:length(ind_burstend)
    starts = find(ind_burststart > ind_burstend(i-1) & ind_burststart < ind_burstend(i));
    if length(starts) > 1
        remove_burststart3 = [remove_burststart3,starts(1:end-1)'];
    end
end
ind_burststart(remove_burststart3) = [];

remove_burstend3 = [];
for i = 2:length(ind_burststart)
    ends = find(ind_burstend > ind_burststart(i-1) & ind_burstend < ind_burststart(i));
    if length(ends) > 1
        remove_burstend3 = [remove_burstend3,ends(1:end-1)'];
    end
end
ind_burstend(remove_burstend3) = [];

%% Remove first or last EJP if end before start and vice versa

while ind_burststart(1) > ind_burstend(1)
    ind_burstend(1) = [];
end
while ind_burststart(end) > ind_burstend(end)
    ind_burststart(end) = [];
end

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

