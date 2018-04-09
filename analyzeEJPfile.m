function [] = analyzeEJPfile(filename,onlyAnalyze,threshold,fitDecay)
% Finds EJP characteristics for files at a single temperature
%   INPUTS:
%       filename: e.g. '900_133_0007.abf'
%       onlyAnalyze: if part of the file is too noisy, put in min and max
%           time to analyze in seconds (e.g. [0 200]); if blank, whole file
%           is analyzed
%       threshold: used for finding burst start and end times (mV); if 
%           blank, automatic values are used
%       fitDecay: not working properly, say 'off' or leave empty for now
%   OUTPUTS:
%       saves an mat file with all values

%% Load file from the specified path

file_path = which(filename);
directory = file_path(1:end-length(filename));
if exist([directory,filename(1:7),'.mat'])>0
    data = load([directory,filename(1:7),'.mat']);
end

abf = LoadAbf(filename);
Vm = abf.data.VM1_2; % membrane potential of muscle fibre
temp = abf.data.Temp; % temperature of bath
time_ms = abf.time'; % time matrix in ms
sampling_freq = 1000/time_ms(2); % sampling frequency of recording

if ~isempty(onlyAnalyze)
    onlyAnalyze = onlyAnalyze*sampling_freq; % convert time to indices
    onlyAnalyze(1) = onlyAnalyze(1)+1;
    Vm = Vm(onlyAnalyze(1):onlyAnalyze(2));
    time_ms = time_ms(onlyAnalyze(1)+1:onlyAnalyze(2));
    temp = temp(onlyAnalyze(1)+1:onlyAnalyze(2));
else
    onlyAnalyze = [1 length(Vm)];
end

time = time_ms/1000; % time matrix in seconds
mean_temp = mean(temp); % average temperature in Celsius

%% Find starts and peaks of EJP

[dVm_dt,ind_startEJP,ind_peakEJP] = findEJPpeaks(Vm,time);
    
time_startEJP = time(ind_startEJP);
time_peakEJP = time(ind_peakEJP);
Vm_startEJP = Vm(ind_startEJP);
Vm_peakEJP = Vm(ind_peakEJP);

%% Find movement artifacts

[ind_artifact] = findEJPartifacts(Vm,time,ind_startEJP,ind_peakEJP,'off');

%% Find EJP burst duration and frequency

Vm_EJP = Vm_peakEJP - Vm_startEJP;

if isempty(threshold)
    if max(Vm) - min(Vm) < 10
        threshold = quantile(Vm_EJP,0.25)/2;
    else
        threshold = 7;
    end
end

try
    [ind_burststart,ind_burstend] = findEJPbursts(Vm,time,threshold,'off');
end

if ~exist('ind_burststart')
    ind_burststart = [];
    ind_burstend = [];
    burst_period = [];
    burst_freq = [];
    burst_duration = [];
else
    burst_period = diff(time(ind_burststart));
    burst_freq = 1/burst_period;
    burst_duration = time(ind_burstend) - time(ind_burststart);
end

%%

% figure
% subplot(1,2,1)
% title('Burst Frequency')
% hist(burst_freq)
% xlabel('Frequency (Hz)')
% ylabel('# of Occurences')
% 
% subplot(1,2,2)
% title('Burst Duration')
% hist(burst_duration)
% xlabel('Duration (s)')
% ylabel('# of Occurences')

%% Plot Vm and dVm/dt

f = figure;
set(gcf,'Position',[50 300 1200 600])
subplot(2,1,1)
hold on
plot(time,Vm,'k','LineWidth',2)
scatter(time_startEJP,Vm_startEJP,'g')
scatter(time_peakEJP,Vm_peakEJP,'r')
scatter(time_peakEJP(ind_artifact),Vm_peakEJP(ind_artifact),'b')
for i = 1:length(ind_burststart)
    plot([time(ind_burststart(i)) time(ind_burstend(i))],...
        [Vm(ind_burststart(i)) Vm(ind_burstend(i))],'c','LineWidth',2)
end
xlim([0 max(time)])
xlabel('Time (s)')
ylabel('V_m (mV)')
t = title([filename(1:12),' EJP Recordings @ temp = ',num2str(mean_temp),'C']);
set(t,'interpreter', 'none')
set(gca,'FontSize',16)

subplot(2,1,2)
hold on
plot(time(1:end-1),dVm_dt,'k','LineWidth',2)
scatter(time(ind_startEJP),dVm_dt(ind_startEJP),'g')
scatter(time(ind_peakEJP),dVm_dt(ind_peakEJP),'r')
xlim([0 max(time)])
xlabel('Time (s)')
ylabel('dV_m/dt (mV/s)')
set(gca,'FontSize',16)

%% Calculate the Rise and Decay Constants

[riseTime] = riseEJP(Vm,time,ind_startEJP,ind_peakEJP,'off');

if isempty(fitDecay)
    fitDecay = 'off';
end

if strcmp(fitDecay,'on')
    [decayTime,decay_gof] = decayEJP(Vm,time,ind_startEJP,ind_peakEJP,'on','off');
end

%% Save Results and Figure

page_num = filename(9:12);

% raw data
data.(['file',page_num]).Vm = Vm';
data.(['file',page_num]).time = time';
data.(['file',page_num]).temp = temp';
data.(['file',page_num]).sampleFrequency = sampling_freq;
data.(['file',page_num]).meanTemp = mean_temp;
data.(['file',page_num]).fileName = filename;
data.(['file',page_num]).indexAnalyzed = onlyAnalyze;

% processed data
data.(['file',page_num]).indexBurstStart = ind_burststart';
data.(['file',page_num]).indexBurstEnd = ind_burstend';
data.(['file',page_num]).indexEJPStarts = ind_startEJP;
data.(['file',page_num]).indexEJPPeaks = ind_peakEJP;
data.(['file',page_num]).indexEJPArtifacts = ind_artifact;

% analyzed data
data.(['file',page_num]).burstDuration = burst_duration';
data.(['file',page_num]).burstFrequency = burst_freq;
data.(['file',page_num]).burstPeriod = burst_period';
data.(['file',page_num]).riseTimes = riseTime';
try
    data.(['file',page_num]).decayTimes = decayTime';
    data.(['file',page_num]).decayGoFit = decay_gof';
end

save([directory,filename(1:7)],'-struct','data'); % save structure
savefig(f,[directory,filename(1:12)]);
saveas(f,[directory,filename(1:12),'.png']);
end

