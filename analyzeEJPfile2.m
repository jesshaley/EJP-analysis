function [] = analyzeEJPfile2(filename,onlyAnalyze,removeArtifact)
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
elseif length(Vm)/sampling_freq > 300
    onlyAnalyze = [1 300*sampling_freq];
else
    onlyAnalyze = [1 length(Vm)];
end

Vm = Vm(onlyAnalyze(1):onlyAnalyze(2));
time_ms = time_ms(onlyAnalyze(1):onlyAnalyze(2));
temp = temp(onlyAnalyze(1):onlyAnalyze(2));

time = time_ms/1000; % time matrix in seconds
mean_temp = mean(temp); % average temperature in Celsius

%% Find starts and peaks of EJP

[dVm_dt,ind_startEJP,ind_peakEJP,Vm_filt] = findEJPpeaks(Vm,time);

%% Find movement artifacts

if strcmp(removeArtifact,'on')
    [artifact] = findEJPartifacts2(Vm,time,ind_startEJP,ind_peakEJP,'off');
else 
    [artifact] = findEJPartifacts3(Vm,time,ind_startEJP,ind_peakEJP,'off');
end
ind_startArtifact = ind_startEJP(artifact);
ind_startEJP(artifact) = [];
ind_peakArtifact = ind_peakEJP(artifact);
ind_peakEJP(artifact) = [];

%% Find dominant frequency using fft

nfft = length(Vm);
fourier = fft(detrend(Vm),nfft);
fourier = fourier(1:nfft/2+1);
increment = sampling_freq/nfft; % frequency increment
freqvec = 0:increment:sampling_freq/2;


fPass = 80; %filter frequency
wFilt = (2 * 1/sampling_freq) * fPass; [B,A] = butter(2, wFilt, 'low');
fourier(1:30) = 0;
fourier_filt = filtfilt(B, A, abs(fourier));
rang = (max(fourier_filt)-min(fourier_filt))/4;

[peaks,best_freq] = max(fourier);

% [peaks,best_freq] = findpeaks(fourier_filt,'MinPeakProminence',rang,...
%     'MinPeakHeight',max(fourier_filt)/3);

% figure
% hold on
% plot(freqvec,abs(fourier),'b')
% xlim([0 10])
% plot(freqvec,fourier_filt,'k')
% scatter(freqvec(best_freq),peaks)

best_freq = freqvec(best_freq);
% best_freq = best_freq(2)/2;

%% Find EJP burst duration and frequency

try
    [ind_burststart,ind_burstend,riseTime,decayTime,numEJPs] = findEJPbursts2(Vm,time,best_freq,ind_startEJP,ind_peakEJP,'off');
end

if ~exist('ind_burststart')
    ind_burststart = [];
    ind_burstend = [];
    riseTime = [];
    decayTime = [];
    numEJPs = [];
    burst_period = [];
    burst_freq = [];
    burst_duration = [];
else
    burst_period = diff(time(ind_burststart));
    burst_freq = 1/burst_period;
    burst_duration = ind_burstend-ind_burststart;
    burst_duration(isnan(burst_duration)) = 1;
    burst_duration = time(burst_duration);
    burst_duration(find(burst_duration == 0)) = NaN;
end

%% Plot Vm and dVm/dt

time_startEJP = time(ind_startEJP);
time_peakEJP = time(ind_peakEJP);
Vm_startEJP = Vm(ind_startEJP);
Vm_peakEJP = Vm(ind_peakEJP);

f = figure;
set(gcf,'Position',[50 300 1200 600])
subplot(2,1,1)
hold on
plot(time,Vm,'k','LineWidth',2)
scatter(time_startEJP,Vm_startEJP,'g')
scatter(time_peakEJP,Vm_peakEJP,'r')
scatter(time(ind_peakArtifact),Vm(ind_peakArtifact),'c')
ind_burststart2 = ind_burststart; ind_burststart2(isnan(ind_burstend)) = [];
ind_burstend2 = ind_burstend; ind_burstend2(isnan(ind_burstend)) = [];
for i = 1:length(ind_burststart2)
    plot([time(ind_burststart2(i)) time(ind_burstend2(i))],...
        [Vm(ind_burststart2(i)) Vm(ind_burstend2(i))],'b','LineWidth',2)
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
data.(['file',page_num]).indexArtifactStarts = ind_startArtifact;
data.(['file',page_num]).indexArtifactPeaks = ind_peakArtifact;

% analyzed data
data.(['file',page_num]).burstDuration = burst_duration';
data.(['file',page_num]).burstFrequency = burst_freq;
data.(['file',page_num]).burstPeriod = burst_period';
data.(['file',page_num]).riseTimes = riseTime'; %only for first EJP
data.(['file',page_num]).decayTimes = decayTime'; % only for last EJP
data.(['file',page_num]).EJPperBurst = numEJPs';

save([directory,filename(1:7)],'-struct','data'); % save structure
savefig(f,[directory,filename(1:12)]);
saveas(f,[directory,filename(1:12),'.png']);
end

