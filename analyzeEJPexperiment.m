function [data] = analyzeEJPexperiment()
%Runs through an entire EJP experiment folder
%   INPUT:
%       directory: gui asks for folder to use
%   OUTPUT:
%       saves an mat file in the folder for all files

directory = uigetdir(); % get directory
files = dir([directory,'/*_*.abf']); % list files
filestoLoad = []; % find files to analyze
for i = 1:length(files)
    if isempty(strfind(files(i).name,'ABF1'))
        filestoLoad = [filestoLoad,{files(i).name}];
    end
end

if exist([directory,directory(end-7:end),'.mat'])
    data = load([directory,directory(end-7:end),'.mat']);
else
    data = struct();
end

for i = 1:length(filestoLoad)
    filename = filestoLoad{i};
    if ~isfield(data,['file',filename(9:12)])
        %analyzeEJPfile(filename,[],[],'off');
        analyzeEJPfile2(filename,[],'off');
    end
end

clear data
data = load([directory,directory(end-7:end),'.mat']);

for i = 1:length(filestoLoad)
    filename = filestoLoad{i};
    data.cumulative.meanBurstPeriod(i) = nanmean(data.(['file',filename(9:12)]).burstPeriod);
    data.cumulative.meanBurstDuration(i) = nanmean(data.(['file',filename(9:12)]).burstDuration);
    data.cumulative.meanriseTime(i) = nanmean(data.(['file',filename(9:12)]).riseTimes);
    data.cumulative.temperature(i) = data.(['file',filename(9:12)]).meanTemp;
    data.cumulative.EJPperBurst(i) = nanmean(data.(['file',filename(9:12)]).EJPperBurst);
end

%%

figure
subplot(2,2,1)
hold on
plot(data.cumulative.temperature(1:5),data.cumulative.meanBurstPeriod(1:5),'ko-')
plot(data.cumulative.temperature(6:10),data.cumulative.meanBurstPeriod(6:10),'bo-')
xlabel('Temperature (\circC)')
ylabel('Mean Burst Period (s)')
legend({'Control','Modulator'},'Location','northeast','box','off')

subplot(2,2,2)
hold on
plot(data.cumulative.temperature(1:5),data.cumulative.meanBurstDuration(1:5),'ko-')
plot(data.cumulative.temperature(6:10),data.cumulative.meanBurstDuration(6:10),'bo-')
xlabel('Temperature (\circC)')
ylabel('Mean Burst Duration (s)')

subplot(2,2,3)
hold on
plot(data.cumulative.temperature(1:5),data.cumulative.meanriseTime(1:5),'ko-')
plot(data.cumulative.temperature(6:10),data.cumulative.meanriseTime(6:10),'bo-')
xlabel('Temperature (\circC)')
ylabel('Mean Rise Constant (s)')

subplot(2,2,4)
hold on
plot(data.cumulative.temperature(1:5),data.cumulative.EJPperBurst(1:5),'ko-')
plot(data.cumulative.temperature(6:10),data.cumulative.EJPperBurst(6:10),'bo-')
xlabel('Temperature (\circC)')
ylabel('Mean EJPs per Burst')

savefig(gcf,[directory,directory(end-7:end)]);

end

