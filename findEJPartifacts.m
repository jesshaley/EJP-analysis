function [artifact] = findEJPartifacts(Vm,time,ind_startEJP,ind_peakEJP,displayArtifact)
% Finds and returns movement artifacts in EJP recordings
%   INPUTS:
%       Vm: membrane potential trace (mV)
%       time: array of time points
%       ind_startEJP: indices of EJP starts
%       ind_peakEJP: indices of EJP peaks
%       displayArtifact: choose 'on' or 'off' to plot artifacts
%   OUTPUTS:
%       ind_artifact: indices of movement artifacts

time_startEJP = time(ind_startEJP);
time_peakEJP = time(ind_peakEJP);
Vm_startEJP = Vm(ind_startEJP);
Vm_peakEJP = Vm(ind_peakEJP);

Vm_EJP = Vm_peakEJP - Vm_startEJP;
time_EJP = time_peakEJP - time_startEJP;

%% Find first and last EJPs

Vm_diff = diff(Vm(ind_startEJP)');
ind_diff = diff(ind_startEJP);
not_decayed = find(ind_diff < sampling_freq/1.9 & abs(Vm_diff) > 0.2 | ind_diff < sampling_freq/5)+1;
firstEJP = [1:length(ind_startEJP)]; firstEJP(not_decayed) = [];
lastEJP = [firstEJP(2:end)-1,length(ind_startEJP)];

%%

artifact = [];
for i = 2:length(ind_peakEJP)-1
    if ((Vm_peakEJP(i) < min(Vm) + 10)... % peak amp < 10 mV above baseline
            || time_peakEJP(i) - time_startEJP(i) > 0.28... % rise time is > 0.4 s
            || time_startEJP(i) - time_startEJP(i-1) > 0.3)... % > 0.3 s lag between last EJP
            && time_peakEJP(i)+1.5 < time_startEJP(i+1)... % >1.5 seconds from peak to next EJP
            && abs(Vm_EJP(i) - quantile(Vm_EJP,0.75)) > std(Vm_EJP) %% significantly larger or smaller than 4th quantile
        artifact = [artifact,i];
    end
end

%| Vm_peakEJP(i) < mean(Vm(ind_peakEJP(i)+300:ind_peakEJP(i)+400)... % Vm is still rising after peak

if strcmp(displayArtifact,'on')
    figure
    subplot(2,2,1)
    hold on
    scatter(time_startEJP,Vm_EJP,'k')
    scatter(time_startEJP(artifact),Vm_EJP(artifact),'b')
    xlim([0 time(end)])
    xlabel('EJP Start Time (s)')
    ylabel('EJP Amplitude (mV)')
    
    subplot(2,2,3)
    hold on
    scatter(time_startEJP,time_EJP,'k')
    scatter(time_startEJP(artifact),time_EJP(artifact),'b')
    xlim([0 time(end)])
    xlabel('EJP Start Time (s)')
    ylabel('EJP Time from Start to Peak (s)')
    
    subplot(2,2,[2,4])
    hold on
    scatter(time_EJP,Vm_EJP,'k')
    s_art = scatter(time_EJP(artifact),Vm_EJP(artifact),'b')
    xlabel('EJP Time from Start to Peak (s)')
    ylabel('EJP Amplitude (mV)')
    legend([s_art],{'Artifact'},'Box','off','Location','Northwest')
end

end

