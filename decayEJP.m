function [decayTime,decay_gof] = decayEJP(Vm,time,ind_startEJP,ind_peakEJP,displayDecay,displayFit)
% Fits EJP decay with exponential function and finds decay constant
%   INPUTS:
%       Vm: membrane potential trace (mV)
%       time: array of time points
%       ind_startEJP: indices of EJP starts
%       ind_peakEJP: indices of EJP peaks
%       displayDecay: choose 'on' or 'off' to plot decay times
%       displayFit: choose 'on' or 'off' to plot exponential decay
%       functions (adds a lot of time to function)
%   OUTPUTS:
%       decayTime: time from peak to 63% decayed (s)
%       decay_gof: r-squared of goodness of fit for decay function

sampling_freq = 1/time(2); % sampling frequency
time_peakEJP = time(ind_peakEJP);
time_startEJP = time(ind_startEJP);
Vm_peakEJP = Vm(ind_peakEJP);

Vm_decay63 = NaN(size(ind_startEJP)); % 63% of amplitude (mV)
Vm_decay = NaN(size(ind_startEJP)); % Vm at 63% decayed (mV)
ind_decay = NaN(size(ind_startEJP)); % index at 63% decayed
time_decay = NaN(size(ind_startEJP)); % time at 63% decayed (s)
decayTime = NaN(size(ind_startEJP)); % decay time (s)
decay_gof = NaN(size(ind_startEJP)); % decay goodness of fit

for i = 1:length(ind_startEJP)-1
    % Define indices for start and end of data to fit
    start = ind_peakEJP(i) + 5; % begin fit 5 indices after EJP peak
    if i == length(ind_startEJP)
        finish = length(time); % if last EJP, use end of file
        baseline = min(Vm(ind_peakEJP(i):end)); %baseline after EJP
    else
        finish = ind_startEJP(i+1); % end fit before next EJP start
        baseline = min(Vm([0:sampling_freq]+ind_peakEJP(i))); % baseline within 1s after EJP
    end
    [min_val,min_index] = min(Vm(start:finish)); % find min Vm between EJP peak and next EJP start
    finish = start + min_index; % make min Vm the end fit
    fit_ind = [start:finish]'; % indices of data to be fit
    
    % Define Vm for data to fit
    fit_Vm = Vm(fit_ind); % Vm of data to be fit
    
    % Fit data, but remove start and offset
    [decay_fun,decay_gofs] = fit(fit_ind-start,fit_Vm-baseline,'exp1'); % General model Exp1: f(a,b,x) = a*exp(b*x)
    decay_gof(i) = decay_gofs.rsquare;
    
    if decay_gof(i) > 0.7 & length(fit_Vm) > 10 & fit_Vm(1) - fit_Vm(end) > 2
        % Extrapolate fit to find decay
        fit_beyond = [0:sampling_freq]'; % find decay constant during next second
        fit_decay = decay_fun(fit_beyond)+baseline; % fitted Vm during EJP decay
        
        % Calculate Vm and index at 63% of decay from EJP peak to baseline
        Vm_decay63(i) = (1-(1/exp(1))) * (Vm_peakEJP(i)-baseline); % 63% change in amplitude of Vm
        Vm_decay(i) = Vm_peakEJP(i) - Vm_decay63(i); % Vm at 63% of decay
        ind_decay(i) = ind_startEJP(i) + max(find(fit_decay < Vm_decay(i))); % index of time for 63% change
        
        % Calculte decay time
        time_decay(i) = time(ind_decay(i)); % time relative to file start for 63% change
        decayTime(i) =  time_peakEJP(i) - time_decay(i); % decay times for 63% change
        
        if strcmp(displayFit,'on')
            figure(1)
            subplot(2,1,1)
            plot(time(fit_beyond+start),fit_decay,'c')
        end
    end
end

% why does last index not work? need to be more careful about start and end
% of file issues

figure
subplot(1,2,1)
scatter(time_startEJP,decayTime,'k')
xlim([0 time(end)])
xlabel('EJP Start Time (s)')
ylabel('EJP Decay Time (63% Amp) (s)')

subplot(1,2,2)
scatter(decayTime,Vm_decay63,'k')
xlabel('EJP Decay Time (63% Amp) (s)')
ylabel('EJP Amplitude (mV)')

end

