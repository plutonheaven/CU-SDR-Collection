function plotAcquisitionMat(acqResults,settings,prn_v)
%Functions plots acquisition matrix
%
%plotAcquisition(acqResults)
%
%   Inputs:
%       acqResults    - Acquisition results from function acquisition.

%% Plot all results =======================================================
% delay vector, in chips
Fs = settings.samplingFreq;
Lc = settings.codeLength;
Fc = settings.codeFreqBasis;
NsamplesPerCodePeriod = Lc/Fc*Fs;
delay_v = (0:NsamplesPerCodePeriod-1)/Fs*Fc;
for indPrn = 1:length(prn_v)
    results = squeeze(acqResults.acqMat(prn_v(indPrn),:,1:NsamplesPerCodePeriod));
    [~, frequencyBinIndex] = max(max(results, [], 2));
    
    % compute SNR
%     Nsize = numel(results);
%     % Number of points corresponding to the main peak : (2.Tc.Fs2) * (2/Tint*F_dop_grid)
%     Npeak = ceil(((2+0.5)/Fc*Fs) * ((2+0.5)/1e-3/(settings.freqBinList(2) - settings.freqBinList(1))));
%     corr_sort = sort(reshape(results,1,Nsize));
%     % var_noise = mean(corr_sort(1:end-Npeak))/(2*settings.acqNonCohTime);
%     peak_height = max(max(results));
%     noise_height = mean(corr_sort(1:end-Npeak));
    
    [peak_height,delayBinIndex] = max(results(frequencyBinIndex,:));
    % mean height of the corr matrix, on the correct frequency bin, and
    % accross all delays except +/- 1 chip around peak
    noise_height = mean(results(frequencyBinIndex,setdiff((1:NsamplesPerCodePeriod),delayBinIndex + (-round(Fs/Fc):round(Fs/Fc)))));
    SNR_dB = 10*log10( peak_height / noise_height);
    fprintf("PRN %02i - estimated SNR = %2.1f dB",prn_v(indPrn),SNR_dB)
    if acqResults.peakMetric(indPrn) < settings.acqThreshold, fprintf(" (sat not acquired)"); end
    fprintf("\n")
    
    figure(101+indPrn);
    surf(delay_v,settings.freqBinList/1e3,results);
    shading interp
    xlabel('Delay (chip)')
    ylabel('Doppler frequency (kHz)')
end; clear indPrn
    