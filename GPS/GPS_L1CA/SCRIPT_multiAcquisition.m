%% Clean up the environment first =========================================
clear all; %close all; clc;
tic

format ('compact');
format ('long', 'g');

%--- Include folders with functions ---------------------------------------
addpath include               % The software receiver functions
addpath Common                % Common functions between differnt SDR receivers

%% Print startup ==========================================================
% fprintf(['\n',...
%     'Welcome to:  softGNSS\n\n', ...
%     'An open source GNSS SDR software project initiated by:\n\n', ...
%     '              Danish GPS Center/Aalborg University\n\n', ...
%     'The code was improved by GNSS Laboratory/University of Colorado.\n\n',...
%     'The software receiver softGNSS comes with ABSOLUTELY NO WARRANTY;\n',...
%     'for details please read license details in the file license.txt. This\n',...
%     'is free software, and  you  are  welcome  to  redistribute  it under\n',...
%     'the terms described in the license.\n\n']);
% fprintf('                   -------------------------------\n\n');

%% Initialize constants, settings =========================================
settings = initSettings_EnacTP();

% Find number of samples per spreading code
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));

% Nb of ms to perform acquisition. Acq is performed on 2 consecutive slices
% of signal, each of length settings.acqNonCohTime*settings.acqCohTime
acqLen_ms = 2*(settings.acqNonCohTime*settings.acqCohTime)+2;

% Number of possible acquisition based on the signal duration
Nacq = floor(settings.msToProcess / acqLen_ms);
% limit Nacq to 100
Nacq = min([100 Nacq]);

%% Generate plot of raw data and ask if ready to start processing =========
disp(['File: ' settings.fileName]);

%% Initialization =========================================================
disp ('Starting processing...');
[fid, message] = fopen(settings.fileName, 'rb');

%Initialize the multiplier to adjust for the data type
if settings.fileType==1, dataAdaptCoeff=1;
else,                    dataAdaptCoeff=2;
end

%If success, then process the data
if (fid > 0)
    for indAcq = 1:Nacq;
       
        
        %% Acquisition ============================================================
        % Read data for acquisition.
        data  = fread(fid, dataAdaptCoeff*acqLen_ms*samplesPerCode, settings.dataType)';
        
        if (dataAdaptCoeff==2)
            data1=data(1:2:end);
            data2=data(2:2:end);
            data=data1 + 1i .* data2;
        end
        
        %--- Do the acquisition -------------------------------------------
        disp(['Acquisition #' num2str(indAcq,'%03i') '/' num2str(Nacq,'%03i') ' - Nb of Coherent Sums: ' num2str(settings.acqCohTime) ' - Nb of Non-Coherent Sums: ' num2str(settings.acqNonCohTime)]);
        acqResults(indAcq) = acquisition_EnacTP(data, settings);

    end
    
    fclose(fid);
end

%% compute acquisition performances
% delay vector, in chips
Fs = settings.samplingFreq;
Lc = settings.codeLength;
Fc = settings.codeFreqBasis;
NsamplesPerCodePeriod = Lc/Fc*Fs;
delay_v = (0:NsamplesPerCodePeriod-1)/Fs*Fc;
% find index of true delay and doppler in name
trueDelay_tc = str2double(settings.fileName(strfind(settings.fileName,'tau=')+4:strfind(settings.fileName,'Tc')-1));
trueDelay_idx = round(trueDelay_tc/settings.codeFreqBasis*settings.samplingFreq)+1;
trueDop_hz = str2double(settings.fileName(strfind(settings.fileName,'dop=')+4:strfind(settings.fileName,'Hz')-1));
trueDop_idx = find(settings.freqBinList == trueDop_hz);
% delay index vector without main peak (for noise level computation)
delayIdxNoPeak = setdiff((1:NsamplesPerCodePeriod),trueDelay_idx + (-round(Fs/Fc):round(Fs/Fc)));

prn_v = settings.acqSatelliteList;
for indPrn = 1:length(settings.acqSatelliteList)
    currentPrn = prn_v(indPrn);
    for indAcq = 1:Nacq
        % detection outcome (1 or 0)
        detection(indPrn,indAcq) = acqResults(indAcq).peakMetric(currentPrn) > settings.acqThreshold;
        % peak height at true delay/doppler
        try % work only for simulated signals
            peak_height(indPrn,indAcq)  = acqResults(indAcq).acqMat(currentPrn,trueDop_idx,trueDelay_idx);
            noise_height(indPrn,indAcq) = mean(acqResults(indAcq).acqMat(currentPrn,trueDop_idx,delayIdxNoPeak));
        end
    end
    % display detection rate
    fprintf("PRN %02i - Detection rate = %.2f\n",currentPrn,sum(detection(indPrn,:))/Nacq)
end; clear indPrn indAcq

figure;
bar(settings.acqSatelliteList,sum(detection,2)/Nacq);
xlabel('PRN')
ylabel('Detection rate');
title(settings.fileName);
grid on

try % works only for simulated signals
    fig = figure;
    subplot(2,1,1); hold all;
    plot(peak_height');
    yline(mean(peak_height,2));
    legend('PRN present','PRN absent');
    ylabel('Peak height');
    subplot(2,1,2); hold all;
    plot(noise_height');
    yline(mean(noise_height,2));
    ylabel('Noise height');
    xlabel('Acquisition index');
catch
    close(fig)
end



toc