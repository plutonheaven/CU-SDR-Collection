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
settings = initSettings();

% Find number of samples per spreading code
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));

% At least 42ms of signal are needed for fine frequency estimation
% % % codeLen = max(42,2*(settings.acqNonCohTime*settings.acqCohTime)+2);
codeLen = 2*(settings.acqNonCohTime*settings.acqCohTime)+2;

% Number of possible acquisition based on the signal duration
Nacq = 100; floor(settings.msToProcess / (settings.acqCohTime*settings.acqNonCohTime));

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
        
% % %         % Move the starting point of processing. Can be used to start the
% % %         % signal processing at any point in the data record (e.g. good for long
% % %         % records or for signal processing in blocks).
% % %         fseek(fid, dataAdaptCoeff*settings.skipNumberOfBytes, 'bof');
        
        %% Acquisition ============================================================
        % Read data for acquisition.
        data  = fread(fid, dataAdaptCoeff*codeLen*samplesPerCode, settings.dataType)';
% % %         fprintf("%i bytes read from file...\n", ftell(fid))
        
        if (dataAdaptCoeff==2)
            data1=data(1:2:end);
            data2=data(2:2:end);
            data=data1 + 1i .* data2;
        end
        
        %--- Do the acquisition -------------------------------------------
        disp(['Acquisition #' num2str(indAcq,'%03i') ' - Nb of Coherent Sums: ' num2str(settings.acqCohTime) ' - Nb of Non-Coherent Sums: ' num2str(settings.acqNonCohTime)]);
        acqResults(indAcq) = acquisition(data, settings);

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

for indPrn = 1:length(settings.acqSatelliteList)
    for indAcq = 1:Nacq
        % detection outcome (1 or 0)
        detection(indPrn,indAcq) = acqResults(indAcq).peakMetric(indPrn) > settings.acqThreshold;
        % peak height at true delay/doppler
        peak_height(indPrn,indAcq)  = acqResults(indAcq).acqMat(indPrn,trueDop_idx,trueDelay_idx);
        noise_height(indPrn,indAcq) = mean(acqResults(indAcq).acqMat(indPrn,trueDop_idx,delayIdxNoPeak));
    end
    % display detection performances
    switch indPrn
        case 1 % satellite present
            fprintf("PRN %02i - Detection rate = %.2f\n",indPrn,sum(detection(indPrn,:))/Nacq)
        case 2
            fprintf("PRN %02i - False Alarm rate = %.2f\n",indPrn,sum(detection(indPrn,:))/Nacq)
    end
end; clear indPrn indAcq

% figure; plot(detection');
% sum(detection,2)/Nacq;
figure; plot(peak_height');
figure; plot(noise_height');


toc