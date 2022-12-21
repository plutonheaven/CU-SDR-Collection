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

%% Initialization =========================================================
[fid, message] = fopen(settings.fileName, 'rb');

%Initialize the multiplier to adjust for the data type
if settings.fileType==1, dataAdaptCoeff=1;
else,                    dataAdaptCoeff=2;
end

%If success, then process the data
if (fid > 0)
    %% Acquisition ============================================================
    disp(['File:                                  ' settings.fileName]);
    disp(['Acquisition - Nb of Coherent Sums:     ' num2str(settings.acqCohTime)]);
    disp(['Acquisition - Nb of Non-Coherent Sums: ' num2str(settings.acqNonCohTime)]);

    % Read data for acquisition.
    data  = fread(fid, dataAdaptCoeff*acqLen_ms*samplesPerCode, settings.dataType)';
    
    if (dataAdaptCoeff==2)
        data1=data(1:2:end);
        data2=data(2:2:end);
        data=data1 + 1i .* data2;
    end
    
    %--- Do the acquisition -------------------------------------------
    acqResults = acquisition_EnacTP(data, settings);
end

fclose(fid);

disp ('   Ploting results...');
if settings.plotAcquisition
    plotAcquisition(acqResults);
%     plotAcquisitionMat(acqResults,settings,settings.acqSatelliteList);
end
disp('Post processing of the signal is over.');
toc