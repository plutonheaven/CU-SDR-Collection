clear all;
close all;
addpath include               % The software receiver functions
addpath Common                % Common functions between differnt SDR receivers

tic

% simulated signal characteristics
prn = 1; % between 1 and 32
delay_true_chips = 511; % in chips, between 0 and 1022
doppler_true_Hz = 1500; % in Hz, typ. between +/- 7 kHz
signal_power_dBm = -140;-130; % in dBm. Typical open sky power is -128.5 dBm

% read settings (mainly to define GPS L1C/A parameters
settings = initSettings();
    
% change signal duration - max duration used in Practical Work is 2*11*11 = 242 ms.
settings.msToProcess = 500;

% Number of Nav message bits
Nbits = settings.msToProcess/20;
% Draw random bits
bits = sign(randn(1,25));
% Determine sign of the correlation function for each code period
codePeriodSign = kron(bits,ones(1,20));

% convert delay in samples
delay_true_samples = round(delay_true_chips/settings.codeFreqBasis*settings.samplingFreq); % in samples
% compute SNR and noise power
noise_psd = -203.9; % dBW/Hz. Typical value for RF front-end with integrated LNA
CN0_dB = signal_power_dBm - 30 - noise_psd;
noise_power = 10^(-CN0_dB/10)*settings.samplingFreq;
% sampling frequency
ts = 1 / settings.samplingFreq;
% number of samples in one code period (1 ms)
N_1ms = 1e-3*settings.samplingFreq; % number of samples in 1 ms
% number of code periods
Nperiod = settings.msToProcess/1e3/(settings.codeLength/settings.codeFreqBasis);
% filename
filename = ['simulatedSignal_tau=' num2str(delay_true_chips) 'Tc_' ...
            'dop=' num2str(doppler_true_Hz) 'Hz_' ...
            'pow=' num2str(signal_power_dBm) 'dBm' ...
            '.bin'];

fprintf("Generating %2i ms of GPS L1 C/A signal for PRN #%02i at Fs = %g MHz and power %3.1f dBm\n", ...
        settings.msToProcess, prn, settings.samplingFreq/1e6,signal_power_dBm);
% initial phase set to 0 radians
phi0 = 0;
% open file
fid = fopen(filename,'w');
for indPeriod = 1:Nperiod
    % create 1ms of GPS L1 C/A signal
    prnCode_1ms = makeCaTable(prn,settings);
    
    % multiply by nav message bit sign
    prnCode_sign = codePeriodSign(indPeriod)*prnCode_1ms;
    
    % introduce delay by circularly shifting the code period by delay
    prnCode_delay = circshift(prnCode_1ms,delay_true_samples);
    
    % introduce doppler
    prnCode_doppler = prnCode_delay.*exp(1j*2*pi*(doppler_true_Hz*(0:N_1ms-1)/settings.samplingFreq)+phi0);
    phi0 = phi0 + N_1ms*doppler_true_Hz/settings.samplingFreq;
    
    % add noise
    noise = sqrt(noise_power/2)*(randn(1,N_1ms)+1j*randn(1,N_1ms));
    prnCode_noise = prnCode_doppler + noise;
    
    % write signal to file (appends value to fid)
    saveSamplesToFile(fid,prnCode_noise,settings.fileType,settings.dataType);
end; clear indPeriod
% close file
fclose(fid);

toc

%% in-line function
function [] = saveSamplesToFile(fid,prnCode_noise,fileType,dataType)
switch fileType
    case 2
        for indSamp = 1:length(prnCode_noise)
            fwrite(fid,real(prnCode_noise(indSamp)),dataType);
            fwrite(fid,imag(prnCode_noise(indSamp)),dataType);
        end
    otherwise
        fprintf("File type '%s' not implemented...\n",string(fileType))
end
end