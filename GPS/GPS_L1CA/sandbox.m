clear all;
close all;

settings = initSettings();
prn = 1;
delay_true_chips = 1;
delay_true_samples = round(delay_true_chips/settings.codeFreqBasis*settings.samplingFreq); % in samples
doppler_true_Hz = 1500; % in Hz
signal_power_dBm = -128.5; % in dBm. Typical open sky power is -128.5 dBm

% change signal duration
settings.msToProcess = 1*1e3;

cohSum = settings.acqCohTime;

% samplesPerCode = round(settings.samplingFreq/(settings.codeFreqBasis/settings.codeLength));   
% samplesPerCoh = cohSum*samplesPerCode;
noise_psd = -203.9; % dBW/Hz. Typical value for RF front-end with integrated LNA
bw_Hz = 4e6;
noise_power_dBm = noise_psd + 10*log10(bw_Hz) + 30; % dBm
ts = 1 / settings.samplingFreq;
N_1ms = 1e-3*settings.samplingFreq; % number of samples in 1 ms
Npoints = settings.samplingFreq*settings.msToProcess/1e3;
phase_v = (0:(Npoints-1)) *2*pi*ts;

Nperiod = settings.msToProcess/1e3/(settings.codeLength/settings.codeFreqBasis);

% fprintf("Generating %2i ms of GPS L1 C/A signal for PRN #%02i at Fs = %4g Hz\n", settings.msToProcess, prn_v(indPrn), settings.samplingFreq);
% create 1ms of GPS L1 C/A signal
phi0 = 0;

% creates 1 ms of signal
prnCode_1ms = makeCaTable(prn,settings);

% introduce delay by circularly shifting the code period by delay
prnCode_delay = circshift(prnCode_1ms,delay_true_samples);

% introduce doppler
prnCode_doppler = prnCode_delay.*exp(1j*2*pi*(doppler_true_Hz*(0:N_1ms-1)/settings.samplingFreq)+phi0);
phi0 = phi0 + N_1ms*doppler_true_Hz/settings.samplingFreq;

% add noise
snr_dB = signal_power_dBm - noise_power_dBm;
std_noise = 10^(-snr_dB/20);
noise = std_noise*(randn(1,N_1ms)+1j*randn(1,N_1ms));
prnCode_noise = prnCode_doppler + noise;

% figure;
% subplot(2,1,1); hold all;
% plot(real(prnCode_noise));
% plot(real(prnCode_doppler));
% subplot(2,1,2); hold all;
% plot(imag(prnCode_noise));
% plot(imag(prnCode_doppler));

%% write signal to file
filename = ['simulatedSignal_tau=' num2str(delay_true_chips/settings.codeFreqBasis*1e6) 'µs_' ...
            'dop=' num2str(doppler_true_Hz) 'Hz_' ...
            'pow=' num2str(signal_power_dBm) 'dBm' ...
            '.bin'];
switch settings.fileType
    case 2
        fid = fopen(filename,'w');
        for indSamp = 1:N_1ms
            fwrite(fid,real(prnCode_noise(indSamp)),settings.dataType);
            fwrite(fid,imag(prnCode_noise(indSamp)),settings.dataType);
        end; clear indSamp
        fclose(fid);
    otherwise
        fprintf("File type '%s' not implemented...\n",string(settings.fileType))
end
