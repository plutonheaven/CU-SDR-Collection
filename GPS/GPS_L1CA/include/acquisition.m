function acqResults = acquisition(longSignal, settings)
%Function performs cold start acquisition on the collected "data". It
%searches for GPS signals of all satellites, which are listed in field
%"acqSatelliteList" in the settings structure. Function saves code phase
%and frequency of the detected signals in the "acqResults" structure.
%
%acqResults = acquisition(longSignal, settings)
%
%   Inputs:
%       longSignal    - 11 ms of raw signal from the front-end 
%       settings      - Receiver settings. Provides information about
%                       sampling and intermediate frequencies and other
%                       parameters including the list of the satellites to
%                       be acquired.
%   Outputs:
%       acqResults    - Function saves code phases and frequencies of the 
%                       detected signals in the "acqResults" structure. The
%                       field "carrFreq" is set to 0 if the signal is not
%                       detected for the given PRN number. 

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
% (C) Updated by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
% Based on the original work by Darius Plausinaitis,Peter Rinder, 
% Nicolaj Bertelsen and Dennis M. Akos
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

%CVS record:
%$Id: acquisition.m,v 1.1.2.12 2006/08/14 12:08:03 dpl Exp $

%% Acquisition initialization =============================================

% Coherent sum number
cohSum = settings.acqCohTime;
% Non coherent sum number
nonCohSum = settings.acqNonCohTime;

% Find number of samples per spreading code
samplesPerCode = round(settings.samplingFreq / ...
                        (settings.codeFreqBasis / settings.codeLength));                   
% Find number of samples per coherent integration
samplesPerCoh = cohSum*samplesPerCode;

% Create matrix of data for the non coherent integration. 2 times the
% number of non coherent sum is needed to avoid bit transition
Signal = zeros(2*nonCohSum,samplesPerCoh);
for ii = 1:nonCohSum
    Signal(2*(ii-1) + 1,:) = longSignal(2*(ii-1)*samplesPerCoh + 1:2*(ii-1)*samplesPerCoh + samplesPerCoh);
    Signal(2*ii,:) = longSignal((2*(ii-1) + 1)*samplesPerCoh + 1:(2*(ii-1) + 1)*samplesPerCoh + samplesPerCoh);   
end

signal0DC = longSignal - mean(longSignal);  

% Find sampling period
ts = 1 / settings.samplingFreq;

% Find phase points of the local carrier wave 
phasePoints = (0 : (samplesPerCoh-1)) * 2 * pi * ts;

%--- List of frequency bins to search for ---------------------------------
freqBinList =  settings.freqBinList;

%--- Initialize acqResults ------------------------------------------------
% Carrier frequencies of detected signals
acqResults.carrFreq     = zeros(1, 32);
% PRN code phases of detected signals
acqResults.codePhase    = zeros(1, 32);
% Correlation peak ratios of the detected signals
acqResults.peakMetric   = zeros(1, 32);
% Number of frequency bin searched to obtain above threshold peak
acqResults.freqBin      = zeros(1, 32);
% acquisition matrix
acqResults.acqMat       = zeros(32,length(freqBinList),samplesPerCoh);

fprintf('(');

% Perform search for all listed PRN numbers ...
for PRN = settings.acqSatelliteList

    %--- Search results of all frequency bins and code shifts (for one satellite)
    results = zeros(length(freqBinList), samplesPerCoh); 
    %--- Generate CA code and make it X ms long ---------------------------
    caCode = repmat(makeCaTable(PRN,settings),1,cohSum);
    %--- Perform DFT of PRN code ------------------------------------------
    caCodeFreqDom = conj(fft(caCode));

    %--- Make the correlation for whole frequency band (for all freq. bins)
    for freqBinIndex = 1:length(freqBinList)
        
        freqBin = settings.IF + freqBinList(freqBinIndex);
        %--- Generate local sine and cosine -------------------------------
        sigCarr = exp(-1i * freqBin * phasePoints);
        
        %--- Initialization of non coherent results
        nonCohResultsA = zeros(1,samplesPerCoh);
        nonCohResultsB = zeros(1,samplesPerCoh);
        
        %For all the non coherent sum
        for ii = 1:nonCohSum
            %--- "Remove carrier" from the signal -----------------------------
            I1      = real(sigCarr .* Signal(2*(ii-1) + 1,:));
            Q1      = imag(sigCarr .* Signal(2*(ii-1) + 1,:));
            I2      = real(sigCarr .* Signal(2*ii,:));
            Q2      = imag(sigCarr .* Signal(2*ii,:));
            
            %--- Convert the baseband signal to frequency domain --------------
            IQfreqDomA = fft(I1 + 1i*Q1);
            IQfreqDomB = fft(I2 + 1i*Q2);
            
            %--- Multiplication in the frequency domain (correlation in time
            %domain)
            convCodeIQA = IQfreqDomA .* caCodeFreqDom;
            convCodeIQB = IQfreqDomB .* caCodeFreqDom;
            
            %--- Perform inverse DFT and store correlation results ------------
            cohResultsA = abs(ifft(convCodeIQA));
            cohResultsB = abs(ifft(convCodeIQB));
            
            %--- Add this value for the non coherent sum
            nonCohResultsA = nonCohResultsA + cohResultsA;
            nonCohResultsB = nonCohResultsB + cohResultsB;                      
        end
                    
        %--- Check which non coh had the greater power and save that, will
        %"blend" 1st and 2nd msec but will correct data bit issues
        if (max(nonCohResultsA) > max(nonCohResultsB))
            results(freqBinIndex, :) = nonCohResultsA;
        else
            results(freqBinIndex, :) = nonCohResultsB;
        end
        
    end % frqBinIndex = 1:numberOfFrqBins

%% Look for correlation peaks in the results ==============================
    % Find the highest peak and compare it to the second highest peak
    % The second peak is chosen not closer than 1 chip to the highest peak
    
    %--- Find the correlation peak and the carrier frequency --------------
    [~, frequencyBinIndex] = max(max(results, [], 2));

    %--- Find code phase of the same correlation peak ---------------------
    [peakSize, codePhase] = max(max(results));

    %--- Find 1 chip wide C/A code phase exclude range around the peak ----
    samplesPerCodeChip   = round(settings.samplingFreq / settings.codeFreqBasis);
    excludeRangeIndex1 = codePhase - samplesPerCodeChip;
    excludeRangeIndex2 = codePhase + samplesPerCodeChip;

    %--- Correct PRN code phase exclude range if the range includes array
    %boundaries
    if excludeRangeIndex1 < 1
        codePhaseRange = excludeRangeIndex2 : ...
                         (samplesPerCode + excludeRangeIndex1);
                         
    elseif excludeRangeIndex2 >= samplesPerCode
        codePhaseRange = (excludeRangeIndex2 - samplesPerCode) : ...
                         excludeRangeIndex1;
    else
        codePhaseRange = [1:excludeRangeIndex1, ...
                          excludeRangeIndex2 : samplesPerCode];
    end

    %--- Find the second highest correlation peak in the same freq. bin ---
    secondPeakSize = max(results(frequencyBinIndex, codePhaseRange));

    %--- Store result -----------------------------------------------------
    acqResults.peakMetric(PRN) = peakSize/secondPeakSize;
    acqResults.acqMat(PRN,:,:) = results;
    
    % If the result is above threshold, then there is a signal ...
    if (peakSize/secondPeakSize) > settings.acqThreshold

%% Fine resolution frequency search =======================================
        
        %--- Indicate PRN number of the detected signal -------------------
        fprintf('%02d ', PRN);
        %--- Save acquisition results
        acqResults.codePhase(PRN) = mod(codePhase,samplesPerCode);
        acqResults.carrFreq(PRN) = settings.IF +  freqBinList(frequencyBinIndex);
        acqResults.freqBin(PRN) = frequencyBinIndex;        
    else
        %--- No signal with this PRN --------------------------------------
        fprintf('. ');
    end   % if (peakSize/secondPeakSize) > settings.acqThreshold
    
end    % for PRN = satelliteList

%=== Acquisition is over ==================================================
fprintf(')\n');
