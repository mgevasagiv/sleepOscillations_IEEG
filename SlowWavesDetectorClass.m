classdef SlowWavesDetectorClass < handle
    
    properties
        
        samplingRate = 1000;
        
        lowLimitPhase = 0.16;
        highLimitPhase = 1.25;
       
        thresholdAmplitude = nan;
        usePredefinedThresholdAmplitude = false;
        predefinedThreshForAmp = 50; %from Lafon et al
        percentileForAmplitude = 40; %the percentile to leave (i.e. if it's set to 40 we set the threshold such that 40% of the waves will pass it)
        percentileForAmplitudeStaresina = 25;
        
        desiredPercentOfSlowWaves = 20; % This is a lower amplitude percentile for SHAM analyis
        
        
        minLengthThresh = 0.8; %seconds
        maxLengthThresh = 2; %seconds
        
        %filtering constants
        defaultFilterOrder = 1;
        nanWarning = 0.01;
        
        %IIS removal constants
        windowAroundIIS = 500; %ms
        
        %sleep scoring parameters
        scoringEpochDuration = 0.001; % How many seconds represented by one individual value in the scoring vector [scalar].
        sleepEpochs = [1]; % all the values in the scoring vector which represent sleep stages for which we want to perform the analysis (like NREM/REM/transitions) [1D vector].
        
        %staresina detection constants
        isPosToNeg = false;
        nanThresh = 0.3;
        
        %Maingret detection constants
        zscoreThreshPeak = 2;
        zscoreThreshPeakWithEnd = 1;
        zscoreThreshEnd = -1.5;
    end
    
    methods
        
        function [slowWavesTimes, slowWavesPeaks] = findSlowWaves(obj, data, sleepScoring, IIStimes)
            
            %The method detects slow waves in the data based on their
            %maximal amplitude and duration.
            %
            %Input - 
            %data - in which we want to detect slow waves
            %Output - 
            %slowWavesTimes - an array where the number of rows is the
            %number of detected slow waves, in column 1 appears the
            %first index of the slow wave, in column 2 the last index
            %
            %The method implemented here is different than the Staresina method:
            %A. Finds candidates differently - by finding change of phase
            %from pi to -pi
            %B. Amplitude is absolute and not relative to trough
            %C. Threshold can be set apriorly and not just relative to
            %input data.
            
            if nargin < 3
                sleepScoring = [];
            end
            
            removeIIS = true;
            if nargin < 4 || isempty(IIStimes)
                removeIIS = false;
            end
            
            %filter the data at the slow waves frequency band
            dataFilteredFPhase = obj.bandpass(data, obj.lowLimitPhase, obj.highLimitPhase);
            
            % if sleepScoring is nonempty: leave only the
            % segments in which there was sleep at the desired stage for
            % the analysis
            if ~isempty(sleepScoring)
                segLength = obj.scoringEpochDuration*obj.samplingRate;
                isSleep = zeros(1,length(sleepScoring)*segLength);
                for iEpoch = 1:length(sleepScoring)
                    if ismember(sleepScoring(iEpoch),obj.sleepEpochs)
                        isSleep((iEpoch-1)*segLength+1:iEpoch*segLength) = ones(1,segLength);
                    end
                end
                %match the length of data with the length of sleepScoring -
                %might get rid of some data points at the end if required, assuming it's
                %negligible
                if length(isSleep)>length(data)
                    isSleep = isSleep(1:length(data));
                else if length(isSleep)<length(data)
                        dataFilteredFPhase = dataFilteredFPhase(1:length(isSleep));
                    end
                end
                %only leave segments of "real" sleep in data
                dataFilteredFPhase(~isSleep) = nan;
            end
            
            %remove windowAroundIIS ms before and after every IIS as
            %provided as input parameter
            if removeIIS
                winAroundIIS = obj.windowAroundIIS*obj.samplingRate/1000;
                for iTime = 1:length(IIStimes)
                    pointsBefore = min(IIStimes(iTime),winAroundIIS);
                    pointsAfter = min(length(dataFilteredFPhase)-IIStimes(iTime),winAroundIIS);
                    dataFilteredFPhase(IIStimes(iTime)-pointsBefore+1:IIStimes(iTime)+pointsAfter) = nan;
                end
            end
            
            %remove nan inds before the transform
            nanIndsP = isnan(dataFilteredFPhase);
            dataFilteredFPhase(nanIndsP) = 0;
            
            %calculate the time series of the phases of dataFilteredFPhase
            phiFP = angle(hilbert(dataFilteredFPhase));
            %Make original NaN indices NaN again
            phiFP(nanIndsP) = nan;
            
            %find points where new cycles start (a jump from phase pi to
            %phase -pi)
            diffPhiFP = diff(phiFP);
            newCyclePoints = find(diffPhiFP<-6)+1;
            nCycles = length(newCyclePoints)-1;
            
            maxAmpSlowPerCycle = nan(1, nCycles);
            lsCycles = nan(1, nCycles);
            
            %calculate the peak amplitude of the LF data per cycle and the
            %lengths
            for iCycle = 1:nCycles
                currCycleSlow = dataFilteredFPhase(newCyclePoints(iCycle):newCyclePoints(iCycle+1)-1);
                
                maxAmpSlowPerCycle(iCycle) = nanmax(currCycleSlow);
                lsCycles(iCycle) = length(currCycleSlow)/obj.samplingRate;
            end
            
            %if there is no pre-defined threshold for slow waves amplitude, use the current data to
            %find a threshold - at the percentileForSlowAmpOutlier
            %percentile
            if isnan(obj.thresholdAmplitude) && ~obj.usePredefinedThresholdAmplitude
                cyclesPassedAmpThresh = maxAmpSlowPerCycle>=prctile(maxAmpSlowPerCycle,100-obj.percentileForAmplitude);
            else
                %if there is an absolute predefined threshold, use
                %the predefined threshold
                if obj.usePredefinedThresholdAmplitude
                    currThresh = obj.predefinedThreshForAmp;
                else
                    %if the threshold is not predefined but was calculated
                    %previously (on a larger data set for example), use
                    %the previously claculated threshold
                    currThresh = obj.thresholdAmplitude;
                end
                cyclesPassedAmpThresh = maxAmpSlowPerCycle >= currThresh;
            end
            
            cyclesPassedLengthThresh = lsCycles >= obj.minLengthThresh & lsCycles <= obj.maxLengthThresh;
            cyclesPassedThresh = find(cyclesPassedAmpThresh & cyclesPassedLengthThresh);
            
            nCycles = length(cyclesPassedThresh);
            slowWavesTimes = zeros(nCycles,2);
            slowWavesPeaks = zeros(1,nCycles);
            for iCycle = 1:nCycles
                slowWavesTimes(iCycle,:) = [newCyclePoints(cyclesPassedThresh(iCycle)),newCyclePoints(cyclesPassedThresh(iCycle)+1)-1];
                currCycleSlow = dataFilteredFPhase(slowWavesTimes(iCycle,1):slowWavesTimes(iCycle,2));
                [~,peakInd] = max(currCycleSlow);
                slowWavesPeaks(iCycle) = slowWavesTimes(iCycle,1)+peakInd-1;
            end
            
        end
        
        function slowWavesTimes = findSlowWavesStaresina(obj, data, sleepScoring, IIStimes)
            %The method detects slow waves in the data based on their
            %maximal amplitude and duration, as described in Staresina et
            %al 2015
            %Input - 
            % data - in which we want to detect slow waves
            % Output -
            % slowWavesTimes - slow waves times corresponding to peak time
            % Slow wave candidates are between zero crossings of the
            % filtered data, only candidates with the highest 25% amplitude
            % (peak-trough) and within duration limits are kept
            
            minLengthThresh = obj.minLengthThresh*obj.samplingRate;
            maxLengthThresh = obj.maxLengthThresh*obj.samplingRate;
            
            if nargin < 3
                sleepScoring = [];
            end
            
            removeIIS = true;
            if nargin < 4 || isempty(IIStimes)
                removeIIS = false;
            end
            
            %filter the data at the slow waves frequency band
            dataFiltered = obj.bandpass(data, obj.lowLimitPhase, obj.highLimitPhase);
            
            % if sleepScoring is nonempty: leave only the
            % segments in which there was sleep at the desired stage for
            % the analysis
            if ~isempty(sleepScoring)
                segLength = obj.scoringEpochDuration*obj.samplingRate;
                isSleep = zeros(1,length(sleepScoring)*segLength);
                for iEpoch = 1:length(sleepScoring)
                    if ismember(sleepScoring(iEpoch),obj.sleepEpochs)
                        isSleep((iEpoch-1)*segLength+1:iEpoch*segLength) = ones(1,segLength);
                    end
                end
                %match the length of data with the length of sleepScoring -
                %might get rid of some data points at the end if required, assuming it's
                %negligible
                if length(isSleep)>length(data)
                    isSleep = isSleep(1:length(data));
                else if length(isSleep)<length(data)
                        dataFiltered = dataFiltered(1:length(isSleep));
                    end
                end
                %only leave segments of "real" sleep in data
                dataFiltered(~isSleep) = nan;
            end
            
            %remove windowAroundIIS ms before and after every IIS as
            %provided as input parameter
            if removeIIS
                winAroundIIS = obj.windowAroundIIS*obj.samplingRate/1000;
                for iTime = 1:length(IIStimes)
                    pointsBefore = min(IIStimes(iTime),winAroundIIS);
                    pointsAfter = min(length(dataFiltered)-IIStimes(iTime),winAroundIIS);
                    dataFiltered(IIStimes(iTime)-pointsBefore+1:IIStimes(iTime)+pointsAfter) = nan;
                end
            end
            
            %find zero crossings of the filtered data
            zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
            zeroCross = zci(dataFiltered);
            
            % can't have a zero-crossing in the first sample
            zeroCross(zeroCross == 1) = [];
            
            %candidates can be from a positive to negative change until
            %negative to positive to change or the other way around
            if obj.isPosToNeg
                cycleInds = zeroCross(dataFiltered(zeroCross-1)>0);
            else
                cycleInds = zeroCross(dataFiltered(zeroCross-1)<0);
            end
            
            %find durations
            cycDurations = [cycleInds(2:end)-cycleInds(1:end-1) ; 0];
            %leave only cycles who pass the duration conditions
            cycStartInds = cycleInds(cycDurations >= minLengthThresh & cycDurations <= maxLengthThresh);
            cycEndInds = cycleInds(find(cycDurations >= minLengthThresh & cycDurations <= maxLengthThresh)+1);
            
            %find cycles' amplitudes
            nCycles = length(cycStartInds);
            ampCycles = zeros(1,nCycles);
            peakInds = zeros(1,nCycles);
            troughInds = zeros(1,nCycles);
            tooManyNaN = zeros(1,nCycles); 
            for iCycle = 1:nCycles
                currData = dataFiltered(cycStartInds(iCycle):cycEndInds(iCycle));
                %remove candidates with too many NaNs
                if sum(isnan(currData)) > obj.nanThresh
                    tooManyNaN(iCycle) = 1;
                    continue;
                end
                [peakVal, peakInds(iCycle)] = max(currData);
                [troughVal, troughInds(iCycle)] = min(currData);
                
                %the amplitude is defined as peak-trough
                ampCycles(iCycle) = peakVal-troughVal;
                
                peakInds(iCycle) = peakInds(iCycle)+cycStartInds(iCycle)-1;

            end
            
            %slow waves are the candidates in the high
            %obj.percentileForAmplitude percentile of amplitude
            cycsPassedThresh = ampCycles>prctile(ampCycles(~tooManyNaN),100-obj.percentileForAmplitudeStaresina);
            %slow wave times are the peak times
            slowWavesTimes = peakInds(cycsPassedThresh);
%             troughInds = troughInds(cycsPassedThresh);
        end
        
        function slowWavesTimes = findSlowWavesMaingret(obj, data, sleepScoring, IIStimes)
            %The method detects slow waves in the data based on their
            %maximal amplitude and duration, as described in Maingret et
            %al 2016
            %Input - 
            %data - in which we want to detect slow waves
            %Output -
            %slowWavesTimes - slow waves times corresponding to peak time
            %The method for finding candidates: they
            %calculate the derivative of the filtered data and find zero
            %crossing in the derivative that are a change from negative to
            %positive (i.e. troughs in the filtered data). Amplitudes are
            %defined as the max value of the zscored filtered data. The
            %amplitude threshold (zscored) is set a-priori.
            
            minLengthThresh = obj.minLengthThresh*obj.samplingRate;
            maxLengthThresh = obj.maxLengthThresh*obj.samplingRate;
            
            if nargin < 3
                sleepScoring = [];
            end
            
            removeIIS = true;
            if nargin < 4 || isempty(IIStimes)
                removeIIS = false;
            end
            
            %filter the data at the slow waves frequency band
            dataFiltered = obj.bandpass(data, obj.lowLimitPhase, obj.highLimitPhase);
            
            % if sleepScoring is nonempty: leave only the
            % segments in which there was sleep at the desired stage for
            % the analysis
            if ~isempty(sleepScoring)
                segLength = obj.scoringEpochDuration*obj.samplingRate;
                isSleep = zeros(1,length(sleepScoring)*segLength);
                for iEpoch = 1:length(sleepScoring)
                    if ismember(sleepScoring(iEpoch),obj.sleepEpochs)
                        isSleep((iEpoch-1)*segLength+1:iEpoch*segLength) = ones(1,segLength);
                    end
                end
                %match the length of data with the length of sleepScoring -
                %might get rid of some data points at the end if required, assuming it's
                %negligible
                if length(isSleep)>length(data)
                    isSleep = isSleep(1:length(data));
                else if length(isSleep)<length(data)
                        dataFiltered = dataFiltered(1:length(isSleep));
                    end
                end
                %only leave segments of "real" sleep in data
                dataFiltered(~isSleep) = nan;
            end
            
            %remove windowAroundIIS ms before and after every IIS as
            %provided as input parameter
            if removeIIS
                winAroundIIS = obj.windowAroundIIS*obj.samplingRate/1000;
                for iTime = 1:length(IIStimes)
                    pointsBefore = min(IIStimes(iTime),winAroundIIS);
                    pointsAfter = min(length(dataFiltered)-IIStimes(iTime),winAroundIIS);
                    dataFiltered(IIStimes(iTime)-pointsBefore+1:IIStimes(iTime)+pointsAfter) = nan;
                end
            end
            
            %find zerocrossings
            zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
            %find the derivative of the data using splines
            fdata = spline([1:length(dataFiltered)],dataFiltered);
            fdataDer = ppval(fnder(fdata,1),[1:length(dataFiltered)]);
            fdataDer(isnan(dataFiltered)) = nan;
            %find zeros crossings of the derivative
            zeroCross = zci(fdataDer);     
            
            %find change of sign from neg to pos (troughs)
            cycleInds = zeroCross(fdataDer(zeroCross-1)<0);
            
            %find durations
            cycDurations = [cycleInds(2:end)-cycleInds(1:end-1) ; 0];
            %leave only cycles who pass the duration conditions
            cycStartInds = cycleInds(cycDurations >= minLengthThresh & cycDurations <= maxLengthThresh);
            cycEndInds = cycleInds(find(cycDurations >= minLengthThresh & cycDurations <= maxLengthThresh)+1);
            
            %z score the data
            zscoreData = (dataFiltered-nanmean(dataFiltered))/nanstd(dataFiltered);
            
            %calculate amplitudes of the candidates - absolute max value of the zscore
            nCycles = length(cycStartInds);
            peakCycles = zeros(1,nCycles);
            endCycles = zeros(1,nCycles);
            peakInds = zeros(1,nCycles);
            
            for iCycle = 1:nCycles
                [peakCycles(iCycle), peakInds(iCycle)] = max(zscoreData(cycStartInds(iCycle):cycEndInds(iCycle)));
                endCycles(iCycle) = zscoreData(cycEndInds(iCycle));                
                peakInds(iCycle) = peakInds(iCycle)+cycStartInds(iCycle)-1;
            end
            
            %candidates pass the threshold either if their max is above
            %zscoreThreshPeak or if their max is above
            %zscoreThreshPeakWithEnd and their end value is above zscoreThreshEnd
            cycsPassedThresh = peakCycles>obj.zscoreThreshPeak | (peakCycles>obj.zscoreThreshPeakWithEnd & endCycles>obj.zscoreThreshEnd);
            slowWavesTimes = peakInds(cycsPassedThresh);
            
        end
        
        function setSlowWaveThresh(obj, data, sleepScoring, IIStimes)
            %The method calculated the threshold for the amplitude for slow
            %wave detection based on the input data. It sets the threshold
            %to be such that obj.percentileForAmplitude percents of the
            %slow waves (of the waves that are at the correct duration
            %range) will pass this threshold. It sets the property
            %thresholdAmplitude to the calculated value.
            
            
            if nargin < 3
                sleepScoring = [];
            end
            
            removeIIS = true;
            if nargin < 4 || isempty(IIStimes)
                removeIIS = false;
            end
            
            %filter the data at the slow waves frequency band
            dataFilteredFPhase = obj.bandpass(data, obj.lowLimitPhase, obj.highLimitPhase);
            
            % if sleepScoring is nonempty: leave only the
            % segments in which there was sleep at the desired stage for
            % the analysis
            if ~isempty(sleepScoring)
                segLength = obj.scoringEpochDuration*obj.samplingRate;
                isSleep = zeros(1,length(sleepScoring)*segLength);
                for iEpoch = 1:length(sleepScoring)
                    if ismember(sleepScoring(iEpoch),obj.sleepEpochs)
                        isSleep((iEpoch-1)*segLength+1:iEpoch*segLength) = ones(1,segLength);
                    end
                end
                %match the length of data with the length of sleepScoring -
                %might get rid of some data points at the end if required, assuming it's
                %negligible
                if length(isSleep)>length(data)
                    isSleep = isSleep(1:length(data));
                else if length(isSleep)<length(data)
                        dataFilteredFPhase = dataFilteredFPhase(1:length(isSleep));
                    end
                end
                %only leave segments of "real" sleep in data
                dataFilteredFPhase(~isSleep) = nan;
            end
            
            %remove windowAroundIIS ms before and after every IIS as
            %provided as input parameter
            if removeIIS
                winAroundIIS = obj.windowAroundIIS*obj.samplingRate/1000;
                for iTime = 1:length(IIStimes)
                    pointsBefore = min(IIStimes(iTime),winAroundIIS);
                    pointsAfter = min(length(dataFilteredFPhase)-IIStimes(iTime),winAroundIIS);
                    dataFilteredFPhase(IIStimes(iTime)-pointsBefore+1:IIStimes(iTime)+pointsAfter) = nan;
                end
            end
            
            %remove nan inds before the transform
            nanIndsP = isnan(data);
            
            dataFilteredFPhase(nanIndsP) = 0;
            
            %calculate the time series of the phases of dataFilteredFPhase
            phiFP = angle(hilbert(dataFilteredFPhase));
            %Make original NaN indices NaN again
            phiFP(nanIndsP) = nan;
            
            %find points where new cycle start (a jump from phase pi to
            %phase -pi)
            diffPhiFP = diff(phiFP);
            newCyclePoints = find(diffPhiFP<-6)+1;
            nCycles = length(newCyclePoints)-1;
            
            maxAmpSlowPerCycle = nan(1, nCycles);
            lsCycles = nan(1, nCycles);
            
            %calculate the mean squared amplitude of the HF data per cycle
            for iCycle = 1:nCycles
                currCycleSlow = dataFilteredFPhase(newCyclePoints(iCycle):newCyclePoints(iCycle+1)-1);
                maxAmpSlowPerCycle(iCycle) = nanmax(currCycleSlow);
                lsCycles(iCycle) = length(currCycleSlow)/obj.samplingRate;
            end
           
            %calculate the threshold only on the cycles that are at the
            %fitting slow waves length
            cyclesPassedLengthThresh = lsCycles >= obj.minLengthThresh & lsCycles <= obj.maxLengthThresh;
            maxAmpSlowPerCycle = maxAmpSlowPerCycle(cyclesPassedLengthThresh);
            
            %set the threshold to be at the 100-obj.percentileForSlowAmpOutlier percentile
            obj.thresholdAmplitude = prctile(maxAmpSlowPerCycle,100-obj.percentileForAmplitude);
            
        end
        
        % Modified from detect_SW_positive_dealWithSAW.m    
        function [upStateTimes, downStateTimes] = detect_SW_positive(obj, data, peakTimes)
            
            n_points_block_size_filtfilt = 10^6; % I just chose a block size that is small enough NOT to cause "Out Of Memory" problems
            if n_points_block_size_filtfilt > length(data)
                n_points_block_size_filtfilt = length(data);
            end
            numOfSegments = floor( length(data) / n_points_block_size_filtfilt );
            for ii_block = 1 : numOfSegments,
                
                idx_block = (1:n_points_block_size_filtfilt)+(ii_block-1)*n_points_block_size_filtfilt ;
                eegSegments_BP(ii_block,:) = data(idx_block);
            end % The "problematic discontinuities" in data_filt_ripples - on the borders of 2 blocks - are SMALL in size
            if length( data( idx_block(end)+1 : end ) ) > 1
                eegSegments_BP(ii_block+1,:) = NaN(1,n_points_block_size_filtfilt);
                lastVec = data(idx_block(end)+1 : end);
                eegSegments_BP(ii_block+1,1:length(lastVec))  = lastVec;
            end
            
            numOfSegments = size(eegSegments_BP,1);
            
            % So far work on one channel only call results 'waves1'
            uvValue = -7777;
            
            %% Resample to 100Hz for detection
            fs = 100;
            ss_factor = 1; % default no submsampling is necessary
            if (obj.samplingRate > 100)
                ss_factor = obj.samplingRate / fs;
            end
            eegSegments_BP_ss = zeros(numOfSegments, size(eegSegments_BP,2)/ss_factor);
            if (obj.samplingRate > 100)
                for i=1:numOfSegments
                    eegSegments_BP_ss(i, :) = resample(eegSegments_BP(i, :), fs, obj.samplingRate);
                end
            else
                eegSegments_BP_ss = eegSegments_BP;
            end
            
            
            %%%%%% Brady's code from here on
            wavesAllSegments = [];
            cntBadsegments = 0;
            for currentSegment = 1:numOfSegments
                % EEG is the data vector sampled at 100Hz
                if sum(isnan(eegSegments_BP_ss(currentSegment, :))) > 1000
                    % if it is a consecutive piece, it won't tamper with SW detection
                    % we don't want a scenario where there's a fixed rate of missing
                    % timestamps
                    % Count how scattered the missing-data are
                    if sum(diff(find(isnan(eegSegments_BP_ss(currentSegment, :)))) ~= 1) > (0.05*n_points_block_size_filtfilt/100) % 0.05% of size
                        disp(sprintf('SW detection - seg %d/%d - too many missing samples, abort',currentSegment,numOfSegments))
                        cntBadsegments = cntBadsegments + 1;
                        continue
                    end
                    if sum(isnan(eegSegments_BP_ss(currentSegment, :))) > 0.5*n_points_block_size_filtfilt/ss_factor
                        disp(sprintf('SW detection - seg %d/%d - too many missing samples, abort',currentSegment,numOfSegments))
                        cntBadsegments = cntBadsegments + 1;
                        continue
                    end
                end
                pos_index=zeros(length(eegSegments_BP_ss(currentSegment, :)),1);
                pos_index(find(eegSegments_BP_ss(currentSegment, :)>0))=1; %index of all positive points for EEG
                difference=diff(pos_index); poscross=find(difference==1) ; negcross=find(difference==-1); %find neg ZX and pos ZX
                EEGder=obj.meanfilt(diff(eegSegments_BP_ss(currentSegment, :)),5); %meanfilt is a function that uses a 5 sample moving window to smooth derivative
                pos_index=zeros(length(EEGder),1);
                pos_index(find(EEGder>0.1))=1; %index of all positive points above minimum threshold
                difference=diff(pos_index);
                peaks=find(difference==-1)+1; troughs=find(difference==1)+1; %find pos ZX and neg ZX of the derivative (the peaks & troughs)
                peaks(eegSegments_BP_ss(currentSegment, peaks)<0 | isnan(eegSegments_BP_ss(currentSegment, peaks)))=[]; % rejects peaks below zero and troughs above zero
                troughs(eegSegments_BP_ss(currentSegment, troughs)>0 | isnan(eegSegments_BP_ss(currentSegment, troughs)))=[]; % rejects peaks below zero and troughs above zero
                
                if negcross(1)<poscross(1);start=1;else start=2;end %makes negcross and poscross same size to start
                if start==2;poscross(1)=[];end
                
                lastpk=NaN; %way to look at Peak to Peak parameters if needed
                
                waves = zeros(length(negcross)-start, 28);
                uvValueLine = ones(1, 28) * uvValue;
                
                for wndx=start:length(negcross)-1
                    
                    wavest=negcross(wndx);   %only used for neg/pos peaks
                    wavend=negcross(wndx+1); %only used for neg/pos peaks
                    mxdn=abs(nanmin(obj.meanfilt(diff(eegSegments_BP_ss(currentSegment, wavest:poscross(wndx))),5)))*fs;     % matrix (27) determines instantaneous positive 1st segement slope on smoothed signal, (name not representative)
                    mxup=nanmax(obj.meanfilt(diff(eegSegments_BP_ss(currentSegment, wavest:poscross(wndx))),5))*fs;     % matrix (28) determines maximal negative slope for 2nd segement (name not representative)
                    negpeaks=troughs(troughs>wavest&troughs<wavend);
                    
                    % In case a peak is not detected for this wave (happens rarely)
                    if (size(negpeaks,1) == 0)
                        waves(wndx, :) = uvValueLine;
                        continue;
                    end
                    
                    pospeaks=peaks(peaks>wavest&peaks<=wavend);
                    if isempty(pospeaks);pospeaks=wavend; end %if negpeaks is empty set negpeak to pos ZX
                    period=wavend-wavest; %matrix(11) /fs
                    poszx=poscross(wndx); %matrix(10)
                    b=nanmin(eegSegments_BP_ss(currentSegment, negpeaks)); % matrix (12) most pos peak /abs for matrix
                    if b>0;b=b(1);end;
                    bx=negpeaks(eegSegments_BP_ss(currentSegment, negpeaks)==b); %matrix (13) max pos peak location in entire night
                    c=nanmax(eegSegments_BP_ss(currentSegment, pospeaks)); % matrix (14) most neg peak
                    if c>0;c=c(1);end;
                    cx=pospeaks(eegSegments_BP_ss(currentSegment, pospeaks)==c); %matrix (15) max neg peak location in entire night
                    maxb2c=c-b; % %matrix (16) max peak to peak amp
                    nump=length(negpeaks); %matrix(24) now number of positive peaks
                    n1=abs(eegSegments_BP_ss(currentSegment, negpeaks(1))); %matrix(17) 1st pos peak amp
                    n1x=negpeaks(1); %matrix(18) 1st pos peak location
                    nEnd=abs(eegSegments_BP_ss(currentSegment, negpeaks(end))); %matrix(19) last pos peak amp
                    nEndx=negpeaks(end);%matrix(20) last pos peak location
                    p1=eegSegments_BP_ss(currentSegment, pospeaks(1)); %matrix(21) 1st neg peak amp
                    p1x=pospeaks(1); %matrix(22) 1st pos peak location
                    meanAmp=abs(mean(eegSegments_BP_ss(currentSegment, negpeaks))); %matrix(23)
                    nperiod=poszx-wavest; %matrix (25)neghalfwave period
                    mdpt=wavest+ceil(nperiod/2); %matrix(9)
                    %         epoch=ceil(bx/(fs*epochsize)); %matrix(1)
                    epoch = uvValue;
                    %         smepoch=ceil(bx/(fs*withinsize)); %matrix(2)
                    smepoch = uvValue;
                    p2p=(cx-lastpk)/fs; %matrix(26) 1st peak to last peak period
                    lastpk=cx;
                    cycle = uvValue; qcycle = uvValue; session = uvValue;
                    % UV: Indicate the 10sec segment # in the first column of the waves data structure
                    waves(wndx, :) = [currentSegment smepoch uvValue cycle qcycle session wavest wavend mdpt poszx period/fs abs(b) bx c cx maxb2c n1 n1x nEnd nEndx p1 p1x meanAmp nump nperiod/fs p2p mxdn mxup];
                    
                end %end wndx loop
                wavesAllSegments = [wavesAllSegments; waves];
            end % of loop through segments
            if cntBadsegments > 0.5*numOfSegments
                disp('bad channel,aborting SW detection');
                upStateTimes = NaN;
                downStateTimes = NaN;
                finalSlowWaves = [];
                return
            end
            slowWaves = wavesAllSegments((wavesAllSegments(:, 25)<1) & (wavesAllSegments(:, 25)>0.25), :); % choose slow waves based on their period
            
            % IDX 29 is for staging - currently not used
            slowWaves(:, 29) = 999; % place a random number for non-staged naps

            clear finalSlowWaves; 
            beforeMax = 500;
            beforeMin = 1;
            corruptedUpStates = [];     goodUpStates      = [];
            goodRows = []; badRows = [];
            % Use only up state times to detect waves preceded by SAW complexes
            for currentSegment = 1:numOfSegments
                %% First get up state times within this segment in 1000Hz so
                %% we can compare them to the already stored SAW file:
                idx_for_this_segment = find(slowWaves(:, 1) == currentSegment);
                if (isempty(idx_for_this_segment))
                    continue;
                end
                % get timings of slow waves in the decimated/subsampled vectors
                upStateTimes_100Hz   = slowWaves(idx_for_this_segment, 13);
                % get timinings of slow waves in the original EEG vector
                subsampledTimeline = 1:ss_factor:n_points_block_size_filtfilt;
                upStateTimesForThisSegment   = subsampledTimeline(upStateTimes_100Hz);
                
                rowsWithSpikeTimesForThisSegment = find(peakTimes(:) < n_points_block_size_filtfilt*currentSegment);
                spikeTimesForThisSegment = peakTimes(rowsWithSpikeTimesForThisSegment) - n_points_block_size_filtfilt*(currentSegment - 1);
                
                %% Separating up states into good and bad
                positionOfBadUpStatesInRowForThisSegment  = [];
                positionOfGoodUpStatesInRowForThisSegment  = [];
                for currentUpState = 1:length(upStateTimesForThisSegment)
                    currentUpTimeWithinSegment = upStateTimesForThisSegment(currentUpState);
                    adjacentSpikes = find( (spikeTimesForThisSegment > (currentUpTimeWithinSegment-beforeMax)) & (spikeTimesForThisSegment < (currentUpTimeWithinSegment-beforeMin)) );
                    if (~isempty(adjacentSpikes))
                        positionOfBadUpStatesInRowForThisSegment = [positionOfBadUpStatesInRowForThisSegment; currentUpState];
                    else
                        positionOfGoodUpStatesInRowForThisSegment = [positionOfGoodUpStatesInRowForThisSegment; currentUpState];
                    end
                end
                goodRows = [goodRows; idx_for_this_segment(positionOfGoodUpStatesInRowForThisSegment) ];
                badRows  = [badRows;  idx_for_this_segment(positionOfBadUpStatesInRowForThisSegment) ];
            end
            slowWavesClean      = slowWaves(goodRows, :);
            slowWavesCorrupted  = slowWaves(badRows, :);
            
            %% Now work separately on good or bad waves -
            %% First - GOOD WAVES
            %% select a subset of slow waves with highest amplitudes
            numOfSlowWaves = size(slowWavesClean,1);
            allAmplitudes = slowWavesClean(:, 12);
            [sortedAmplitudes,sortedIndices] = sort(allAmplitudes, 'descend');
            slowWavesSortedByAmplitude = slowWavesClean(sortedIndices, :);
            cutoffNumber = round((obj.desiredPercentOfSlowWaves/100) * numOfSlowWaves);
            finalSlowWaves = slowWavesSortedByAmplitude(1:cutoffNumber, :);
            
            %% Define downStateTimes and upStateTimes - each is a two-column matrix
            %% where first column is the segment and second is time in ms within that
            %% segment in original (1000Hz) sampling frequency
            downStateTimes = zeros(size(finalSlowWaves,1), 2);
            upStateTimes   = zeros(size(finalSlowWaves,1), 2);
            for currentSegment = 1:numOfSegments
                idx_for_this_segment = find(finalSlowWaves(:, 1) == currentSegment);
                % get timings of slow waves in the decimated/subsampled vectors
                downStateTimes_100Hz = finalSlowWaves(idx_for_this_segment, 15);
                upStateTimes_100Hz   = finalSlowWaves(idx_for_this_segment, 13);
                
                % get timinings of slow waves in the original EEG vector
                subsampledTimeline = 1:ss_factor:n_points_block_size_filtfilt;
                
                downStateTimes(idx_for_this_segment, 1) = currentSegment;
                downStateTimes(idx_for_this_segment, 2) = subsampledTimeline(downStateTimes_100Hz);
                
                upStateTimes(idx_for_this_segment, 1) = currentSegment;
                upStateTimes(idx_for_this_segment, 2) = subsampledTimeline(upStateTimes_100Hz);
                
            end
            
            
        end % func
        
        
        function BP = bandpass(obj, timecourse, lowLimit, highLimit, filterOrder)
            
            %bandpass code - from Maya
            
            if (nargin < 5)
                filterOrder = obj.defaultFilterOrder;
            end
            
            % Maya GS - handle NAN values
            indices = find(isnan(timecourse));
            if length(indices) > obj.nanWarning*length(timecourse)
                warning('many NaN values in filtered signal')
            end
            timecourse(indices) = 0;
            %
            
            [b, a] = butter(filterOrder, [(lowLimit/obj.samplingRate)*2 (highLimit/obj.samplingRate)*2]);
            BP = filtfilt(b, a, timecourse );
            BP(indices) = NaN;
        end
        
        % same local meanfilt as used in SW-detector
        function [filtdata] = meanfilt(obj, datatofilt,pts)
            
            if length(datatofilt)>=pts
                filtdata=[];
                ptsaway=floor(pts/2); isEven = (ptsaway*2 == pts);
                filtdata([1:pts])=datatofilt([1:pts]);
                filtdata([length(datatofilt)-(pts-1):length(datatofilt)])=datatofilt([length(datatofilt)-(pts-1):length(datatofilt)]);
                for wndw=pts-ptsaway:length(datatofilt)-(pts-ptsaway)
                    filtdata(wndw)=nanmean(datatofilt([wndw-(ptsaway)+ isEven:wndw+(ptsaway)] ));
                end
            else filtdata=datatofilt;
            end
            
        end
        

    end
end
