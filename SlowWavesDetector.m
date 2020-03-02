classdef SlowWavesDetector < handle
    
    properties
        
        samplingRate = 1000;
        
        lowLimitPhase = 0.16;
        highLimitPhase = 1.25;
       
        thresholdAmplitude = nan;
        usePredefinedThresholdAmplitude = false;
        predefinedThreshForAmp = 50; %from Lafon et al
        percentileForAmplitude = 40; %the percentile to leave (i.e. if it's set to 40 we set the threshold such that 40% of the waves will pass it)
        percentileForAmplitudeStaresina = 25;
        
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
            %data - in which we want to detect slow waves
            %Output -
            %slowWavesTimes - slow waves times corresponding to peak time
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
%                 troughInds(iCycle) = troughInds(iCycle)+cycStartInds(iCycle)-1;
                
%                 plot(cycStartInds(iCycle):cycEndInds(iCycle),dataFiltered(cycStartInds(iCycle):cycEndInds(iCycle)));
%                 hold all;
%                 plot(peakInds(iCycle),dataFiltered(peakInds(iCycle)),'*');
%                 hold all;
%                 currPhi = phiFP(cycStartInds(iCycle):cycEndInds(iCycle));
%                 zeroInd = zci(currPhi);
%                 zeroInd = zeroInd(currPhi(zeroInd)<1 & currPhi(zeroInd)>-1);
%                 zeroInd = cycStartInds(iCycle)-1+zeroInd;
%                 plot(zeroInd,dataFiltered(zeroInd),'*','color','g');
%                 hold all;
%                 plot(cycStartInds(iCycle):cycEndInds(iCycle),phiFP(cycStartInds(iCycle):cycEndInds(iCycle)));
%                 close all;
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
        

    end
end
