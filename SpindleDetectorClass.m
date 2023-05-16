classdef SpindleDetectorClass < handle
    
    %The class implements spindle detection and verification of channels
    %using the method described in 'Sleep Spindles in Humans: Insights from Intracranial EEG
    %and Unit Recordings' Andrillon, Nir, et al, j of neuroscience, 2011
    %(the default parameters are also based on the paper)
    
    properties
        samplingRate = 1000;

        spindleRangeMin = 9; %Hz % replace to 11Hz for high-freq spindles
        spindleRangeMax = 16; %Hz
        
        %verification step 1 parameters
        segmentsLengthStep1 = 10; %in seconds
        pvalThreshStep1 = 0.001;
        nanThresh = 0.01;
        minFitRange = 3;
        maxFitRange = 30;
        restrictFitPowerToNegative = true;
        
        %spindle detectione parameters
        detectionThresholdSD = 3; % number of SD
        detectionThresholdRejectionSD = 5; % number of SD
        detectionThresholdStartEndSD = 0.7; % number of SD
        eventMinDuration = 0.5; %sec
        eventMaxDuration = 2; % sec
        minDistBetweenEvents = 1; %sec. Any two events with less than this distance apart are merged
        rejectionRangeMin = 20;
        rejectionRangeMax = 30;
        minThreshForStageClassification = 0.7; %a spindle should have at least minThreshForStageClassification of its duration within a sleep stage to be classified as such
        
        %verification step 3 parameters
        timeWindowAroundSpindle = 1; % sec
        controlsPerSpindle = 2;
        minDistControlSpindle = 2; %sec
        maxDistControlSpindle = 5; %sec
        pvalThreshStep3 = 0.0001;
        
        %spindle detectione parameters - Staresina
        detectionThresholdPerc = 75; % percentile - top percentile set the threshold
        eventMinDurationStar = 0.5; %sec
        eventMaxDurationStar = 3; % sec
        RMSWindowDuration = 200; %ms

        %IIS removal constants
        windowAroundIIS = 500; %ms
        
        %scoring params
        scoringEpochDuration = 0.001; % How many seconds represented by one individual value in the scoring vector [scalar].
        sleepEpochs = [1]; % all the values in the scoring vector which represent sleep stages for which we want to perform the analysis (like NREM/REM/transitions) [1D vector].
        
        %constants for bandpass
        defaultFilterOrder = 2;
        nanWarning = 0.01;
        
        % constants for population analysis
        
        %STIMULATION removal constants
        windowAroundSTIM = 200; %ms
        
        winFromLastSpike = 1000; %ms
        shortTimeRangeAfterStim = 3;%seconds
        midTimeRangeAfterStim = 60; %seconds
        stimulusDuration = 50; %ms
        minSpikeRateToIncludeUnit = 1;
        
        avgRippleBeforeAfter = 1; %second
        freqoiForAvgSpec = [0:0.5:10];
        freqRangeForAvgSpec = [5:30];
        timeBeforeAfterEventRipSpec = 1; %second
        timeForBaselineRip = 1; %second
        minNCycles = 5;
        minWinSizeSpec = 100; %ms
        minNripples = 5;
        
        freqRangeForShowingSpindles = [5 30];
        
        %plotting constants
        blockSizePlot = 1;
        plotBeforeAfter = 500;
        windowSpec = 100;
        noverlapSpec = 80;
        nfftSpec = 1028;
        minThresholdSpec = -50;
        scalingFactorDeltaLog = 0; 
        sigmaImgaussfilt = 3;
        ylimPlot = 50;
        
        %simple plot constants
        beforeAfterSimple = 2;
        subplotSizeX = 5;
        subplotSizeY = 5;
        
    end
    
    methods
        
        function isVerified = verifyChannelStep1(obj,data,sleepScoring,IIStimes)
            
            % step 1, first verification
            % Input - data - data recorded from the channel we wish to verify
            % sleepScoring - a vector in which each element represents the
            % sleep stage for an epoch of length obj.sleepEpochs (e.g. 1
            % second). The values which represent the required sleep stages
            % are kept in obj.sleepEpochs
            % output - isVerified - a boolean, true if the channel is
            % verified as a channel with robust spindle activity
            
            if nargin < 3
                sleepScoring = [];
            end
            
            removeIIS = true;
            if nargin < 4 || isempty(IIStimes)
                removeIIS = false;
            end

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
                        data = data(1:length(isSleep));
                    end
                end
                %only leave segments of "real" sleep in data
                data(~isSleep) = nan;
            end
            
            %remove windowAroundIIS ms before and after every IIS as
            %provided as input parameter
            if removeIIS
                winAroundIIS = obj.windowAroundIIS*obj.samplingRate/1000;
                for iTime = 1:length(IIStimes)
                    pointsBefore = min(IIStimes(iTime),winAroundIIS);
                    pointsAfter = min(length(data)-IIStimes(iTime),winAroundIIS);
                    data(IIStimes(iTime)-pointsBefore+1:IIStimes(iTime)+pointsAfter) = nan;
                end
            end
            
            %divide the data to segments and calculate power spectrum and
            %fit for each segment
            
            segLength = obj.segmentsLengthStep1*obj.samplingRate;
            freq = 0:obj.samplingRate/segLength:obj.samplingRate/2;
            fitRangeInds = find(freq>=obj.minFitRange & freq<=obj.maxFitRange);
            freq = freq(fitRangeInds);
            
            nSegments = floor(length(data)/segLength);
            powerSpectrums = [];
            fitSpectrums = [];
            if obj.restrictFitPowerToNegative
                fitOpts = fitoptions('power2');
                fitOpts.Upper = [inf 0 inf];
            else
                fitOpts = fitoptions('power2');
            end
                
%             bs = [];
            
            for iSeg = 1:nSegments
                currSegment = data((iSeg-1)*segLength+1:iSeg*segLength);
                %if there are nans - we assume it means this is not a sleep
                %segment
                if sum(isnan(currSegment))/length(currSegment)>obj.nanThresh
                    continue;
                end
                psdx = obj.getPS(currSegment);
                psdx = psdx(fitRangeInds);
                psdx(isnan(psdx)) = 0; % zero out if we have a small percentage of NaN values
                %                 try
                f = fit(freq',psdx','power2',fitOpts);
                %                 catch
                %                     a=1;
                %                 end
                y = feval(f,freq);
                %                 bs(end+1) = f.b;
                
                powerSpectrums = [powerSpectrums;psdx];
                fitSpectrums = [fitSpectrums; y'];
            end
            
            %find the point of maximal difference inside the spindle range
            meanPS = nanmean(powerSpectrums,1);
            meanFS = nanmean(fitSpectrums, 1);
            
            spindleRangeInds = find(freq>=obj.spindleRangeMin & freq<=obj.spindleRangeMax);
            [~,maxInd] = max(meanPS(spindleRangeInds)-meanFS(spindleRangeInds));
            maxInd = spindleRangeInds(maxInd);
            
            %one tailed ttest for the difference between the power spectrum
            %and the fit
            [~,p] = ttest(powerSpectrums(:,maxInd),fitSpectrums(:,maxInd),'tail','right');
            
            %if the p-values is smaller than the threshold, the channel is
            %verified
            if p < obj.pvalThreshStep1
                isVerified = true;
            else
                isVerified = false;
            end
        end
        
        function [spindleTimes, spindleStats, startEndTimes] = detectSpindles(obj,data,sleepScoring, IIStimes, returnStats)
            
            % Step 2, spindle detection:
            % Input - data - data recorded from the channel we wish to verify
            % sleepScoring - a vector in which each element represents the
            % sleep stage for an epoch of length obj.sleepEpochs (e.g. 1
            % second). The values which represent the required sleep stages
            % are kept in obj.sleepEpochs
            % output - spindleTimes - a vector with spindle indices (the
            % estimated middle of the spindle event)
            
            if nargin < 3
                sleepScoring = [];
            end
            
            removeIIS = true;
            if nargin < 4 || isempty(IIStimes)
                removeIIS = false;
            end
            
            if nargin < 5
                returnStats = false;
            end
            
            %band pass in spindle range
            bandpassSignal  = obj.bandpass(data, obj.samplingRate, obj.spindleRangeMin, obj.spindleRangeMax);
            %band pass in rejection range (anything that will pass
            %threshold in that range will not be considered as detection)
            bandpassRejection  = obj.bandpass(data, obj.samplingRate, obj.rejectionRangeMin, obj.rejectionRangeMax);
            
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
                        bandpassSignal = bandpassSignal(1:length(isSleep));
                        bandpassRejection = bandpassRejection(1:length(isSleep));
                    end
                end
                %only leave segments of "real" sleep in data
                bandpassSignal(~isSleep) = nan;
                bandpassRejection(~isSleep) = nan;
            end
            
            %remove windowAroundIIS ms before and after every IIS as
            %provided as input parameter
            if removeIIS
                winAroundIIS = obj.windowAroundIIS*obj.samplingRate/1000;
                for iTime = 1:length(IIStimes)
                    pointsBefore = min(IIStimes(iTime),winAroundIIS);
                    pointsAfter = min(length(bandpassSignal)-IIStimes(iTime),winAroundIIS);
                    bandpassSignal(IIStimes(iTime)-pointsBefore+1:IIStimes(iTime)+pointsAfter) = nan;
                    bandpassRejection(IIStimes(iTime)-pointsBefore+1:IIStimes(iTime)+pointsAfter) = nan;
                end
            end
            
            %hilbert transform doesn't support NaN
            NaNInds = isnan(bandpassSignal);
            bandpassSignal(NaNInds) = 0;
            bandpassRejection(NaNInds) = 0;
            
            %calculate the envelops
            envelope = abs(hilbert(bandpassSignal));
            envelopeRejection = abs(hilbert(bandpassRejection));
            envelope(NaNInds) = NaN;
            envelopeRejection(NaNInds) = NaN;
            bandpassSignal(NaNInds) = nan;
            bandpassRejection(NaNInds) = nan;
            
            %find points which pass the peak threshold and the start/end
            %threshold, only among points which are sleep according to the
            %sleep scoring (if given)
            meanSleep = nanmean(envelope);
            stdSleep = nanstd(envelope);
            
            meanSleepRejection = nanmean(envelopeRejection);
            stdSleepRejection = nanstd(envelopeRejection);
            
            pointsPassedThresh = ((envelope-meanSleep)/stdSleep > obj.detectionThresholdSD) & ((envelopeRejection-meanSleepRejection)/stdSleepRejection < obj.detectionThresholdRejectionSD);
            pointsPassedThreshStartEnd = ((envelope-meanSleep)/stdSleep > obj.detectionThresholdStartEndSD);
            %calculate the diff for pointsPassedThreshStartEnd in order to
            %detect start and end of spindle events
            diffStartEnd = diff(pointsPassedThreshStartEnd);
            
            %events are defined by sequences which have a peak above
            %detectionThresholdSD and a duration within the limits of
            %required duration. The duration is considered to be from the
            %first point above pointsPassedThreshStartEnd until the last
            %point above it
            nextInd = find(pointsPassedThresh,1);
            spindleEventsMinLimit = [];
            spindleEventsMaxLimit = [];
            while ~isempty(nextInd)
                %last pass above the start threshold before current spindle
                startCurrEvent = find(diffStartEnd(1:nextInd-1)==1, 1, 'last')+1;
                %last pass above the end threshold after current spindle
                endCurrEvent = find(diffStartEnd(nextInd:end)==-1,1,'first')+nextInd-1;
                if ~isempty(startCurrEvent) && ~isempty(endCurrEvent)
                    %calcualte the duration
                    currDuration = (endCurrEvent-startCurrEvent)/obj.samplingRate;
                else
                    currDuration = nan;
                end
                %if the duration is within required limits, add to spindles
                %list
                if currDuration > obj.eventMinDuration && currDuration < obj.eventMaxDuration
                    spindleEventsMinLimit = [spindleEventsMinLimit; startCurrEvent];
                    spindleEventsMaxLimit = [spindleEventsMaxLimit; endCurrEvent];
                end
                nextInd = find(pointsPassedThresh(endCurrEvent+1:end),1,'first')+endCurrEvent;
            end
            
            %merge events which are too close apart, and set the spindle
            %timing to be the middle of the event
            eventDiffs = spindleEventsMinLimit(2:end)-spindleEventsMaxLimit(1:end-1);
            %find only events for which enough far apart
            differentEvents = [0 find(eventDiffs'>obj.minDistBetweenEvents*obj.samplingRate) length(spindleEventsMinLimit)];
            spindleTimes = zeros(1,length(differentEvents)-1);
            
            spindleStats = struct('startTime',[],'endTime',[],'duration',[],'peakEnergy',[],'peakTime',[],'sigmaPower',[],'freqSpindle',[],'currentStage',[]);
            %build final list of spindle times and statistics
            startEndTimes = zeros(length(differentEvents)-1,2);
            for iSpindle = 1:length(differentEvents)-1
                currStartTime = spindleEventsMinLimit(differentEvents(iSpindle)+1);
                currEndTime = spindleEventsMaxLimit(differentEvents(iSpindle+1));
%                 spindleTimes(iSpindle) = mean(currStartTime, currEndTime);
                
                currentSpindle = bandpassSignal(currStartTime:currEndTime) - mean(bandpassSignal(currStartTime:currEndTime));
                [maxSpindle,maxIndSpindle] = max(currentSpindle); 
                spindleTimes(iSpindle) = currStartTime+maxIndSpindle;
                startEndTimes(iSpindle,1) = currStartTime;
                startEndTimes(iSpindle,2) = currEndTime;
                
                %all the statistics are in ms
                if returnStats
                    msConversionConst = (1000/obj.samplingRate);
                    
                    %start time, end time and duration in ms
                    spindleStats(iSpindle).startTime = currStartTime*msConversionConst;
                    spindleStats(iSpindle).endTime = currEndTime*msConversionConst;
                    spindleStats(iSpindle).duration = spindleStats(iSpindle).endTime - spindleStats(iSpindle).startTime;
                    
                    %find peak time (ms) and energy (std relative to mean
                    %and std of sleep)
                    currentSpindle = bandpassSignal(currStartTime:currEndTime);
                    [maxSpindle,maxIndSpindle] = max(currentSpindle);
                    spindleStats(iSpindle).peakTime = spindleStats(iSpindle).startTime+maxIndSpindle*msConversionConst;
                    spindleStats(iSpindle).peakEnergy = (maxSpindle-meanSleep)/stdSleep;
                    
                    %find frequency with maximal energy within the spindle
                    %frequencies range, and the sum of the power spectrum
                    %within the range
                    [freq, pow] = obj.hereFFT(data(currStartTime:currEndTime));
                    rangeInds = find(freq>=obj.spindleRangeMin & freq<=obj.spindleRangeMax);
                    [~,maxPowInd] = max(pow(rangeInds));
                    maxPowInd = freq(rangeInds(maxPowInd));
                    spindleStats(iSpindle).freqSpindle = maxPowInd;
                    spindleStats(iSpindle).sigmaPower = sum(pow(rangeInds));
                    
                    %current sleep stage
                    if ~isempty(sleepScoring)
                        %this supports sleep stages in milliseconds
                        % TBD - change here if adding finer scoring (stage
                        % 1,2 SWS)
                        currSleepStages = sleepScoring(currStartTime:currEndTime);
                        uss = unique(sleepScoring);
                        histStages = hist(currSleepStages,uss);
                        histStages = histStages./sum(histStages);
                        currStage = find(histStages >= obj.minThreshForStageClassification);
                        if length(currStage)==1
                            spindleStats(iSpindle).currentStage = uss(currStage);
                        else
                            %mixed
                            spindleStats(iSpindle).currentStage = NaN;
                        end
                    end
                end
                
            end
            
        end
        
        function isVerified = verifyChannelStep3(obj,data,spindleTimes)
            % step 3, final verification
            % Input - data - data recorded from the channel we wish to verify
            % sleepScoring - a vector in which each element represents the
            % sleep stage for an epoch of length obj.sleepEpochs (e.g. 1
            % second). The values which represent the required sleep stages
            % are kept in obj.sleepEpochs
            % output - isVerified - a boolean, true if the channel is
            % verified as a channel with robust spindle activity
            
            nSpindles = length(spindleTimes);
            
            %find segments of spindles and control segments of non-spindle
            %near the spindle
            spindlesPS = [];
            controlPS = [];
            timeWindow = obj.timeWindowAroundSpindle*obj.samplingRate;
            %calculate the range in which control segments can selected
            rng = obj.maxDistControlSpindle-obj.minDistControlSpindle-obj.timeWindowAroundSpindle;
            
            for iSpindle = 1:nSpindles
                %data in and around the spindle
                currSpindleData = data(spindleTimes(iSpindle)-round(timeWindow/2):spindleTimes(iSpindle)+round(timeWindow/2)-1);
                spindlesPS = [spindlesPS; obj.getPS(currSpindleData)];
                
                
                for iControl = 1:obj.controlsPerSpindle
                    %rand sign
                    randSign = (-1)^randi(2);
                    %randomly choose an index for the control data
                    controlInd = round(spindleTimes(iSpindle)+(rand*rng+obj.minDistControlSpindle)*randSign*obj.samplingRate);
                    
                    if (controlInd+round(timeWindow/2)-1) > length(data)
                        currControlData = data(end-timeWindow+1:end);
                    elseif (controlInd-round(timeWindow/2)) < 0
                        currControlData = data(1:timeWindow);
                    else
                        currControlData = data(controlInd-round(timeWindow/2):controlInd+round(timeWindow/2)-1);
                    end
                    
                    controlPS = [controlPS; obj.getPS(currControlData)];
                end
            end
            
            freq = 0:obj.samplingRate/timeWindow:obj.samplingRate/2;
            
            %calcualte the mean of the power spectrum for the spindles and
            %the control and find the point of maximal difference within
            %spindle frequencies range
            meanSpindle = mean(spindlesPS,1);
            meanControl = mean(controlPS,1);
            
            spindleRangeInds = find(freq>=obj.spindleRangeMin & freq<=obj.spindleRangeMax);
            [~,maxInd] = max(meanSpindle(spindleRangeInds)-meanControl(spindleRangeInds));
            maxInd = spindleRangeInds(maxInd);
            
            %one tailed unpaired ttest for the difference between the power spectrum
            %of the spindles and the controls
            [~,p] = ttest2(spindlesPS(:,maxInd),controlPS(:,maxInd),'tail','right');
            
            %if the pvalue is below threshold teh channel is verified
            if p < obj.pvalThreshStep3
                isVerified = true;
            else
                isVerified = false;
            end
        end
        
        function plotSpindlesSimple(obj, data, spindlesTimes, folderToSave)
            
            %plots single simples and saves the figures as jpg in the
            %folder folderToSave if provided

            if nargin < 3 || isempty(folderToSave)
                toSave = false;                
            else
                toSave = true;
            end
            
            %convert window size to sampling points
            secondBefAfter = obj.beforeAfterSimple*obj.samplingRate;
            
            %filter data to required range - uncomment if filtered data
            %should also be plotted
%             filteredData = obj.bandpass(data, obj.samplingRate, obj.spindleRangeMin, obj.spindleRangeMax);
            
            nSpindles = length(spindlesTimes);
            nInPlot = obj.subplotSizeX*obj.subplotSizeY;
            nPlots = ceil(nSpindles/nInPlot);
            
            indSpindle = 1;
            figInd = 1;
            for iPlot = 1:nPlots-1
                f = figure;
                
                for iInPlot = 1:nInPlot
                    subplot(obj.subplotSizeY,obj.subplotSizeX,iInPlot);
                    minInd = max(spindlesTimes(indSpindle)-secondBefAfter,1);
                    maxInd = min(spindlesTimes(indSpindle)+secondBefAfter,length(data));
                    
%                     plot([-secondBefAfter:secondBefAfter]/obj.samplingRate,filteredData(minInd:maxInd),'color','r');
%                     hold all;
                    plot([-secondBefAfter:secondBefAfter]/obj.samplingRate,data(minInd:maxInd),'color','k');
%                     xlim([1 secondBefAfter*2]);
                    
                    title(['Spindle time = ', num2str(spindlesTimes(indSpindle)/obj.samplingRate/60),' mins']);
                    indSpindle = indSpindle+1;
                end
                if toSave
                        set(f, 'Position', get(0, 'Screensize'));
                        saveas(f, [folderToSave,'\all_spindles_' num2str(figInd),'.jpg']);
                        close(f);
                else
                    pause;
                end
                figInd = figInd+1;
            end
            
        end
        
        function plotSpindles(obj, data, spindleTimes, blockSize)
            %plots the spindles and their spectrogram. 
            %blockSize sets the amount of spindles to display per figure
            hasImGaussFilt = ~isempty(which('imgaussfilt'));
            if ~hasImGaussFilt
                warning('The matlab function imgaussfilt (available in the image processing toolbox) is not found on your path. Spectrograms will not be smoothed.')
            end
            if nargin < 4 || isempty(blockSize)
                blockSize = obj.blockSizePlot;
            end
            
            %obj.plotBeforeAfter is in ms - translate to number of data
            %points
            beforeAfterPoints = round(obj.plotBeforeAfter*obj.samplingRate/1000);
            
            nSpindles = length(spindleTimes);
            indBlock = 1;
            for iSpindle = 1:nSpindles
                
                %time frequncy
                ax(1) = subplot(2, blockSize, indBlock);
                %in case it's an edge (start/end of data) and the spike can be presented int he
                %middle of the block
                if spindleTimes(iSpindle)>beforeAfterPoints
                    minPoint = spindleTimes(iSpindle)-beforeAfterPoints;
                    spPoint = beforeAfterPoints;
                else
                    minPoint = 1;
                    spPoint = spindleTimes(iSpindle);
                end
                
                if spindleTimes(iSpindle)+beforeAfterPoints>length(data)
                    maxPoint = length(data);
                else
                    maxPoint = spindleTimes(iSpindle)+beforeAfterPoints;
                end
                
                currData = data(minPoint:maxPoint);
                currTS = 0:1/obj.samplingRate:length(currData)/obj.samplingRate;
                if length(currTS)>length(currData), currTS = currTS(1:length(currData));end
                plot(currTS,currData);
                hold all;
                plot(currTS(spPoint),min(currData)*2,'marker','*','color','r');                
                hold off;
                title(['Spindle #',num2str(iSpindle)]);
                
                %spectogram
                ax(2) = subplot(2, blockSize, indBlock+blockSize);
                [S,F,T,P] = spectrogram(currData,obj.windowSpec,obj.noverlapSpec,obj.nfftSpec,obj.samplingRate,'yaxis','MinThreshold',obj.minThresholdSpec);
                P = P/max(max(P));
                P1 = (10*log10(abs(P+obj.scalingFactorDeltaLog)))';
                P1 = [P1(:,1) P1 P1(:,end)];
                T = [0 T T(end)+median(diff(T))];
                if hasImGaussFilt
                    imagesc(T,F,imgaussfilt(P1',obj.sigmaImgaussfilt),[obj.minThresholdSpec,0]);
                else
                    imagesc(T,F,P1',[obj.minThresholdSpec,0]);
                end
                axis xy;
%                 imagesc(T,F,P1',[obj.minThresholdSpec,0]);axis xy;
                set(ax(2),'ylim',[0, obj.ylimPlot]);
                
                linkaxes(ax,'x')
                
                indBlock = indBlock+1;
                if indBlock > blockSize || iSpindle == nSpindles
                    indBlock = 1;
                    pause;
                end
            end
        end
        
        function [spindleTimes, startEndTimes] = detectSpindlesStaresina(obj,data,sleepScoring, IIStimes)
            
            %Spindles detection based on Staresina et al 2015
            
            %convert sizes according to the sampling rate
            RMSWindowDuration = obj.RMSWindowDuration*obj.samplingRate/1000;
            minDurationAboveThresh = obj.eventMinDurationStar*obj.samplingRate;
            maxDurationAboveThresh = obj.eventMaxDurationStar*obj.samplingRate;
            minDistBetweenSpindles = obj.minDistBetweenEvents*obj.samplingRate;
            
            if nargin < 3
                sleepScoring = [];
            end
            
            removeIIS = true;
            if nargin < 4 || isempty(IIStimes)
                removeIIS = false;
            end
            
            %filter data to spindle range
            filteredData = obj.bandpass(data, obj.samplingRate, obj.spindleRangeMin, obj.spindleRangeMax);
            %filter data to rejection range - note this is not part of the
            %method described in Staresina but an addition
            bandpassRejection  = obj.bandpass(data, obj.samplingRate, obj.rejectionRangeMin, obj.rejectionRangeMax);
            
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
                if length(isSleep)>length(filteredData)
                    isSleep = isSleep(1:length(filteredData));
                else if length(isSleep)<length(filteredData)
                        filteredData = filteredData(1:length(isSleep));
                        data = data(1:length(isSleep));
                        bandpassRejection = bandpassRejection(1:length(isSleep));
                    end
                end
                %only leave segments of "real" sleep in data
                filteredData(~isSleep) = nan;
                data(~isSleep) = nan;
                bandpassRejection(~isSleep) = nan;
            end
            
            %remove windowAroundIIS ms before and after every IIS as
            %provided as input parameter
            if removeIIS
                winAroundIIS = obj.windowAroundIIS*obj.samplingRate/1000;
                for iTime = 1:length(IIStimes)
                    pointsBefore = min(IIStimes(iTime),winAroundIIS);
                    pointsAfter = min(length(filteredData)-IIStimes(iTime),winAroundIIS);
                    filteredData(IIStimes(iTime)-pointsBefore+1:IIStimes(iTime)+pointsAfter) = nan;
                    data(IIStimes(iTime)-pointsBefore+1:IIStimes(iTime)+pointsAfter) = nan;
                    bandpassRejection(IIStimes(iTime)-pointsBefore+1:IIStimes(iTime)+pointsAfter) = nan;
                end
            end
            
            NaNInds = isnan(bandpassRejection);
            bandpassRejection(NaNInds) = 0;
            
            %calculate the envelope of the rejection range
            envelopeRejection = abs(hilbert(bandpassRejection));
            envelopeRejection(NaNInds) = NaN;
%             bandpassRejection(NaNInds) = nan;
            
            meanSleepRejection = nanmean(envelopeRejection);
            stdSleepRejection = nanstd(envelopeRejection);
            
            %find points where the zscored envelope of the rejection range passed the rejection threshold
            passedThreshRejection = (envelopeRejection-meanSleepRejection)/stdSleepRejection >= obj.detectionThresholdRejectionSD;
                        
            %calculate the root mean squared signal of the filtered (spindle range) data for windows of length 
            %RMSWindowDuration
            rmsSignal = zeros(1,length(filteredData)-RMSWindowDuration+1);
            for iPoint = 1:length(filteredData)-RMSWindowDuration+1
                rmsSignal(iPoint) = rms(filteredData(iPoint:iPoint+RMSWindowDuration-1));
            end
            
            %calculate the threshold as the rippleThreshPercentile
            %percentile of the rms signal
            spindleThresh = prctile(rmsSignal,obj.detectionThresholdPerc);
            
            %find windows that pass the thresh
            didPassThresh = rmsSignal>=spindleThresh;
            
            %find segments that pass the threshold for a duration longer than the threshold:
            %eventMinDurationStar milliseconds and shorter than eventMaxDurationStar
            spindleSegs = [];
            ind = 1;
            indSpindle = 1;
            while ind <= length(didPassThresh)-minDurationAboveThresh+1
                if all(didPassThresh(ind:ind+minDurationAboveThresh-1))
                    endSeg = ind+find(didPassThresh(ind:end)==0,1) - 2;
                    segLength = endSeg-ind;
                    if segLength<=maxDurationAboveThresh
                        spindleSegs(indSpindle,1) = ind;
                        spindleSegs(indSpindle,2) = endSeg + RMSWindowDuration - 1;
                        indSpindle = indSpindle+1;
                    end
                    ind = endSeg+2;
                else
                    ind = ind+1;
                end
            end
            
            %merge spindles who are close (distance between end of segment
            %and beginning of next segment less than
            %minDistBetweenSpindles)
            spindlesSegsMerged = [];
            spindlesDiffsSmall = (spindleSegs(2:end,1)-spindleSegs(1:end-1,2))<minDistBetweenSpindles;
            indOld = 1;
            indNew = 0;
            
            while indOld < size(spindleSegs,1)
                indNew = indNew+1;
                if spindlesDiffsSmall(indOld)==0
                    spindlesSegsMerged(indNew,:) = spindleSegs(indOld,:);
                    indOld = indOld+1;
                else
                    nextMerge = find(spindlesDiffsSmall(indOld+1:end)==0,1)+indOld;
                    if isempty(nextMerge)
                        nextMerge = size(spindleSegs,1);
                    end
                    spindlesSegsMerged(indNew,:) = [spindleSegs(indOld,1) spindleSegs(nextMerge,2)];
                    indOld = nextMerge+1;
                end
            end
            if sum(spindlesDiffsSmall)==0 || nextMerge<size(spindleSegs,1)
                spindlesSegsMerged(end+1,:) = spindleSegs(end,:);
            end
            
            startEndTimes = spindlesSegsMerged;
            nSpindles = size(startEndTimes,1);
            
            %remove spindles who pass the rejection threshold - where any
            %of the points within the spindle segment passes the threshold
            indsRemove = [];
            for iSpindle = 1:nSpindles
                if any(passedThreshRejection(startEndTimes(iSpindle,1):startEndTimes(iSpindle,2)))
                    indsRemove = [indsRemove iSpindle];
                end
            end
                
            startEndTimes(indsRemove,:) = [];
            nSpindles = size(startEndTimes,1);
            
            spindleTimes = zeros(1,nSpindles);
            
            %return stats is not supported in this method - if required
            %should be added

            for iSpindle = 1:nSpindles
                currStartTime = startEndTimes(iSpindle,1);
                currEndTime = startEndTimes(iSpindle,2);
                
                currentSpindle = filteredData(currStartTime:currEndTime);
                [maxSpindle,maxIndSpindle] = max(currentSpindle);
                %spindle times are the indices of the maximal value
                spindleTimes(iSpindle) = currStartTime+maxIndSpindle;
                
            end
            
        end
        
        %% help functions
        
        function psdx = getPS(obj,segment)
            %an help method to calcualte the power spectrum of a segment
            
            segLength = length(segment);
            xdft = fft(segment);
            xdft = xdft(1:segLength/2+1);
            psdx = (1/(obj.samplingRate*segLength)) * abs(xdft).^2;
            psdx(2:end-1) = 2*psdx(2:end-1);
            psdx = 10*log10(psdx);
            
        end
        
        function BP = bandpass(obj, timecourse, SamplingRate, low_cut, high_cut, filterOrder)
            
            %bandpass code - from Maya
            
            if (nargin < 6)
                filterOrder = obj.defaultFilterOrder;
            end
            
            % Maya GS - handle NAN values
            indices = find(isnan(timecourse));
            if length(indices) > obj.nanWarning*length(timecourse)
                warning('many NaN values in filtered signal')
            end
            timecourse(indices) = 0;
            %
            
            [b, a] = butter(filterOrder, [(low_cut/SamplingRate)*2 (high_cut/SamplingRate)*2]);
            BP = filtfilt(b, a, timecourse );
            BP(indices) = NaN;
        end
        
        function [f, pow] = hereFFT (obj, signal)
            % Calculate the fft of the signal, with the given sampling rate, make
            % normalization(AUC=1). Return requencis and respective powers.
            
            % Matlab sorce code of FFT
            Y = fft(signal);
            power_spec = Y.* conj(Y) / length(signal);
            
            % keep only half of the power spectrum array (second half is irrelevant)
            amp = power_spec(1:ceil(length(power_spec)/2)) ;
            
            % Define the frequencies relevant for the left powers, and cut for same
            % number of values (for each frequency - a power value)
            f = obj.samplingRate*(1:length(amp))/(length(amp)*2);
            pow = amp(1:length(f));
            
            %----- End of Regular fft -----
            
            pow = pow / sum (pow);   % normalize AUC
            
        end
        

        
    end
    
end
