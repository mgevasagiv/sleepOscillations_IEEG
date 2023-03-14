classdef SpindleDetector < handle
    
    %The class implements spindle detection and verification of channels
    %using the method described in 'Sleep Spindles in Humans: Insights from Intracranial EEG
    %and Unit Recordings' Andrillon, Nir, et al, j of neuroscience, 2011
    %(the default parameters are also based on the paper)
    
    properties
        samplingRate = 1000;

        spindleRangeMin = 9; %Hz
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
                subplot(2, blockSize, indBlock);
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
                plot(currData);
                hold all;
                plot(spPoint,min(currData)*2,'marker','*','color','r');                
                hold off;
                title(['Spindle #',num2str(iSpindle)]);
                
                %spectogram
                subplot(2, blockSize, indBlock+blockSize);
                [S,F,T,P] = spectrogram(currData,obj.windowSpec,obj.noverlapSpec,obj.nfftSpec,obj.samplingRate,'yaxis','MinThreshold',obj.minThresholdSpec);
%                 [S,F,T,P] = spectrogram(currData,obj.windowSpec,obj.noverlapSpec,obj.nfftSpec,obj.samplingRate);
                P = P/max(max(P));
                P1 = (10*log10(abs(P+obj.scalingFactorDeltaLog)))';
                P1 = [P1(:,1) P1 P1(:,end)];
                T = [0 T T(end)+1];
                imagesc(T,F,imgaussfilt(P1',obj.sigmaImgaussfilt),[obj.minThresholdSpec,0]);
                axis xy;
%                 imagesc(T,F,P1',[obj.minThresholdSpec,0]);axis xy;
                set(gca,'ylim',[0, obj.ylimPlot]);
                
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
        
        
        % Revision - based on runRippleData from RippleDetector
        function results = runSpindleData(obj, runData, fileNameResults, whatToRun)
            
            % The method produces information about the ripples in a channel that includes:
            % A. Average ripple � before stimulation and during stimulations (short effect).
            % B. Spectrum of average ripple � before stimulation and during stimulations (short effect).
            % C. Average of TFR around ripples � before stimulation and during stimulations (short effect).
            % D. Average of TFR around spindles � before stimatulation.
            % E. Polar histogram of synchronization index between spindles and ripples range.
            %
            % The input runData is a struct in the length of number of patients (for which the analysis is required).
            % In addition it receives the input parameter fileNameResults which includes the file name into which the results
            % will be saved (optional).
            % Each element (=patient) in runData should include the fields:
            % patientName
            % channelsToRunOn � list of channel indices for which to perform the analysis.
            % DataFolder � The folder in which the raw data files are saved (the method assumes the prefix for the files is
            % CSC, can be changed by the property dataFilePrefix).
            % macroMontageFileName - the file name (including path) of the macromontage.
            % RipplesFileNames - name (including path) of the ripple mat files in which the ripple times for the macro
            % channels are saved (the method assumes the name of the file is RipplesFileNames <#channel index>
            % SpindlesFileNames - name (including path) of the spindle mat files in which the spindle times for the macro
            % channels are saved (the method assumes the name of the file is SpindlesFileNames <#channel index>
            % SpikesFileNames - name (including path) of the spikes mat files in which the spikes times for the macro channels
            % are saved (the method assumes the name of the file is SpikesFileNames <#channel index>). If not provided spikes
            % will not be removed from the data for the analysis.
            % sleepScoringFileName � file name (including path) of the sleep scoring mat file. If not provided all the data will be used.
            %
            % The output struct results includes all the results of the analysis, which can then be plotted using
            % plotResultsRipplesData. The output struct is a struct with the length of the number of patients (=the length
            % of runData), where each element includes:
            % patientName
            % resultsPerChan � a struct in the length of the number of channels required for the analysis per the patient.
            % Each element in resultsPerChan includes the fields:
            % channelNum
            % area
            % nRipplesBefore, nRipplesStim - number of ripples before stimulations and after stimulations (short effect)
            % respectively
            % avgBefore, avgStim � average ripples before stimulations and after stimulations (short effect) respectively
            % stdBefore, stdStim � std of ripples before stimulations and after stimulations (short effect) respectively
            % specBefore, specStim � spectrum of ripples average before stimulations and after stimulations (short effect)
            % respectively
            % meanTFRRipBefore, meanTFRRipStim � mean of ripple triggered TFR before stimulations and after stimulations
            % (short effect) respectively
            % SIanglesSpRip � Synchronization Indices of spindles-ripples before stimaulations (an array with the length as
            % number of spindles)
            % R � results of the r-test for the polar histogram (of SIanglesSpRip)
            % V � results of the v-test for the polar histogram (of SIanglesSpRip)
            % meanSpecs � mean spindle-triggered TFR before stimulations.
            % meanEpochs � mean spindle for the spindles for which meanSpecs were calculated.
            % nEpochs � number of spindles in the spindle-ripple analyses (all before stimatulions).
            
            
            if nargin < 3
                fileNameResults = '';
            end
            
            removeIIS = 1; 
            removeSTIM_artifacts = 1;
            
            shortTimeRangeAfterStim = obj.shortTimeRangeAfterStim*obj.samplingRate;
            midTimeRangeAfterStim = obj.midTimeRangeAfterStim*obj.samplingRate;
            stimulusDuration = obj.stimulusDuration*obj.samplingRate/1000;
            avgRippleBeforeAfter = obj.avgRippleBeforeAfter*obj.samplingRate; %ms
            timeWin = min((1./obj.freqRangeForAvgSpec)*obj.minNCycles,ones(size(obj.freqRangeForAvgSpec))*obj.minWinSizeSpec);
            
            nPatients = length(runData);
            
            %go over all required patients
            for iPatient = 1:nPatients
                disp(['Patient ',runData(iPatient).patientName,' ',num2str(iPatient),'/',num2str(nPatients)]);
                results(iPatient).patientName = runData(iPatient).patientName;
                
                %load exp data for stimulation timings
                expData = load(runData(iPatient).ExpDataFileName);
                expData = expData.EXP_DATA;
                stimTimes = expData.stimTiming.validatedTTL_NLX;
                firstStim = stimTimes(1);
                
                %load sleep scoring
                if isfield(runData(iPatient), 'sleepScoringFileName') && ~isempty(runData(iPatient).sleepScoringFileName)
                    try
                        sleepScoring = load(runData(iPatient).sleepScoringFileName);
                        sleepScoring = sleepScoring.sleep_score_vec;
                    catch
                        disp([runData(iPatient).sleepScoringFileName '.mat doesn''t exist']);
                        sleepScoring = [];
                    end
                end
                
                %load macro montage
                macroMontage = load(runData(iPatient).macroMontageFileName);
                macroMontage = macroMontage.MacroMontage;
                
                %go over required channels per patient
                resultsPerChan = [];
                nChans = length(runData(iPatient).channelsToRunOn);
                for iChan = 1:nChans
                    currChan = runData(iPatient).channelsToRunOn(iChan);
                    resultsPerChan(iChan).channelNum = currChan;
                    
                    fprintf(['channel ',num2str(currChan),',',num2str(iChan),'/',num2str(nChans)])
                    
                    currArea = macroMontage(resultsPerChan(iChan).channelNum).Area;
                    resultsPerChan(iChan).area = currArea;
                    
                    %load the data
                    try
                        data = [runData(iPatient).DataFolder '\CSC' num2str(currChan) '.mat'];
                        data = load(data);
                        data = data.data;
                    catch
                        disp([runData(iPatient).DataFolder '\CSC' num2str(currChan) '.mat doesn''t exist']);
                        continue;
                    end
                    
                    %load IIS times
                    if isfield(runData(iPatient), 'SpikesFileNames') && ~isempty(runData(iPatient).SpikesFileNames)
                        try
                            IIStimes = [runData(iPatient).SpikesFileNames num2str(currChan) '.mat'];
                            IIStimes = load(IIStimes);
                            IIStimes = IIStimes.peakTimes;
                        catch
                            disp([runData(iPatient).SpikesFileNames num2str(currChan) '.mat doesn''t exist']);
                            IIStimes = [];
                        end
                    else
                        IIStimes = [];
                    end
                    % Remove IIS from data before using it for ripple
                    % average and TFR:
                    % remove windowAroundSTIM ms before and after every stimulation as
                    % provided as input parameter
                    if removeSTIM_artifacts
                        winAroundIIS = obj.windowAroundSTIM*obj.samplingRate/1000;
                        pointsBefore = winAroundIIS;
                        pointsAfter = winAroundIIS;
                        data(stimTimes-pointsBefore+1:stimTimes+pointsAfter) = nan;
                    end

                    %remove windowAroundIIS ms before and after every IIS as
                    %provided as input parameter
                    if removeIIS
                        winAroundIIS = obj.windowAroundIIS*obj.samplingRate/1000;
                        pointsBefore = winAroundIIS;
                        pointsAfter = winAroundIIS;
                        
                        IIStimes(IIStimes < winAroundIIS) = [];
                        IIStimes(IIStimes > length(data) - winAroundIIS) = [];
                        
                        data(IIStimes-pointsBefore+1:IIStimes+pointsAfter) = nan;
                    end
                    

                    %load spindles
                    if isfield(runData(iPatient), 'SpindlesFileNames') && ~isempty(runData(iPatient).SpindlesFileNames)
                        try
                            if whatToRun.HighFreqSpindles == 0
                                spindlesTimes = [runData(iPatient).SpindlesFileNames num2str(currChan) '.mat'];
                            else
                                spindlesTimes = [runData(iPatient).HighFreqSpindlesFileNames num2str(currChan) '.mat'];
                            end
                            
                            spindlesTimes = load(spindlesTimes);
                            spindlesTimes = spindlesTimes.spindlesTimes;
                        catch
                            disp([runData(iPatient).SpindlesFileNames num2str(currChan) '.mat doesn''t exist']);
                            spindlesTimes = [];
                        end
                    end
                    
                    ripplesTimes = spindlesTimes;
                    %indices of ripples before stimulations
                    ripplesBeforeInds = ripplesTimes<firstStim-obj.windowAroundSTIM;
                    %get inds of ripples that are short and mid effect and not too
                    %close to the stimulus
                    dataDuration = floor(stimTimes(end))+midTimeRangeAfterStim;
                    stimInds = zeros(1,dataDuration);
                    %short effect
                    for iStim = 1:length(stimTimes)
                        stimInds(stimTimes(iStim)+stimulusDuration:stimTimes(iStim)+shortTimeRangeAfterStim) = 1;
                    end
                    %add also mid effect
                    stimDiffs = diff(stimTimes);
                    stimIndsWithMidPauseAfter = [find(stimDiffs >= midTimeRangeAfterStim) length(stimTimes)];
                    stimIndsWithMidPauseAfterTimes = stimTimes(stimIndsWithMidPauseAfter);
                    for iStim = 1:length(stimIndsWithMidPauseAfter)
                        stimInds(stimIndsWithMidPauseAfterTimes(iStim)+shortTimeRangeAfterStim:stimIndsWithMidPauseAfterTimes(iStim)+midTimeRangeAfterStim) = 1;
                    end
                    % MGS - edge effect
                    if (length(stimInds) > dataDuration); stimInds(dataDuration+1:end) = []; end;
                    ripplesIndsLog = zeros(1,dataDuration);
                    ripplesIndsLog(ripplesTimes(ripplesTimes<=dataDuration)) = 1;
                    %inds of short and mid effect ripples
                    ripplesDuringStimInds = ismember(ripplesTimes,find(ripplesIndsLog & stimInds));
                    
                    
                    %calculate average ripples and TFRs
                    ripplesTimesBefore = ripplesTimes(ripplesBeforeInds);
                    nRipplesBefore = length(ripplesTimesBefore);
                    
                    ripplesTimesStim = ripplesTimes(ripplesDuringStimInds);
                    nRipplesStim = length(ripplesTimesStim);
                    
                    %calculate average of ripples
                    ripplesBefore = zeros(nRipplesBefore,length([-avgRippleBeforeAfter:avgRippleBeforeAfter]));
                    ripplesStim = zeros(nRipplesStim,length([-avgRippleBeforeAfter:avgRippleBeforeAfter]));
                    

                    if nRipplesBefore > obj.minNripples
                        %before stimulations
                        rmvInd = ripplesTimesBefore < avgRippleBeforeAfter;
                        ripplesTimesBefore(rmvInd) = [];
                        rmvInd = ripplesTimesBefore > (length(data)-avgRippleBeforeAfter);
                        ripplesTimesBefore(rmvInd) = [];

                        for iRipple = 1:length(ripplesTimesBefore)
                            ripplesBefore(iRipple,:) = data(ripplesTimesBefore(iRipple)-avgRippleBeforeAfter:ripplesTimesBefore(iRipple)+avgRippleBeforeAfter);
                        end
   
                        avgBefore = nanmean(ripplesBefore);
                        stdBefore = nanstd(ripplesBefore);
                        %spectrum of average
                        specBefore = ft_specest_mtmfft(avgBefore,[1:length(avgBefore)]/obj.samplingRate,'freqoi',obj.freqoiForAvgSpec,'taper','hanning');
                        specBefore = abs(squeeze(specBefore(1,1,:,:)));
                    else
                        avgBefore = nan(1,size(ripplesBefore,2));
                        stdBefore = nan(1,size(ripplesBefore,2));
                        specBefore = nan(length(obj.freqoiForAvgSpec),1);
                    end
                    
                    %after stimulations
                    for iRipple = 1:nRipplesStim
                        ripplesStim(iRipple,:) = data(ripplesTimesStim(iRipple)-avgRippleBeforeAfter:ripplesTimesStim(iRipple)+avgRippleBeforeAfter);
                    end
                    if nRipplesStim > obj.minNripples && ~isempty(ripplesStim)
                        avgStim = nanmean(ripplesStim);
                        stdStim = nanstd(ripplesStim);
                        specStim = ft_specest_mtmfft(avgStim,[1:length(avgStim)]/obj.samplingRate,'freqoi',obj.freqoiForAvgSpec,'taper','hanning');
                        specStim = abs(squeeze(specStim(1,1,:,:)));
                    else
                        avgStim = nan(1,size(ripplesStim,2));
                        stdStim = nan(1,size(ripplesStim,2));
                        specStim = nan(1,length(obj.freqoiForAvgSpec));
                    end
                    
                    %get ripple centered TFRs (implemented in
                    %PACCalculator)
                    pacCalc = PACCalculator;
                    pacCalc.freqRange = obj.freqRangeForAvgSpec;
                    pacCalc.timeBeforeAfterEvent = obj.timeBeforeAfterEventRipSpec; %seconds
                    pacCalc.timeForBaseline = obj.timeForBaselineRip; %seconds, from Starestina et al
                    pacCalc.minNCycles = obj.minNCycles;
                    
                    meanTFRRipBefore = pacCalc.plotAvgSpecDiff(data, ripplesTimesBefore);
                    meanTFRRipStim = pacCalc.plotAvgSpecDiff(data, ripplesTimesStim);
                    
                    resultsPerChan(iChan).nRipplesBefore = nRipplesBefore;
                    resultsPerChan(iChan).nRipplesStim = nRipplesStim;
                    
                    resultsPerChan(iChan).avgBefore = avgBefore;
                    resultsPerChan(iChan).avgStim = avgStim;
                    resultsPerChan(iChan).stdBefore = stdBefore;
                    resultsPerChan(iChan).stdStim = stdStim;
                    resultsPerChan(iChan).specBefore = specBefore;
                    resultsPerChan(iChan).specStim = specStim;
                    resultsPerChan(iChan).meanTFRRipBefore = meanTFRRipBefore;
                    resultsPerChan(iChan).meanTFRRipStim = meanTFRRipStim;
                    
                    
                end
                
                results(iPatient).resultsPerChan = resultsPerChan;
            end
            
            if ~isempty(fileNameResults)
                save(fileNameResults,'results');
            end
            
        end
        
        % Based on plotResultsRipplesData from rippleDetector
        function plotResultsSpindlesData (obj, results, folderToSave)
            
            %The method receives as input a struct with the format of the output of runSpindlesData and produces figures.
            %The input folderToSave (optional) sets the folder into which the figures will be saved.
            
            if nargin < 3 || isempty(folderToSave)
                toSave = false;
            else
                toSave = true;
            end
            
            avgRippleBeforeAfter = obj.avgRippleBeforeAfter*obj.samplingRate;
            
            nPatients = length(results);
            
            for iPatient = 1:nPatients
                nChans = length(results(iPatient).resultsPerChan);
                
                for iChan = 1:nChans
                    %spindle summary figure - average spindle and spectogram
                    %before and during stimulation
                    f1 = newA4figure(['Spindles Data patient ',results(iPatient).patientName,' channel ',num2str(results(iPatient).resultsPerChan(iChan).channelNum),' Area ',results(iPatient).resultsPerChan(iChan).area]);
                    
%                     subplot(4, 3, [1 4]);
%                     meanSpecs = results(iPatient).resultsPerChan(iChan).meanSpecs;
%                     meanEpochs = results(iPatient).resultsPerChan(iChan).meanEpochs;
%                     imagesc(obj.xaxisForRipSp/obj.samplingRate, obj.freqRangeSpRip, meanSpecs(:,obj.specStartPointRipSp+1:end-obj.specStartPointRipSp));
%                     set(gca, 'YDir','normal');
%                     xlabel('Time (sec)');
%                     ylabel('Frequency (Hz)');
%                     colorbar;
%                     set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
%                     
%                     yyaxis right;
%                     if length(meanEpochs)>1
%                         plot(obj.xaxisForRipSp/obj.samplingRate,meanEpochs(obj.specStartPointRipSp+1:end-obj.specStartPointRipSp),'color','k','linewidth',1);
%                     end
%                     title({['Power spectrum around spindles, Before'],['nSpindles=',num2str(results(iPatient).resultsPerChan(iChan).nEpochs)]});
%                     ylabel('uV');
%                     
%                     subplot(4,3,[7 10]);
%                     currAngles = results(iPatient).resultsPerChan(iChan).SIanglesSpRip;
%                     if ~isempty(currAngles)
%                         circm = circ_mean(currAngles(~isnan(currAngles))');
%                         tmpplot = polarhistogram(currAngles(~isnan(currAngles)),obj.nBinsPolar,'Normalization','probability');
%                         hold all;
%                         polarplot([circm circm],[0 max(tmpplot.Values)/2],'color','r','linewidth',2);
%                         hold off;
%                         nSpindles = sum(~isnan(currAngles));
%                         title({['pval r = ',num2str(results(iPatient).resultsPerChan(iChan).r)], ['v = ', num2str(results(iPatient).resultsPerChan(iChan).v)], ['nSpindles = ',num2str(nSpindles)]});
%                     else
%                         title('No spindles');
%                     end
                     
                    subplot(4,3,2);
                    shadedErrorBar([-avgRippleBeforeAfter:avgRippleBeforeAfter]/obj.samplingRate,results(iPatient).resultsPerChan(iChan).avgBefore,results(iPatient).resultsPerChan(iChan).stdBefore/sqrt(results(iPatient).resultsPerChan(iChan).nRipplesBefore));
                    xlabel('Time (sec)');
                    ylabel('uV');
                    title(['Average of detected spindle - before, nSpindles=',num2str(results(iPatient).resultsPerChan(iChan).nRipplesBefore)]);
                    
                    subplot(4,3,3);
                    shadedErrorBar([-avgRippleBeforeAfter:avgRippleBeforeAfter]/obj.samplingRate,results(iPatient).resultsPerChan(iChan).avgStim,results(iPatient).resultsPerChan(iChan).stdStim/sqrt(results(iPatient).resultsPerChan(iChan).nRipplesStim));
                    xlabel('Time (sec)');
                    ylabel('uV');
                    title(['Average of detected spindle - stimulation, nSpindles=',num2str(results(iPatient).resultsPerChan(iChan).nRipplesStim)]);
                    
                    subplot(4,3,5);
                    plot(obj.freqoiForAvgSpec,results(iPatient).resultsPerChan(iChan).specBefore);
                    xlabel('Frequency');
                    title('Spectrum of average 0-10 Hz');
                    
                    subplot(4,3,6);
                    plot(obj.freqoiForAvgSpec,results(iPatient).resultsPerChan(iChan).specStim);
                    xlabel('Frequency');
                    title('Spectrum of average 0-10 Hz');
                    
                    
                    cmin = nanmin(nanmin(results(iPatient).resultsPerChan(iChan).meanTFRRipBefore(:)),nanmin(results(iPatient).resultsPerChan(iChan).meanTFRRipStim(:)));
                    cmax = nanmax(nanmax(results(iPatient).resultsPerChan(iChan).meanTFRRipBefore(:)),nanmax(results(iPatient).resultsPerChan(iChan).meanTFRRipStim(:)));
                    if isnan(cmin); cmin = 0; end
                    if isnan(cmax); cmax = inf; end
                                        
                    subplot(4,3,8);
                    if size(results(iPatient).resultsPerChan(iChan).meanTFRRipBefore,1)>1 && cmax > cmin 
                        imagesc([-obj.timeBeforeAfterEventRipSpec*obj.samplingRate:obj.timeBeforeAfterEventRipSpec*obj.samplingRate]/obj.samplingRate,obj.freqRangeForAvgSpec, results(iPatient).resultsPerChan(iChan).meanTFRRipBefore,[cmin cmax]);
                        set(gca, 'YDir','normal');
                        xlabel('Time (sec)');
                        ylabel('frequency');
                        title('Average of spindle triggerd TFR - before');
                        colorbar;
                    else
                        title('No spindles (before)');
                    end
                    
                    subplot(4,3,9);
                    if size(results(iPatient).resultsPerChan(iChan).meanTFRRipStim,1)>1
                        imagesc([-obj.timeBeforeAfterEventRipSpec*obj.samplingRate:obj.timeBeforeAfterEventRipSpec*obj.samplingRate]/obj.samplingRate,obj.freqRangeForAvgSpec, results(iPatient).resultsPerChan(iChan).meanTFRRipStim,[cmin cmax]);
                        set(gca, 'YDir','normal');
                        xlabel('Time (sec)');
                        ylabel('frequency');
                        title('Average of spindle triggerd TFR - stimulation');
                        colorbar;
                    else
                        title('No spindles (stim)');
                    end
                    
                    freqIndsSpindlesRange = find(obj.freqRangeForAvgSpec==obj.freqRangeForShowingSpindles(1)):find(obj.freqRangeForAvgSpec==obj.freqRangeForShowingSpindles(end));
                    if size(results(iPatient).resultsPerChan(iChan).meanTFRRipBefore,1)>1
                        TFRSmallBefore = results(iPatient).resultsPerChan(iChan).meanTFRRipBefore(freqIndsSpindlesRange,:);
                    else
                        TFRSmallBefore = nan;
                    end
                    if size(results(iPatient).resultsPerChan(iChan).meanTFRRipStim,1)>1
                        TFRSmallStim = results(iPatient).resultsPerChan(iChan).meanTFRRipStim(freqIndsSpindlesRange,:);
                    else
                        TFRSmallStim = nan;
                    end
                    
                    
                    suptitle(['Spindles Channel Data ',results(iPatient).patientName,' ',num2str(num2str(results(iPatient).resultsPerChan(iChan).channelNum)),' Area ',results(iPatient).resultsPerChan(iChan).area]);
                    
                    if toSave
                        set(f1, 'Position', get(0, 'Screensize'));
                        if ~isfolder([folderToSave,'\',results(iPatient).patientName])
                            mkdir([folderToSave,'\',results(iPatient).patientName]);
                        end
                        saveas(f1,[folderToSave,'\',results(iPatient).patientName,'\spindle_channel_data_' num2str(results(iPatient).resultsPerChan(iChan).channelNum),'.jpg']);
                        close(f1);
                    end
                end
            end
        end
        
        
        function plotPopulationFig(obj, runData, filename,figSuffix, outputFigureFolder, aspectRatio)
            SAVE_TABLE = 0;
            
            if exist('aspectRatio','var')
                TFR_aspectRatio = aspectRatio;
            else
                TFR_aspectRatio = 0;
            end
                
            mm = matfile(filename);
            
            spindleInfo = mm.spindleInfo;
            ch_area = spindleInfo.area;
            nRip = spindleInfo.nRip;
            nRip_stim = spindleInfo.nRip_stim;
            avgBefore = spindleInfo.avgBefore;
            stdBefore = spindleInfo.stdBefore;
            meanTFRRipBefore = spindleInfo.meanTFRRipBefore;
            meanTFRRipStim = spindleInfo.meanTFRRipStim;
            artifactInd = spindleInfo.artifactInd;
            avgStim = spindleInfo.avgStim;
            nRipThreshold = 20;
            
            nChan = size(avgBefore,1);
            
            for ii_a = 1:3
                
                clear fig_info areaInd
                
                if ii_a == 1
                    for ii = 1:length(ch_area)
                        area = ch_area{ii};
                        if area.isFrontal
                            areaInd(ii) = 1;
                        end
                    end
                    ind = find(areaInd);
                    figName = ['spindleGrandAverage_frontal_',figSuffix];
                elseif ii_a == 2
                    ind = 1:length(ch_area);
                    figName = ['spindleGrandAverage_all_channels_',figSuffix];
                elseif ii_a == 3
                     for ii = 1:length(ch_area)
                        area = ch_area{ii};
                        if ~area.isFrontal
                            areaInd(ii) = 1;
                        end
                    end
                    ind = find(areaInd);
                    figName = ['spindleGrandAverage_nonFrontal_',figSuffix];
                end
                
                ind(ismember(ind,artifactInd)) = [];
                ind(ismember(ind,find(nRip < nRipThreshold))) = [];
                
                NRip = sum(nRip(ind));
                
                % stats
                disp(figName)
                disp(sprintf('%d electrodes, %d pts, %d ripples',length(ind),length(unique(spindleInfo.ptNum_v(ind))),NRip))
                fig_info{1} = sprintf('%d electrodes, %d pts, %d spindles',length(ind),length(unique(spindleInfo.ptNum_v(ind))),NRip);
                disp(sprintf('%d channels in freq profile comparison', length(ind)))

                
                avgWeighted = zeros(1,length(avgBefore)); stdWeighted = zeros(1,length(avgBefore));
                freqN = size(squeeze(meanTFRRipBefore(1,:,:)),1);
                meanTFRRipWeighted = zeros(freqN,length(avgBefore));
                for iiC = ind
                    avgWeighted = avgWeighted + (nRip(iiC)/NRip) * avgBefore(iiC,:);
                    stdWeighted = stdWeighted + (sqrt( (nRip(iiC)-1)* stdBefore(iiC,:).^2 )/(NRip-iiC));
                    meanTFRRipWeighted(:,:) = meanTFRRipWeighted + (nRip(iiC)/NRip) * squeeze(meanTFRRipBefore(iiC,:,:)) ;
                end
                
                NRip_stim = sum(nRip_stim(ind));
                meanTFRRipStimWeighted = zeros(freqN,length(avgBefore));
                for iiC = ind
                     meanTFRRipStimWeighted(:,:) = meanTFRRipStimWeighted + (nRip_stim(iiC)/NRip) * squeeze(meanTFRRipStim(iiC,:,:)) ;
                end
                fig_info{2} = sprintf('%d spindles in TFR-stim',NRip_stim);

                
                f0 = figure('Name', figName,'NumberTitle','off');
                set(gcf,'DefaultAxesFontSize',8);
                set(gcf,'DefaultAxesFontName','arial');
                set(gcf,'PaperUnits','centimeters','PaperPosition',[0.2 0.2 15 9]); % this size is the maximal to fit on an A4 paper when printing to PDF
                set(gcf,'PaperOrientation','portrait');
                set(gcf,'Units','centimeters','Position', get(gcf,'paperPosition')+[1 1 0 0]);
                colormap('jet');
                
                if SAVE_TABLE
                    varNames = {'spindle_mean','spindle_sem'};
                    % data3 = meanTFRRipWeighted;
                    data2 = stdWeighted/sqrt(nChan);
                    data1 = avgWeighted;
                    data_table = table(data1(:),data2(:),'VariableNames',varNames);
                    outputPubFolder = 'C:\Users\mgeva\Documents\GitHub\closedLoop-pub\figureGen';
                    
                    if ii_a == 1
                        filename = 'ExtendedDataFig10b_data_table'; end
                    if ii_a == 2
                        filename = 'Fig2_spindle_data_table'; end
                    save(fullfile(outputPubFolder,filename),'data_table');
                end
                
                % Ripple waveform average
                axes('position',[0.12,0.62,0.3,0.3])
                shadedErrorBar(1:length(avgBefore),avgWeighted,stdWeighted/sqrt(nChan))
                axis tight
                set(gca,'xlim',[750 1250]);
                XLIM = get(gca,'xlim');
                
                line([XLIM(1),XLIM(1)]+10,[20 70],'color','k','linewidth',7)
                line([XLIM(1), XLIM(1)+100] ,[-80 -80],'color','k','linewidth',7)
                axis off
                
                % Average TFR
                if ~TFR_aspectRatio
                    axes('position',[0.12,0.1,0.4,0.4])
                else
                   axes('position',[0.12,0.1,0.25,0.4])
                end
                meanTFRRipWeighted(isnan(meanTFRRipWeighted)) = 0;
                M = ceil(max(max(meanTFRRipWeighted)));
                h = imagesc(meanTFRRipWeighted, [-M,M]);
                set(gca,'xtick',[1,1000,2000],'xticklabels',[-1,0,1])
                freq_labels = [obj.freqRangeForAvgSpec(1),15, obj.freqRangeForAvgSpec(end)];
                yaxis_base = 5;
                set(gca,'ytick',[1,15-yaxis_base,freqN],'YTickLabel',freq_labels,'TickDir','out','TickLength',[0.0100 0.030])
                hold all; 
                LW = 1.2
                plot(get(gca,'xlim'),[9 9]-yaxis_base,'k--','linewidth',LW)
                plot(get(gca,'xlim'),[16 16]-yaxis_base,'k--','linewidth',LW)
                axis xy
                cc = colorbar;
                cc.Ticks = [-M,0,M];
                cc.TickLabels = [-M,0,M]*100;
                box off
                
                % Average TFR during stimulation blocks (immediate effect)
                if ~TFR_aspectRatio
                    axes('position',[0.58,0.1,0.4,0.4])
                else
                    axes('position',[0.58,0.1,0.25,0.4])
                end
                
                meanTFRRipStimWeighted(isnan(meanTFRRipStimWeighted)) = 0;
                M = ceil(max(max(meanTFRRipStimWeighted)));
                h = imagesc(meanTFRRipStimWeighted, [-M,M]);
                set(gca,'xtick',[1,1000,2000],'xticklabels',[-1,0,1])
                freq_labels = [obj.freqRangeForAvgSpec(1),15, obj.freqRangeForAvgSpec(end)];
                set(gca,'ytick',[1,15-yaxis_base,freqN],'YTickLabel',freq_labels,'TickDir','out','TickLength',[0.0100 0.030])
                hold all; 
                plot(get(gca,'xlim'),[9 9]-yaxis_base,'k--','linewidth',LW)
                plot(get(gca,'xlim'),[16 16]-yaxis_base,'k--','linewidth',LW)
                axis xy                
                cc = colorbar;
                cc.Ticks = [-M,0,M];
                cc.TickLabels = [-M,0,M]*100;
                
                
                [spectrum,ntaper,freqoi]  = ft_specest_mtmfft(avgBefore(ind,500:1500),[1:length(avgBefore)]/obj.samplingRate,'freqoi',[0.5:0.1:20],'taper','dpss','tapsmofrq',8);
                specBefore = 10*log10(abs(squeeze(spectrum(1,1,:)).^2));
                % specBefore = specBefore/max(specBefore);
                
                [spectrum_stim,ntaper_stim,freqoi_stim]  = ft_specest_mtmfft(avgStim(ind,500:1500),[1:length(avgStim)]/obj.samplingRate,'freqoi',[0.5:0.1:20],'taper','dpss','tapsmofrq',8);
                specStim = 10*log10(abs(squeeze(spectrum_stim(1,1,:)).^2));
              
                axes('position',[0.62,0.62,0.15,0.35])
                plot(freqoi, specBefore', 'k','linewidth',2)
                set(gca,'ylim',[-10 45],'ytick',[-10 0 45])    
                box off
                [val, ind] = max(specBefore);
                fig_info{2} = sprintf('max at %2.2f(Hz)',freqoi(ind));
                hold on
                plot(freqoi, specStim', 'r','linewidth',2)
                
                XLIM = get(gca,'xlim');
                YLIM = get(gca,'ylim');

                text(XLIM(2)+diff(xlim)/10,YLIM(1)+diff(YLIM)/2,fig_info,'fontsize',5)
                
                res =  600;
                a = gcf;
                eval(['print ', [outputFigureFolder,'\',a.Name], ' -f', num2str(a.Number),sprintf(' -dtiff  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!
            end
            
            
        end % plot func 
        
        
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
