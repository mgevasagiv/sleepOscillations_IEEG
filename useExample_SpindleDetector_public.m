%the struct runData holds data about patients and where the different event
%types are stored

data_p_path = 'E:\Data_p\';

patients = {'p1'};
expNames = {'EXP4'};
sleepScoreFileName = {'sleepScore_p1'};
channelsPerPatient = {[22]};


%noisy micro channels
noisyChannelsPerPatient = {[]};

%building the run data, not that for all file names of detections the
%methods assume the name is <provided name according to runData><channel
%num>, e.g. if runData(iPatient).SWStaresinaFileName='c:\slow_wave', then
%the slow waves file for channel 1 is 'c:\slow_wave1.mat'.
runData = [];
nPatients = length(patients);
for iPatient = 1:nPatients 
    runData(iPatient).patientName = patients{iPatient};
    runData(iPatient).expNames = expNames{iPatient};

    %The folder where the raw data is stored - you will need to change it
    runData(iPatient).DataFolder = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO'];
    
    runData(iPatient).MicroDataFolder = [data_p_path, patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MICRO'];    
    runData(iPatient).microChannelsFolderToLoad = runData(iPatient).MicroDataFolder;
  
    
    %The folder+filename into which the spikes results is going to be stored or is
    %already stored if the spikes detection was already run (the folder should
    %exist)
    macroSpikeFolder = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO\'];
    runData(iPatient).SpikesFileNames = fullfile(macroSpikeFolder,...
            sprintf('MacroInterictalSpikeTimesFor_%s_%s_',patients{iPatient},expNames{iPatient}));
    
    SW_folder = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO\SWStaresinaResults'];
    runData(iPatient).SWStaresinaFileName = fullfile(SW_folder,'SWTimes');
    if isempty(dir(SW_folder))
        mkdir(SW_folder)
    end
    
    %The folder+filename into which the staresina spindle detections
    %results is going to be stored (the folder should exist)
    spindleFolder = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO\spindleStaresinaResults'];
    runData(iPatient).SpindlesStaresinaFileNames = fullfile(spindleFolder,'spindleTimes');
    if isempty(dir(spindleFolder))
        mkdir(spindleFolder)
    end
    
    spindleFolder = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO\spindleResults'];
    runData(iPatient).HighFreqSpindlesFileNames = fullfile(spindleFolder,'highFreqSpindleTimes');
    if isempty(dir(spindleFolder))
        mkdir(spindleFolder)
    end

    %The folder+filename from which spindles are going to be loaded (should be the same
    %as the line above if the detections are first run and saved into SpindlesStaresinaFileNames
    spindleFolder = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO\spindleResults'];
    runData(iPatient).SpindlesFileNames = fullfile(spindleFolder,'spindleTimes');
    if isempty(dir(spindleFolder))
        mkdir(spindleFolder)
    end
     
    %name of the EXP data for the patient
    runData(iPatient).ExpDataFileName = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\',patients{iPatient},'_',expNames{iPatient},'_dataset.mat'];
    
    % Extract stimulation times
    mm = matfile(runData(iPatient).ExpDataFileName); EXP_DATA = mm.EXP_DATA;
    runData(iPatient).stimulation_times = EXP_DATA.stimTiming.validatedTTL_NLX;
    clear mm EXP_DATA
    
    %name of the sleep scoring mat file for the patient
    runData(iPatient).sleepScoringFileName = [runData(iPatient).DataFolder,'\',sleepScoreFileName{iPatient},'.mat'];
    
    %extra fields for the micro coupling analysis
    runData(iPatient).MicroDataFolder = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MICRO\'];
    %The folder+filename from which micro ripples are going to be loaded
    runData(iPatient).MicroRipplesFileNames = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MICRO\rippleResults\rippleTimes'];
    runData(iPatient).microMontageFileName = [data_p_path,'MACRO_MONTAGE','\',patients{iPatient},'\',expNames{iPatient},'\montage.mat'];
    runData(iPatient).macroMontageFileName = [data_p_path,'MACRO_MONTAGE','\',patients{iPatient},'\',expNames{iPatient},'\MacroMontage.mat'];
    runData(iPatient).noisyChannels = noisyChannelsPerPatient{iPatient};
    
    runData(iPatient).spikeData = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\averagedRef\',patients{iPatient},'_spike_timestamps_post_processing.mat'];

    %channels that the detections will be performed on
    runData(iPatient).channelsToRunOn = channelsPerPatient{iPatient};   

end


%% an example for detecting spindles directly using SpindleDetector (it's the same thing the wrapper does inside)
sd = SpindleDetector;

%an example of using the ripple detection directly and not with the wrapper
%(on the first channel of the first patient for this example)
iPatient = 1;
% here focusing on channels selected forcouple analysis
currChan = runData(iPatient).channelsToRunOn(1);

%loading - sleep scoring, IIS, data
sleepScoring = load(runData(1).sleepScoringFileName);
sleepScoring = sleepScoring.sleep_score_vec;
peakTimes = load([runData(iPatient).SpikesFileNames,num2str(currChan),'.mat']);
peakTimes = peakTimes.peakTimes;
currData = load([runData(iPatient).DataFolder,'\CSC',num2str(currChan),'.mat']);
currData = currData.data;
%detecting the spindles
returnStats = 1;
sd.spindleRangeMin = 11;
[spindlesTimes,spindleStats,spindlesStartEndTimes] = sd.detectSpindles(currData, sleepScoring, peakTimes, returnStats);

%plotting the single spindles and saving the figures
sd.plotSpindles(currData,spindlesTimes);
outputFolder = 'E:\Data_p\ClosedLoopDataset\spindleDetResults';

%% spindle related analyses

%% This analysis is focused on channels specifically used for triple synchrony analysis - to generate sup fig 8 - 
% batch option

%% an example for saving ripples using the wrapper AnalyzeSleepOsc.saveDetectionResults
as = AnalyzeSleepOsc;

%setting which detections to run - 
whatToRun.runSpikes = false;
whatToRun.runRipples = false;
whatToRun.runRipplesBiPolar = false;
whatToRun.runSpindles = true;
whatToRun.runSpindlesStaresina = false;
whatToRun.runSWStaresina = false;
whatToRun.runSWMaingret = false;
whatToRun.HighFreqSpindles = true;

%saving detections (in this example, bipolar ripples detection)
parfor ii = 1:length(runData)
    as.saveDetectionResults(runData(ii), whatToRun);
end

%%
sd = SpindleDetector;

figSuffix = 'highFreqRange';
filename = fullfile(HighFreqSpindlesFolder,'spindleInfo_highFreq');
outputFigureFolder = 'E:\Figures\';
sd.plotPopulationFig(runData, filename, figSuffix, outputFigureFolder)

outputFigureFolder = 'E:\Figures\';
figSuffix = 'fullFreqRange_allCh';
filename = fullfile(fullFreqRangeSpindles_allChannelsFolder,'spindleInfo_fullFreq_allChanels');
aspectRatio = true;
sd.plotPopulationFig(runData, filename, figSuffix, outputFigureFolder, aspectRatio)


