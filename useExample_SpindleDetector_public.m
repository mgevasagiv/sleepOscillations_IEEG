%the struct runData holds data about patients and where the different event
%types are stored

% ---- UPDATE this part -

% the main path for extracted data here -
% for the given example, it's in the same folder as this code:
data_p_path = fullfile(fileparts(which('useExample_SpindleDetector_public')),'example');

% the code assumes extracted data will be found under
% runData(iPatient).DataFolder = fullfile(data_p_path,patients{iPatient},expNames{iPatient},'MACRO');

% We will save figures in the same example folder:
outputFolder = fullfile(data_p_path,'output','figures'); 
if ~exist(outputFolder,'dir'),mkdir(outputFolder);end

% subject name
patients = {'p1'};
% session name
expNames = {'EXP1'};
% sleep-scoring vector ('1' for NREM epochs to analyze)
sleepScoreFileName = {'sleepScore_p1'};
% channel id to analyze
channelsPerPatient = {22};

% macroMontageFileName contains channel ids and area names per subject

% ----------------------

%building the run data, not that for all file names of detections the
%methods assume the name is <provided name according to runData><channel
%num>, e.g. if runData(iPatient).SWStaresinaFileName='c:\slow_wave', then
%the slow waves file for channel 1 is 'c:\slow_wave1.mat'.
runData = [];
nPatients = length(patients);
for iPatient = nPatients :-1:1
    runData(iPatient).patientName = patients{iPatient};
    runData(iPatient).expNames = expNames{iPatient};

    %The folder where the raw data is stored - you will need to change it
    runData(iPatient).DataFolder = fullfile(data_p_path,patients{iPatient},expNames{iPatient},'MACRO');
    
    runData(iPatient).MicroDataFolder = fullfile(data_p_path, patients{iPatient},expNames{iPatient},'MICRO');    
    runData(iPatient).microChannelsFolderToLoad = runData(iPatient).MicroDataFolder;
  
    
    %The folder+filename into which the spikes results is going to be stored or is
    %already stored if the spikes detection was already run (the folder should
    %exist)
    macroSpikeFolder = fullfile(data_p_path,patients{iPatient},expNames{iPatient},'MACRO');
    if ~exist(macroSpikeFolder,'dir'),mkdir(macroSpikeFolder);end
    runData(iPatient).SpikesFileNames = fullfile(macroSpikeFolder,...
            sprintf('MacroInterictalSpikeTimesFor_%s_%s_',patients{iPatient},expNames{iPatient}));
    
    SW_folder = fullfile(data_p_path,patients{iPatient},expNames{iPatient},'MACRO','SWStaresinaResults');
    runData(iPatient).SWStaresinaFileName = fullfile(SW_folder,'SWTimes');
    if ~exist(SW_folder,'dir')
        mkdir(SW_folder)
    end
    
    %The folder+filename into which the staresina spindle detections
    %results is going to be stored (the folder should exist)
    spindleFolder = fullfile(data_p_path,patients{iPatient},expNames{iPatient},'MACRO','spindleStaresinaResults');
    runData(iPatient).SpindlesStaresinaFileNames = fullfile(spindleFolder,'spindleTimes');
    if ~exist(spindleFolder,'dir')
        mkdir(spindleFolder)
    end
    
    spindleFolder = fullfile(data_p_path,patients{iPatient},expNames{iPatient},'MACRO','spindleResults');
    runData(iPatient).HighFreqSpindlesFileNames = fullfile(spindleFolder,'highFreqSpindleTimes');
    if isempty(dir(spindleFolder))
        mkdir(spindleFolder)
    end

    %The folder+filename from which spindles are going to be loaded (should be the same
    %as the line above if the detections are first run and saved into SpindlesStaresinaFileNames
    spindleFolder = fullfile(data_p_path,patients{iPatient},expNames{iPatient},'MACRO','spindleResults');
    runData(iPatient).SpindlesFileNames = fullfile(spindleFolder,'spindleTimes');
    if ~exist(spindleFolder,'dir')
        mkdir(spindleFolder)
    end
     
    %name of the EXP data for the patient
    runData(iPatient).ExpDataFileName = fullfile(data_p_path,patients{iPatient},expNames{iPatient},[patients{iPatient},'_',expNames{iPatient},'_dataset.mat']);

    %name of the sleep scoring mat file for the patient
    runData(iPatient).sleepScoringFileName = fullfile(runData(iPatient).DataFolder,[sleepScoreFileName{iPatient},'.mat']);
    
    runData(iPatient).macroMontageFileName = fullfile(data_p_path,'MACRO_MONTAGE',patients{iPatient},expNames{iPatient},'MacroMontage.mat');
    
    runData(iPatient).spikeData = fullfile(data_p_path,patients{iPatient},expNames{iPatient},'averagedRef',patients{iPatient},'_spike_timestamps_post_processing.mat');

    %channels that the detections will be performed on
    runData(iPatient).channelsToRunOn = channelsPerPatient{iPatient};   

end


%% an example for detecting spindles directly using SpindleDetectorClass (it's the same thing the wrapper below does in batch)
sd = SpindleDetectorClass;

%an example of using the ripple detection directly and not with the wrapper
%(on the first channel of the first patient for this example)
iPatient = 1;
% here focusing on channels selected forcouple analysis
currChan = runData(iPatient).channelsToRunOn(1);

%loading - sleep scoring, IIS, data
sleepScoring = load(runData(1).sleepScoringFileName);
sleepScoring = sleepScoring.sleep_score_vec;

%% load or perform interictal Spikes Detection
try
    peakTimes = load([runData(iPatient).SpikesFileNames,num2str(currChan),'.mat']);
    peakTimes = peakTimes.peakTimes;
catch
    peakTimes = [];
end

currData = load(fullfile(runData(iPatient).DataFolder,['CSC',num2str(currChan),'.mat']));
currData = currData.data;
%detecting the spindles
returnStats = 1;
sd.spindleRangeMin = 11;
[spindlesTimes,spindleStats,spindlesStartEndTimes] = sd.detectSpindles(currData, sleepScoring, peakTimes, returnStats);

%plotting the single spindles and saving the figures
sd.plotSpindlesSimple(currData, spindlesTimes, outputFolder)

% scroll through spindles and their spectrograms using any key
blockSize = 4;
sd.plotSpindles(currData,spindlesTimes,blockSize);

%% an example for saving ripples using the wrapper AnalyzeSleepOsc.saveDetectionResults
% Can be downloaded from https://github.com/mgevasagiv/rippleDetection_IEEG
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

