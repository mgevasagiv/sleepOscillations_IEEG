%the struct runData holds data about patients and where the different event
%types are stored

data_p_path = 'E:\Data_p\';

patients = {'p488','p487','p489','p486','p499','p485','p490','p496','p498','p510','p505','p515','p510','p520','p497','p538','p541','p544','p545'};
expNames = {'EXP4','EXP3','EXP3','EXP8','EXP8','EXP8','EXP3','EXP8','EXP3','EXP7','EXP1','EXP9','EXP1','EXP3','EXP2','EXP3','EXP8','EXP8','EXP3'};
sleepScoreFileName = {'sleepScore_manualValidated_p488_4_LPHG2',...
    'sleepScore_manualValidated_p487_3_ROF2',...
    'sleepScore_manualValidated_p489_3_RPHG3',...
    'sleepScore_manualValidated_p486_8_RA1',...
    'sleepScore_manualValidated_p499_8_RPT2',...
    'sleepScore_manualValidated_p485_8_ROF7_0016',...
    'sleepScore_manualValidated_p490_3_RSTG1',...
    'sleepScore_manualValidated_p496_8_RSO4',...
    'sleepScore_manualValidated_p498_3_ROF1',...
    'sleepScore_manualValidated_p510_7_RPSM2',...
    'sleepScore_manualValidated_p505_1_RSTG1',...
    'sleepScore_manualValidated_p515_9_LOF3',...
    'sleepScore_manualValidated_p510_1_RPSM2',...
    'sleepScore_manualValidated_p520_3_CSC26',...
    'sleepScore_manualValidated_p497_2_RPHG4',...
    'sleepScore_manualValidated_p538_3_LOPR1',...
    'sleepScore_manualValidated_p541_8_LMI5',...
    'sleepScore_manualValidated_p544_8_LOF4',...
    'sleepScore_manualValidated_p545_3_LOF-AC8'};


%note macro montage is not necessary to run AnalyzeCoupling
% macroMontageNames = {'MM488','MM487','MM489','MM486','MM499','MM485','MM490','MM496','MM498','MM510','MM505','MM515','MM510_1'};

%for bipolar ripple detection - in every row the first index is the channel in which ripple
%detection is required and the second is the reference channel
% patients = {'p488','p487','p489','p486','p499','p485','p490','p496','p498','p510','p505','p515','p510','p520','p497'};
biPolarCouplesPerPatient = {[78 81;1 4;23 26],... % 488
                            [16 18; 10 11; 39 41; 31 32; 60 62],... % 487
                            [1 4; 8 11; 36 39; 29 33; 44 46],...% 489
                            [8 11; 46 47],... % p486
                            [46 49; 38 40; 24 25],... % p499
                            [52 54],... % p485
                            [1 3; 8 10; 50 53],... % 490 - RANH - excluded from analysis
                            [1 3],...% 496
                            [32 35; 33 35; 39 42; 40 42; 47 50; 48 50],...%498
                            [9 13; 47 51],... %510
                            [2 4],... 505
                            [9 12; 39 42; 40 42; 65 67; 66 67],... % 515
                            [9 13; 47 51],... %510
                            [17 18; 9 10],...% 520
                            [9 10; 16 19; 1 2; 40 41; 48 49],...% 497
                            [9 12],...%538 - verified - very nice, 
                            [1 2],... % Too noisy - [10 12; 47 50; 100 101],...%541 - verified
                            [1 2],... % No ripples on the other - [9 11, 29 30, 36 38],... % 544 - verified
                            [1 2;48 49; 32 33]}; % 1 looks artifact, 32 IED ;  % 545 - verified

%the couples on which the analysis will be run (MTL - 1st index, frontal -
%2nd index)
channelCouplesPerPatient = {[78 55; 1 10; 23 10; 1 62; 23 62],... % 488
                            [16 25; 10 25; 39 47; 39 67; 31 47; 31 67; 60 47; 60 67],... % 487
                            [1 18; 8 18; 29 53; 44 53],... % 489 % LMH removed - high IED levels, cortical areas - choosing away from IED contacts
                            [8 27; 8 22; 8 35; 46 58; 46 72; 46 77],... % 486
                            [46 55; 38 55; 24 17; 24 65],... % 499
                            [52 71],... % 485
                            [1 43],[1 54],... % 490, 496
                            [32 57; 40 57; 48 57],...%498
                            [9 21; 47 55; 47 79],... % 510
                            [2 16],... % 505
                            [9 24; 39 77; 39 47; 65 73; 65 47],...% 515
                            [9 21; 47 55; 47 79],... % 510
                            [17 25; 9 25],... % 520 
                            [9 24; 1 24; 16 24; 48 56],... %p497 
                            [9 22],... % 538
                            [1 28],... % 541                            
                            [1 17],... % 544
                            [48, 68]}; % 545

%channels on which detections will be performed (just an example)
% channelsPerPatient = {[2 23], [9 17 32 61], [2 8 29 37],[9],[39 23 47],[53],[],[],[],[47 48],[],[],[],[]};
channelsPerPatient = {[], [], [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]};

%probe chan is required for up vs down comparison
probeChans = [78, 16, 1, 46, 46, 8, nan, 1, 32, 9, 1,9,9, 17,9,...
                1,48,9,1];

%noisy micro channels
noisyChannelsPerPatient = {[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]};

INCLUDE_BOTH_HEM = 1;
%building the run data, not that for all file names of detections the
%methods assume the name is <provided name according to runData><channel
%num>, e.g. if runData(iPatient).SWStaresinaFileName='c:\slow_wave', then
%the slow waves file for channel 1 is 'c:\slow_wave1.mat'.
runData = [];
nPatients = length(patients);
for iPatient = 1:nPatients 
    runData(iPatient).patientName = patients{iPatient};
    runData(iPatient).expNames = expNames{iPatient};

    %the probe channel is required only for the up vs down comparison (can
    %be left empty otherwise)
    runData(iPatient).probeChan = probeChans(iPatient);
    %The folder where the raw data is stored - you will need to change it
    runData(iPatient).DataFolder = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO'];
    
    runData(iPatient).MicroDataFolder = [data_p_path, patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MICRO'];    
    runData(iPatient).microChannelsFolderToLoad = runData(iPatient).MicroDataFolder;
    % runData(iPatient).DataFolder = runData(iPatient).MicroDataFolder;
  
    
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
    
    %The folder+filename into which the bipolar ripples detections results is going
    %to be stored
    rippleFolder = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO\rippleBipolarResults'];
    runData(iPatient).RipplesBipolarFileNames = fullfile(rippleFolder,'rippleTimes');
    if isempty(dir(rippleFolder))
        mkdir(rippleFolder)
    end
    
    %-% USing bipolar ripples for analysis %-%
    %The folder+filename from which ripples are going to be loaded (should be the same
    %as the line above if the bipolar detections are first run and saved into
    %RipplesBipolarFileNames)
    %     rippleFolder = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO\rippleResults'];
    %     runData(iPatient).RipplesFileNames = fullfile(rippleFolder,'rippleTimes');
    %     if isempty(dir(rippleFolder))
    %         mkdir(rippleFolder)
    %     end
    
    runData(iPatient).RipplesFileNames = runData(iPatient).RipplesBipolarFileNames;
    
    %list of couples for bipolar ripple detection - where each row has the channel
    %in which the detection is performed in the first index, and the
    %reference channel in the second
    runData(iPatient).biPolarCouples = biPolarCouplesPerPatient{iPatient};
   
    %name of the EXP data for the patient
    runData(iPatient).ExpDataFileName = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\',patients{iPatient},'_',expNames{iPatient},'_dataset.mat'];
    
    % Extract stimulation times
    mm = matfile(runData(iPatient).ExpDataFileName); EXP_DATA = mm.EXP_DATA;
    runData(iPatient).stimulation_times = EXP_DATA.stimTiming.validatedTTL_NLX;
    clear mm EXP_DATA
    
    %name of the sleep scoring mat file for the patient
    runData(iPatient).sleepScoringFileName = [runData(iPatient).DataFolder,'\',sleepScoreFileName{iPatient},'.mat'];
    %couples for the coupling analysis
    runData(iPatient).channelCouples = channelCouplesPerPatient{iPatient};
    
    %extra fields for the micro coupling analysis
    runData(iPatient).MicroDataFolder = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MICRO\'];
    %The folder+filename from which micro ripples are going to be loaded
    runData(iPatient).MicroRipplesFileNames = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MICRO\rippleResults\rippleTimes'];
    runData(iPatient).microMontageFileName = [data_p_path,'MACRO_MONTAGE','\',patients{iPatient},'\',expNames{iPatient},'\montage.mat'];
    runData(iPatient).macroMontageFileName = [data_p_path,'MACRO_MONTAGE','\',patients{iPatient},'\',expNames{iPatient},'\MacroMontage.mat'];
    runData(iPatient).noisyChannels = noisyChannelsPerPatient{iPatient};
    
    runData(iPatient).spikeData = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\averagedRef\',patients{iPatient},'_spike_timestamps_post_processing.mat'];

    %channels that the detections will be performed on
    runData(iPatient).channelsToRunOn = unique(channelCouplesPerPatient{iPatient}(:,2)');   
    
    % Leave cortical/MTL channels only for the same hemisphere
    if sum(strcmp(runData(iPatient).patientName,{'p485','p490'}))
        continue
    end
    
    mm = matfile(runData(iPatient).macroMontageFileName);
    MacroMontage = mm.MacroMontage;
    probe_str = MacroMontage(runData(iPatient).probeChan).Area;
    hemisphere_str = probe_str(1);
    rmv_idx = []; 
    for ss = 1:size(runData(iPatient).channelCouples,1)
        area_str = MacroMontage(runData(iPatient).channelCouples(ss,1)).Area;
        if ~strcmp(area_str(1),hemisphere_str)
            rmv_idx = [rmv_idx ss];
        end
    end
    
    if ~INCLUDE_BOTH_HEM
        runData(iPatient).channelCouples(rmv_idx,:) = [];
    else
        runData(iPatient).otherHemCouples = rmv_idx;
    end

end


%% an example for saving spindles using the wrapper AnalyzeCoupling.saveDetectionResults
ac = AnalyzeCoupling;

%setting which detections to run - in this example, standard range spindles
%and high-freq spindle would be saved for each pt
whatToRun.runSpikes = false;
whatToRun.runRipples = false;
whatToRun.runRipplesBiPolar = false;
whatToRun.runSpindles = false;
whatToRun.runSpindlesStaresina = false;
whatToRun.runSWStaresina = false;
whatToRun.runSWMaingret = false;
whatToRun.HighFreqSpindles = true;

%saving detections (in this example, spindles detection)
parfor ii = 1:length(runData)
    ac.saveDetectionResults(runData(ii), whatToRun);
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

% Gen spindle data summary per patient\channel - separately for full range
% of spindle detection and high-freq detection
clear whatToRun;
whatToRun.HighFreqSpindles = 1;
HighFreqSpindlesFolder = 'E:\Data_p\ClosedLoopDataset\spindleDetResults\MACRO_fastSpindles\';
figureFolder = 'E:\Data_p\ClosedLoopDataset\spindleDetResults\spindlesFiguresFastSpindles';
parfor ii = 1:length(runData)

    resultsFile = fullfile(HighFreqSpindlesFolder,...
        sprintf('resultsSpindlesMACRO_%s_%s.mat',runData(ii).patientName,expNames{ii} ));
    runData(ii).channelsToRunOn = unique(channelCouplesPerPatient{ii}(:,2)');   
    if isempty(dir(resultsFile))
        resultsData = sd.runSpindleData(runData(ii),resultsFile,whatToRun);
    end
    
end
for ii = 1:length(runData)
    resultsFile = fullfile(HighFreqSpindlesFolder,...
        sprintf('resultsSpindlesMACRO_%s_%s.mat',runData(ii).patientName,expNames{ii} ));
    
    disp(['plotting results ',runData(ii).patientName])
    mm = matfile(resultsFile);
    resultsData = mm.results;
    
    sd.plotResultsSpindlesData(resultsData,figureFolder);
    
end


clear whatToRun;
whatToRun.HighFreqSpindles = 0;
fullFreqRangeSpindlesFolder = 'E:\Data_p\ClosedLoopDataset\spindleDetResults\MACRO_allSpindles\';
figureFolder = 'E:\Data_p\ClosedLoopDataset\spindleDetResults\spindlesFiguresAllSpindles';
parfor ii = 1:length(runData)
    
    resultsFile = fullfile(fullFreqRangeSpindlesFolder,...
        sprintf('resultsSpindlesMACRO_%s_%s.mat',runData(ii).patientName,expNames{ii} ));
    runData(ii).channelsToRunOn = unique(channelCouplesPerPatient{ii}(:,2)');   
    if isempty(dir(resultsFile))
        resultsData = sd.runSpindleData(runData(ii),resultsFile,whatToRun);
    end
    
end
for ii = 1:length(runData)
    
    resultsFile = fullfile(fullFreqRangeSpindlesFolder,...
        sprintf('resultsSpindlesMACRO_%s_%s.mat',runData(ii).patientName,expNames{ii} ));
    
    disp(['plotting results ',runData(ii).patientName])
    mm = matfile(resultsFile);
    resultsData = mm.results;
    
    sd.plotResultsSpindlesData(resultsData,figureFolder);
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load SW/spindle indices - used for whole brain effects (Fig 2) to create
% supplementary figure for this channel group
summaryFile = 'E:\Data_p\ClosedLoopDataset\stimEffectResults\allContactsStimResults.mat';
mm = matfile(summaryFile);
corticalStimEffectIndices = mm.corticalStimEffectIndices;
runData = mm.runData;
clear whatToRun;
whatToRun.HighFreqSpindles = 0;
fullFreqRangeSpindlesFolder = 'E:\Data_p\ClosedLoopDataset\spindleDetResults\MACRO_allFreqSpindles_allChannels\';
figureFolder = 'E:\Data_p\ClosedLoopDataset\spindleDetResults\spindlesFigures_allFreqSpindles_allChannels';
pt_list = unique(corticalStimEffectIndices.pt_vec);

for iii = 1:length(corticalStimEffectIndices.exp_vec)
    A = corticalStimEffectIndices.exp_vec;
    exp_vec(iii) =  str2num(A(iii,end));
end
for ii = 1:length(runData)
    
    resultsFile = fullfile(fullFreqRangeSpindlesFolder,...
        sprintf('resultsSpindlesMACRO_%s_%s.mat',runData(ii).patientName,runData(ii).EXP ));
    
    pt_num = str2num(runData(ii).patientName(2:end));
    exp_num = str2num(runData(ii).EXP);
    pt_ind = find(corticalStimEffectIndices.pt_vec == pt_num);
    exp_ind = find(exp_vec == exp_num);
    entry_ind = exp_ind(ismember(exp_ind,pt_ind));
    
    channel_ind = corticalStimEffectIndices.chNum_vec(entry_ind);
    runData(ii).channelsToRunOn = unique(channel_ind);   
    if isempty(dir(resultsFile))
        resultsData = sd.runSpindleData(runData(ii),resultsFile,whatToRun);
    end
    
end

FULL_NIGHT_PT = 19;
for ii = 1:length(runData(1:FULL_NIGHT_PT)) % remove nap patients
    
    resultsFile = fullfile(fullFreqRangeSpindlesFolder,...
        sprintf('resultsSpindlesMACRO_%s_%d.mat',runData(ii).patientName,str2num(runData(ii).EXP) ));
    
    disp(['plotting results ',runData(ii).patientName])
    mm = matfile(resultsFile);
    resultsData = mm.results;
    
    sd.plotResultsSpindlesData(resultsData,figureFolder);
    
end



%% collect population figure for all patients - this is for channels use for triple analysis

% These folders hold only channels included in couple-analysis for triple
% ripple-spindle-sw analysis a la Maingret)
HighFreqSpindlesFolder = 'E:\Data_p\ClosedLoopDataset\spindleDetResults\MACRO_fastSpindles\';
fullFreqRangeSpindlesFolder = 'E:\Data_p\ClosedLoopDataset\spindleDetResults\MACRO_allSpindles\';

% This folder holds all channels participating in stimulation-effect
% analyis (inclusive list)
fullFreqRangeSpindles_allChannelsFolder = 'E:\Data_p\ClosedLoopDataset\spindleDetResults\MACRO_allFreqSpindles_allChannels\';

for ii_a = 1:3
    if ii_a == 1
        dataFolder = fullFreqRangeSpindlesFolder;
        filename = 'spindleInfo_fullFreq';
    elseif ii_a == 2
        dataFolder = HighFreqSpindlesFolder;
        filename = 'spindleInfo_highFreq';
    elseif ii_a == 3
        dataFolder = fullFreqRangeSpindles_allChannelsFolder;
        filename = 'spindleInfo_fullFreq_allChanels';
        
    end
    
    cnt = 0;
    nRipThreshold = -1;  % do not require ripples bfr stim blocks because for some pts this block is out of nrem sleep
    nRip = [];
    ptNum_v = []; chNum_v = [];
    vtest_v = []; rtest_v = [];
    chMTL_area = cell(1,1);
    avgBefore = zeros(1,2001);
    stdBefore = zeros(1,2001);
    nSU = 0;
    clear avgBefore stdBefore meanTFRRipBefore

    for iiP = 1:length(runData)
        if ii_a == 3
        resultsFile = fullfile(dataFolder,...
            sprintf('resultsSpindlesMACRO_%s_%s.mat',runData(iiP).patientName,expNames{iiP}(end) ));
        else
            resultsFile = fullfile(dataFolder,...
            sprintf('resultsSpindlesMACRO_%s_%s.mat',runData(iiP).patientName,expNames{iiP} ));
        end
        resultsData = load(resultsFile, 'results');
        resultsData = resultsData.results;
        ptNum = str2num(resultsData.patientName(2:end));
        
        for iiC = 1:length(resultsData.resultsPerChan)
            if resultsData.resultsPerChan(iiC).nRipplesBefore > nRipThreshold
                if mean(resultsData.resultsPerChan(iiC).stdBefore) > 500
                    continue
                end
                
                if resultsData.resultsPerChan(iiC).nRipplesBefore ~= 0
                    A = resultsData.resultsPerChan(iiC).avgBefore(995:1005);
                else
                    A = resultsData.resultsPerChan(iiC).avgStim(995:1005);
                end
                [pks,locs,w,p]  = findpeaks(A,'SortStr','descend','NPeaks',1);
                
                if isempty(locs)
                    continue
                end
                
                if resultsData.resultsPerChan(iiC).nRipplesBefore ~= 0
                    cnt = cnt + 1;

                    D1 = resultsData.resultsPerChan(iiC).avgBefore;
                    D2 = resultsData.resultsPerChan(iiC).stdBefore;
                    
                    D3 = resultsData.resultsPerChan(iiC).avgStim;
                    D4 = resultsData.resultsPerChan(iiC).stdStim;
                    
                    L = length(D1);
                    if locs == 6
                        avgBefore(cnt,:) = D1;
                        stdBefore(cnt,:) = D2;
                        avgStim(cnt,:) = D3;
                        stdStim(cnt,:) = D4;
                        
                    elseif locs == 7
                        avgBefore(cnt,:) = [D1(2:L),0];
                        stdBefore(cnt,:) = [D2(2:L),0];
                        avgStim(cnt,:) = [D3(2:L),0];
                        stdStim(cnt,:) = [D4(2:L),0];
                    elseif locs == 8
                        avgBefore(cnt,:) = [D1(3:L),0,0];
                        stdBefore(cnt,:) = [D2(3:L),0,0];
                        avgStim(cnt,:) = [D3(3:L),0,0];
                        stdStim(cnt,:) = [D4(3:L),0,0];
                        
                    elseif locs == 9
                        avgBefore(cnt,:) = [D1(4:L),0,0,0];
                        stdBefore(cnt,:) = [D2(4:L),0,0,0];
                        avgStim(cnt,:) = [D3(4:L),0,0,0];
                        stdStim(cnt,:) = [D4(4:L),0,0,0];
                    elseif locs == 10
                        avgBefore(cnt,:) = [D1(5:L),0,0,0,0];
                        stdBefore(cnt,:) = [D2(5:L),0,0,0,0];
                        avgStim(cnt,:) = [D3(5:L),0,0,0,0];
                        stdStim(cnt,:) = [D4(5:L),0,0,0,0];
                    elseif locs == 5
                        avgBefore(cnt,:) = [0,D1(1:L-1)];
                        stdBefore(cnt,:) = [0,D2(1:L-1)];
                        avgStim(cnt,:) = [0,D3(1:L-1)];
                        stdStim(cnt,:) = [0,D4(1:L-1)];
                    elseif locs == 4
                        avgBefore(cnt,:) = [0,0,D1(1:L-2)];
                        stdBefore(cnt,:) = [0,0,D2(1:L-2)];
                        avgStim(cnt,:) = [0,0,D3(1:L-2)];
                        stdStim(cnt,:) = [0,0,D4(1:L-2)];
                    elseif locs == 3
                        avgBefore(cnt,:) = [0,0,0,D1(1:L-3)];
                        stdBefore(cnt,:) = [0,0,0,D2(1:L-3)];   
                        avgStim(cnt,:) = [0,0,0,D3(1:L-3)];
                        stdStim(cnt,:) = [0,0,0,D4(1:L-3)];   
                    elseif locs == 2
                        avgBefore(cnt,:) = [0,0,0,0,D1(1:L-4)];
                        stdBefore(cnt,:) = [0,0,0,0,D2(1:L-4)];   
                        avgStim(cnt,:) = [0,0,0,0,D3(1:L-4)];
                        stdStim(cnt,:) = [0,0,0,0,D4(1:L-4)];   
                    else
                        disp(['maxima ', num2str(locs)])
                        disp(['ch #', num2str(cnt)])
                    end
                    
                    
                    ptNum_v = [ptNum_v, ptNum];
                    chNum_v = [chNum_v, resultsData.resultsPerChan(iiC).channelNum];
                    area = classifyArea(resultsData.resultsPerChan(iiC).area);
                    chMTL_area{cnt} = area;
                    nRip(cnt) = resultsData.resultsPerChan(iiC).nRipplesBefore;
                    nRip_stim(cnt) = resultsData.resultsPerChan(iiC).nRipplesStim;
                    
                    meanTFRRipBefore(cnt,:,:) = resultsData.resultsPerChan(iiC).meanTFRRipBefore;          
                                        
                    meanTFRRipStim(cnt,:,:) = resultsData.resultsPerChan(iiC).meanTFRRipStim;          

                end
       
            end
        end
    end
    
    
    spindleInfo.ptNum_v = ptNum_v;
    spindleInfo.chNum_v = chNum_v;
    spindleInfo.area = chMTL_area;
    spindleInfo.meanTFRRipBefore = meanTFRRipBefore;
    spindleInfo.meanTFRRipStim = meanTFRRipStim;
    spindleInfo.avgBefore = avgBefore;
    spindleInfo.stdBefore = stdBefore;
    spindleInfo.nRip = nRip;
    spindleInfo.avgStim = avgStim;
    spindleInfo.stdStim = stdStim;
    spindleInfo.nRip_stim = nRip_stim;
    
    ind = [];
    spindleInfo.artifactInd = ind; % channels where clear artifacts are contaminating average
    
    save(fullfile(dataFolder,filename),'spindleInfo')
    
end

%%
sd = SpindleDetector;

% outputFigureFolder = 'E:\Dropbox\Nir_Lab\closedLoopRevision\Figures\links_fig2';
% figSuffix = 'fullFreqRange';
% filename = fullfile(fullFreqRangeSpindlesFolder,'spindleInfo_fullFreq');
% sd.plotPopulationFig(runData, filename, figSuffix, outputFigureFolder)

figSuffix = 'highFreqRange';
filename = fullfile(HighFreqSpindlesFolder,'spindleInfo_highFreq');
outputFigureFolder = 'E:\Dropbox\Nir_Lab\closedLoopRevision\Figures\links_sup_fig10';
sd.plotPopulationFig(runData, filename, figSuffix, outputFigureFolder)

outputFigureFolder = 'E:\Dropbox\Nir_Lab\closedLoopRevision\Figures\links_fig2';
figSuffix = 'fullFreqRange_allCh';
filename = fullfile(fullFreqRangeSpindles_allChannelsFolder,'spindleInfo_fullFreq_allChanels');
aspectRatio = true;
sd.plotPopulationFig(runData, filename, figSuffix, outputFigureFolder, aspectRatio)


