%% Prepare data for gait analysis
clear all 

if ispc
    PATH = ['D:\Data\grailswitch'];
elseif ismac
    PATH = ['/Users/julianreiser/ownCloud/grailswitch'];
elseif isunix
    PATH = '/mnt/data/reiserOn8TB/grailswitch';
end

PATH_RAW = [PATH '/data/2_raw-EEGLAB'];
PATH_COND = [PATH '/lists/conditions'];
PATH_DFLOW = [PATH '/lists/dflow'];
PATH_STAT = [PATH '/lists/stats'];
PATH_PLOT = [PATH '/plots'];
subjlist = {'02','03','04','05','06','07','09','10','11','12','13','15','16','17','18','19','20','21','22','23','24','25','27','28','29'};
subjects = [2,3,4,5,6,7,9,10,11,12,13,15,16,17,18,19,20,21,22,23,24,25,27,28,29];

movelist = {'ST','WA','PE'};
modalist = {'V','A'};
tasklist = {'G','R','MIX'};


for subj = 1:length(subjlist)
    %%
    close all
    eeglab
    SUBJ = subjlist{subj};
    EEG = pop_loadset(['gs_vp' SUBJ '_prep.set'],[PATH_RAW '/gs_vp' SUBJ]);

    %% get rid of GRAIL datachannels and steppin events
    % delete chans
    EEG = pop_select(EEG,'rmchannel',[68:length(EEG.chanlocs)]);
    % find stepping events
    stepidx = find(ismember({EEG.event.type},{'98','99'}));
    EEG.event(stepidx) = [];
    EEG = eeg_checkset(EEG,'eventconsistency');

    %% load GRAIL files
    PATH_dflowsub = [PATH_DFLOW '/vp' SUBJ '/'];

	% find start and end points and check for sanity (should be 15)
    b = find(strcmpi({EEG.event.type},'B  1'));
    e = find(strcmpi({EEG.event.type},'E  1'));

    % add event information
    boundIdx = find(strcmpi({EEG.event.type},'boundary'));
    if size([b;e],2) == 11
        EEG.event(end+1).type = 'E  1';
        EEG.event(end).latency = EEG.event(end-1).latency+2;
        EEG.event(end).code = 'EndMotek';
        EEG.event(end+1).type = 'B  1';
        EEG.event(end).latency = ceil(EEG.event(boundIdx(1)).latency+1);
        EEG.event(end).code = 'BeginMotek';
    end % if size
    EEG = eeg_checkset(EEG,'eventconsistency');

  	% find start and end points and check for sanity (should be 15)
    b = find(strcmpi({EEG.event.type},'B  1'))
    e = find(strcmpi({EEG.event.type},'E  1'))
    [b;e]
    e-b
    
    %% now sync with dflow data
    subfolder = ['subj' SUBJ];
    
    % Read and clean EEG data
    evt = {EEG.event.code}';
    % evtFrames = vertcat(eegData.event.latency);

    % Match each begin trigger with an end trigger and remove extras
    blockMarkerIdx = find(strcmp(evt,'BeginMotek') | strcmp(evt,'EndMotek'));
    removeIdx = [];
    expectedMarker = 'BeginMotek';
    for b = 1:length(blockMarkerIdx)
        if ~strcmp(evt(blockMarkerIdx(b)),expectedMarker)
            if strcmp(evt(blockMarkerIdx(b)),'BeginMotek')
                removeIdx = [removeIdx; b-1];
            else
                removeIdx = [removeIdx; b];
            end
        else
            if strcmp(evt(blockMarkerIdx(b)),'BeginMotek')
                expectedMarker = 'EndMotek';
            else
                expectedMarker = 'BeginMotek';
            end
        end
    end
    blockMarkerIdx(removeIdx) = [];
    fullBlockIdx = (1:12)';
    blockMarkerIdx = blockMarkerIdx([(fullBlockIdx.*2)-1 fullBlockIdx.*2]);

    %% Organize GRAIL files
    % Get seperate lists for forceplate and trigger files
    fpFiles = dir([fullfile(PATH_DFLOW, subfolder, ['*forceplate*.txt'])]);
    triggerFiles = dir(fullfile(PATH_DFLOW, subfolder, ['*treadmill*.txt']));
    % Find latest recording for each condition (files for each condition are
    % in chronological order because of the timestamp in the name)

    % Sort collected files by timestamp
    [~,tmpIdx] = sort(string(regexp({fpFiles.name},'\d+-\d+-\d+_\d+','match')));
    fpFiles = fpFiles(tmpIdx);
    [~,tmpIdx] = sort(string(regexp({triggerFiles.name},'\d+-\d+-\d+_\d+','match')));
    triggerFiles = triggerFiles(tmpIdx);

    % exclude baseline / training blocks
    fpBaIdxTmp = strfind({fpFiles.name},'BA');
    fpBaIdx = find(~cellfun(@isempty,fpBaIdxTmp));
    fpFiles(fpBaIdx) = [];

    trigBaIdxTmp = strfind({triggerFiles.name},'BA');
    trigBaIdx = find(~cellfun(@isempty,trigBaIdxTmp));
    triggerFiles(trigBaIdx) = [];

    syncedEEG = EEG;
    %% Add force data to EEG
    markerChanNames = {};
    fpChanNames = {};
    triggerNames = {};
    currentBlock = 1;
    lastRowIdx = size(syncedEEG.data,1);
    for i = 1:length(fpFiles)
        disp(['starting block ' num2str(i)])
        currentEEGBlockMarkerIdx = blockMarkerIdx(currentBlock,:);
        eegStartIdx = syncedEEG.event(currentEEGBlockMarkerIdx(1)).latency;
        eegEndIdx = syncedEEG.event(currentEEGBlockMarkerIdx(2)).latency;

        eegTime = syncedEEG.times(ceil(eegStartIdx):eegEndIdx) - syncedEEG.times(eegStartIdx);

        fid = fopen(fullfile(PATH_DFLOW,subfolder,fpFiles(i).name),'r');
        tid = fopen(fullfile(PATH_DFLOW,subfolder,triggerFiles(i).name),'r');
        if fid == -1
            fprintf('\nThe source file ''%s'' was not found.\n',currentFile);
            continue
        else
            fpHeader = textscan(fid,'%s',1,'Delimiter','\r\n', 'EndOfLine','\r\n');
            fclose(fid);
            fpHeader = fpHeader{1};
            fpHeader = fpHeader{1};

            fpHeader = regexp(fpHeader,'\t','split');

            triggerHeader = textscan(tid,'%s',1,'Delimiter','\r\n', 'EndOfLine','\r\n');
            fclose(tid);
            triggerHeader = triggerHeader{1};
            triggerHeader = triggerHeader{1};

            triggerHeader = regexp(triggerHeader,'\t','split');
        end

        %% now write mocap data and forceplate data
        fpData = dlmread(fullfile(PATH_DFLOW,subfolder,fpFiles(i).name),'\t',1,0);
        fpTime = fpData(:,1) - fpData(1,1);
        % find markers of interest
        idx = ~cellfun('isempty',regexp(fpHeader,'((Right)|(Left))((Heel)|(Toe)|(Ankle))|(Head)|(Chest)|(Hip)|(Leg)|(Knee)|(Arm)|(Shoulder)|(Wrist)|(Elbow)'));
        %exclude rotational markers
        idx2 = ~contains(fpHeader,'Rot');
        % exclude VICON precomputed (empty) datachannels:
        % Hip: 81:83
        % Head: 129:131
        idx3 = ones(1,length(fpHeader)); idx3([81:83,129:131]) = 0;
        
        % select correct channels
        fpHeaderRedux = fpHeader(idx & idx2 & idx3);
        markerData = fpData(:,idx & idx2 & idx3);
    
        % create empty data arrays
        if isempty(markerChanNames)
            % mocapDataBeg = find(strcmpi(fpHeader,'CentralHip.PosX'));
            % mocapDataEnd = find(strcmpi(fpHeader,'LeftLowerArm.RotZ'));
            markerChanNames = fpHeader(idx & idx2 & idx3)';
            syncedEEG.data = [syncedEEG.data; zeros(size(markerData(:,length(markerChanNames)),2),size(syncedEEG.data,2))];
        end

        % only interpolate real marker Data, not already aggregated ones (hip, upper body, etc.)
        clear mocapData
        mocapData = markerData(:, :);
 
        % first cleaning for 0-values, then cleaning for outliers
        nullIdx = zeros(size(mocapData));
        nullIdx2 = zeros(size(mocapData));

        % go through every marker (xyz-coordinates, so 3 channels for 1 marker)
        for chani = 1:3:size(mocapData,2)-2

            % go through all time-points and look at flat-lines
            for nulli = 1:size(mocapData,1)
                
                % 1) look for periods of 0 values in xyz channels, set them to nan
                if all(mocapData(nulli,chani:chani + 2)) == 0
                    nullIdx(nulli,chani:chani+2) = 1;
                end
            end %for nulli
          
            % before outlier-detection, set all flatline timepoints for marker chani to nan (so it outlier-identification is not misguided by zero-values
            mocapDataOutl = mocapData;
            mocapDataOutl(find(nullIdx == 1)) = nan;

            % define a time-point as a nan-value, if any of the xyz-data-points of a channel is defined as an outlier using a moving average and a MAD-criterion of 20
            nullIdx2(:,chani:chani+2) = isoutlier(mocapDataOutl(:,chani:chani+2),'movmedian',size(mocapDataOutl,1)/100,1,'ThresholdFactor',100);

            for nulli = 1:size(nullIdx2,1)
                
                % 2) if there is an outlier in any of the xyz-channels, set them to nan
                if any(nullIdx2(nulli,chani:chani+2) == 1)
                    nullIdx2(nulli,chani:chani+2) = 1;
                end % if any
            end % for nulli
            
            % optional plotting for debugging
            % plot(1:size(mocapData,1),mocapDataOutl(:,chani)); hold on; plot(1:size(mocapData,1),mocapDataOutl(:,chani+1)); plot(1:size(mocapData,1),mocapDataOutl(:,chani+2)); plot(1:size(mocapData,1),nullIdx(:,chani)); plot(1:size(mocapData,1),nullIdx2(:,chani)-1); hold off;
            % close gcf
        end
        
        % exchange outlier values for nan
        mocapDataOutl(find(nullIdx2 == 1)) = nan;

        % use data-driven interpolation method from Gloersen & Federolf (2016)
        markerDataInterp = PredictMissingMarkers(mocapDataOutl,'Algorithm',2);

        % optional for debugging: plot interpolation results
        figure('Visible','off');
        for chani = 1:3:size(markerDataInterp,2)
            plot(markerDataInterp(:,chani)-1)
            hold on;
            plot(markerData(:,chani))
        end
        hold off;
        set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 50, 35], 'PaperUnits', 'centimeters', 'PaperSize', [29.7,21])
        saveas(gcf,[PATH_PLOT '/vicon/subj' SUBJ '/gs_' SUBJ '_block' num2str(i) '.png']); close gcf;

        %% add data to EEG dataset
        for m = 1:size(markerDataInterp,2)
            markerDataRsmp = interp1(fpTime.*1000,markerDataInterp(:,m),eegTime,'nearest');
            syncedEEG.data(lastRowIdx+m,eegStartIdx:eegEndIdx) = markerDataRsmp';
        end

        idx = ~cellfun('isempty',regexp(fpHeader,'FP[12].Cop[XZ]'));
        copData = fpData(:,idx);
        if isempty(fpChanNames)
            fpChanNames = fpHeader(idx);
            syncedEEG.data = [syncedEEG.data; zeros(size(copData,2),size(syncedEEG.data,2))];
        end

        for c = 1:size(copData,2)
           copDataRsmp = interp1(fpTime.*1000,copData(:,c),eegTime,'spline');
           syncedEEG.data(lastRowIdx+size(markerDataInterp,2)+c,eegStartIdx:eegEndIdx) = copDataRsmp';
        end

        % get rid off all the mocap data to save some RAM
        clear markerData mocapData mocapDataOutl nullIdx nullIdx2

        %% now write vgait pitch data
        triggerData = dlmread(fullfile(PATH_DFLOW,subfolder,triggerFiles(i).name),'\t',1,0);
        triggerTime = triggerData(:,1) - triggerData(1,1);
        idx = ~cellfun('isempty',regexp(triggerHeader,'(VGait_Pitch)'));
        pitchData = triggerData(:,idx);

        if isempty(triggerNames)
            triggerNames = triggerHeader(idx);
            syncedEEG.data = [syncedEEG.data; zeros(size(pitchData,2),size(syncedEEG.data,2))];
        end

        for t = 1:size(pitchData,2)
            pitchDataRsmp = interp1(triggerTime.*1000,pitchData(:,t),eegTime,'spline');
            syncedEEG.data(lastRowIdx+size(markerDataInterp,2)+ size(copData,2) + t,eegStartIdx:eegEndIdx) = pitchDataRsmp';
        end

        currentBlock = currentBlock + 1;
    end % for i = fpFiles

    for i = 1:length(markerChanNames)
        syncedEEG.chanlocs(end+1).labels = markerChanNames{i};
    end

    for i = 1:length(fpChanNames)
        syncedEEG.chanlocs(end+1).labels = fpChanNames{i};
    end

    for i = 1:length(triggerNames)
        syncedEEG.chanlocs(end+1).labels = triggerNames{i};
    end

    clear pitchDataRsmp fpTime fpData pitchData pitchDataRsmp  markerDataRsmp markerDataInterp triggerData triggerTime copData copDataRsmp eegTime
   
    %% write steps
    % look again for start and end events
    b = find(strcmpi({syncedEEG.event.type},'B  1'));
    e = find(strcmpi({syncedEEG.event.type},'E  1'));

    % now calculate the steps :)
    startpoint = b(1);
    detectrange = [syncedEEG.event(startpoint).latency:syncedEEG.event(end).latency];
    [stepsyncEEG boundtmpl] = get_stepslaviegrail(syncedEEG,{'LeftHeel.PosY'},detectrange,'98',-1,[],1);
    set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 50, 35], 'PaperUnits', 'centimeters', 'PaperSize', [29.7,21])
    saveas(gcf,[PATH_PLOT '/stepdetect/leftfoot_gs_' SUBJ '_new.png']); close gcf;
    [rightfoot boundtmpr] = get_stepslaviegrail(syncedEEG,{'RightHeel.PosY'},detectrange,'99',-1,[],1);
    set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 50, 35], 'PaperUnits', 'centimeters', 'PaperSize', [29.7,21])
    saveas(gcf,[PATH_PLOT '/stepdetect/rightfoot_gs_' SUBJ '_new.png']); close gcf;

    % copy right toe events to stepsyncEEG (until now only left toe events)
    rightfootidx = find(strcmpi({rightfoot.event.type},'99'));
    oldevents = length(stepsyncEEG.event);
    eventcount = 1;
    for i = 1:length(rightfootidx)
        stepsyncEEG.event(oldevents+eventcount) = rightfoot.event(rightfootidx(i));
        eventcount = eventcount + 1;
    end
    stepsyncEEG = eeg_checkset(stepsyncEEG,'eventconsistency');
    
    % save stepmats
    stepmatl = sortrows(boundtmpl,[1,2,3]);
    stepmatr = sortrows(boundtmpl,[1,2,3]);
    save(['/Volumes/Work2TB/Seafile/grailswitch/gs_vp' SUBJ '_stepmatl.mat'],'stepmatl');
    save(['/Volumes/Work2TB/Seafile/grailswitch/gs_vp' SUBJ '_stepmatr.mat'],'stepmatr');

    % save the prepped set for preprocessing
    stepsyncEEG = pop_saveset(stepsyncEEG, 'filename', ['gs_vp' SUBJ '_corrvicon.set'], 'filepath', ['/Volumes/Work2TB/Seafile/grailswitch/']);

end % subj