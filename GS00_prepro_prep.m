%% Preparation for preprocessing
clear all 

if ispc
    PATH = ['D:\Data\grailswitch'];
elseif ismac
    PATH = ['/Users/julianreiser/owncloud/grailswitch'];
elseif isunix
    PATH = '/mnt/data/reiserOn8TB/grailswitch';
end

PATH_RAW = [PATH '/data/1_unprocessed'];
PATH_COND = [PATH '/lists/conditions'];
PATH_DFLOW = [PATH '/lists/dflow'];
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
    if strcmpi(SUBJ,'12')
        EEG = pop_loadbv(PATH_RAW, 'gs_vp12_1.vhdr', [1 2053301], [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67]);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off'); 
        EEG = pop_loadbv(PATH_RAW, 'gs_vp12_2.vhdr', [1 1350602], [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67]);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
        EEG = eeg_checkset( EEG );
        EEG = pop_mergeset( ALLEEG, [1  2], 0);
    else
        EEG = pop_loadbv(PATH_RAW,['gs_vp' SUBJ '.vhdr']);
    end

    % load GRAIL files
    PATH_dflowsub = [PATH_DFLOW '/vp' SUBJ '/'];
    
    % filter from 0.1 Hz to 100 Hz
%    EEG  = pop_basicfilter( EEG, 1:length({EEG.chanlocs.labels})-3 , 'Boundary', 'boundary', 'Cutoff',  0.1, 'Design', 'butter', 'Filter', 'highpass', 'Order',  2, 'RemoveDC', 'on' );

	% find start and end points and check for sanity (should be 15)
    b = find(strcmpi({EEG.event.type},'B  1'));
    e = find(strcmpi({EEG.event.type},'E  1'));

    % add events that are missing and delete segments that are wrong 
	if strcmpi(SUBJ,'01')
        EEG = eeg_eegrej( EEG, [1350728 1379021;1594235 1631185;2041669 2046460;3001161 3044451;3825527 3847326]);
        EEG.event(length(EEG.event) + 1).latency = EEG.event(1057).latency + 500;
        EEG.event(length(EEG.event)).type = 'E  1';
        EEG.event(length(EEG.event)).latency = EEG.event(1057).latency + 500;
        EEG.event(length(EEG.event)).duration = 1;
        EEG.event(length(EEG.event)).channel = 0;
        EEG.event(length(EEG.event)).code = 'EndMotek';
        EEG = eeg_checkset(EEG,'eventconsistency');    
    elseif strcmpi(SUBJ,'02')
        EEG = eeg_eegrej( EEG, [1403487 1485661;1672891 1729786;1925596 1995170;2281735 2291397]);
        EEG.event(1545) = [];
        EEG = eeg_checkset(EEG,'eventconsistency');    
    elseif strcmpi(SUBJ,'03')
        EEG = eeg_eegrej( EEG, [2639281 2822145] );
        EEG = eeg_checkset(EEG,'eventconsistency');    
    elseif strcmpi(SUBJ,'04')
        EEG = eeg_eegrej( EEG, [1365106 1449522;2750389 2753634]);
    elseif strcmpi(SUBJ,'07')
        EEG = eeg_eegrej( EEG, [EEG.event(b(3)).latency - 5 EEG.event(e(3)).latency + 5;...
                                EEG.event(b(12)).latency - 5 EEG.event(e(12)).latency + 5]);
    elseif strcmpi(SUBJ,'09')
        EEG = eeg_eegrej( EEG, [EEG.event(b(11)).latency - 5 EEG.event(e(11)).latency + 5;...
                               EEG.event(b(13)).latency - 5 EEG.event(e(14)).latency + 5;...
                               EEG.event(b(16)).latency - 5 EEG.event(e(16)).latency + 5;...
                               EEG.event(b(18)).latency - 5 EEG.event(e(18)).latency + 5]);
	elseif strcmpi(SUBJ,'10')
        EEG = eeg_eegrej( EEG, [EEG.event(b(4)).latency - 5 EEG.event(e(4)).latency + 5;...
                               EEG.event(b(7)).latency - 5 EEG.event(e(7)).latency + 5;...
                               EEG.event(b(10)).latency - 5 EEG.event(e(10)).latency + 5]);
    elseif strcmpi(SUBJ,'13')
        EEG = eeg_eegrej( EEG, [EEG.event(b(4)).latency - 5 EEG.event(e(4)).latency + 5;...
                               EEG.event(b(14)).latency - 5 EEG.event(e(14)).latency + 5]);
    elseif strcmpi(SUBJ,'14')
        EEG.event(end+1).type = 'E  1';
        EEG.event(end).latency = EEG.event(b(2)).latency-5;
        EEG = eeg_checkset(EEG,'eventconsistency');
        b = find(strcmpi({EEG.event.type},'B  1'));
        e = find(strcmpi({EEG.event.type},'E  1'));
        EEG.event(e([7,8])) = [];
        EEG = eeg_eegrej( EEG, [EEG.event(b(6)).latency - 5 EEG.event(e(6)).latency + 5;...
                                EEG.event(b(9)).latency - 5 EEG.event(e(9)).latency + 5]);
%         EEG = eeg_eegrej( EEG, [EEG.event(e(11)).latency - 5 EEG.event(b(12)).latency]);
        EEG = eeg_checkset(EEG,'eventconsistency');
    elseif strcmpi(SUBJ,'15')
        EEG = eeg_eegrej( EEG, [EEG.event(b(4)).latency - 5 EEG.event(e(4)).latency + 5]);
    elseif strcmpi(SUBJ,'17')
        % add end-event for standing baseline
        EEG.event(length(EEG.event) + 1).latency = EEG.event(b(1)).latency + 180*EEG.srate/1000;
        EEG.event(length(EEG.event)).type = 'E  1';
        EEG.event(length(EEG.event)).duration = 1;
        EEG.event(length(EEG.event)).channel = 0;
        EEG.event(length(EEG.event)).code = 'EndMotek';
        EEG = eeg_checkset(EEG,'eventconsistency'); 
        b = find(strcmpi({EEG.event.type},'B  1'));
        e = find(strcmpi({EEG.event.type},'E  1'));
    elseif strcmpi(SUBJ,'18')
        EEG = eeg_eegrej( EEG, [EEG.event(b(1)).latency - 5, EEG.event(e(1)).latency + 5;...
                                EEG.event(b(7)).latency - 5, EEG.event(e(10)).latency + 5]);
        EEG = eeg_checkset(EEG,'eventconsistency');
    elseif strcmpi(SUBJ,'20')
        EEG = eeg_eegrej( EEG, [EEG.event(e(6)).latency - 5, EEG.event(e(6)).latency + 5;...
                                EEG.event(b(10)).latency - 5, EEG.event(e(11)).latency + 5]);
    elseif strcmpi(SUBJ,'23')
        EEG = eeg_eegrej( EEG, [EEG.event(e(10)).latency - 5, EEG.event(e(10)).latency + 5;...
                                EEG.event(b(11)).latency - 5, EEG.event(e(12)).latency + 5]);
    elseif strcmpi(SUBJ,'24')
        EEG = eeg_eegrej( EEG, [EEG.event(e(4)).latency - 5, EEG.event(e(6)).latency + 5;...
                                EEG.event(b(10)).latency - 5, EEG.event(e(12)).latency + 5;...
                                EEG.event(b(14)).latency - 5, EEG.event(e(15)).latency + 5]);
    elseif strcmpi(SUBJ,'25')
        EEG = eeg_eegrej( EEG, [EEG.event(b(8)).latency - 5, EEG.event(e(10)).latency + 5]);
    elseif strcmpi(SUBJ,'26')
        % add end-event for walking baseline
        EEG.event(length(EEG.event) + 1).latency = EEG.event(b(2)).latency - 500;
        EEG.event(length(EEG.event)).type = 'E  1';
        EEG.event(length(EEG.event)).duration = 1;
        EEG.event(length(EEG.event)).channel = 0;
        EEG.event(length(EEG.event)).code = 'EndMotek';
        EEG = eeg_checkset(EEG,'eventconsistency'); 
        b = find(strcmpi({EEG.event.type},'B  1'));
        e = find(strcmpi({EEG.event.type},'E  1'));
        % add start-event for fifth experimental condition
        EEG.event(length(EEG.event) + 1).latency = EEG.event(e(8)).latency + 500;
        EEG.event(length(EEG.event)).type = 'B  1';
        EEG.event(length(EEG.event)).duration = 1;
        EEG.event(length(EEG.event)).channel = 0;
        EEG.event(length(EEG.event)).code = 'BeginMotek';
        EEG = eeg_checkset(EEG,'eventconsistency'); 
        b = find(strcmpi({EEG.event.type},'B  1'));
        e = find(strcmpi({EEG.event.type},'E  1'));
    elseif strcmpi(SUBJ,'27')
        EEG = eeg_eegrej( EEG, [EEG.event(b(12)).latency - 5, EEG.event(e(12)).latency + 5]);
    elseif strcmpi(SUBJ,'28')
        EEG = eeg_eegrej( EEG, [EEG.event(e(1)).latency - 5, EEG.event(e(1)).latency + 5]);
        EEG = eeg_checkset(EEG,'eventconsistency');
        
        % add end-event for walking baseline
        EEG.event(length(EEG.event) + 1).latency = EEG.event(b(1)).latency + (180*1000*EEG.srate/1000);
        EEG.event(length(EEG.event)).type = 'E  1';
        EEG.event(length(EEG.event)).duration = 1;
        EEG.event(length(EEG.event)).channel = 0;
        EEG.event(length(EEG.event)).code = 'EndMotek';
        EEG = eeg_checkset(EEG,'eventconsistency'); 
        b = find(strcmpi({EEG.event.type},'B  1'));
        e = find(strcmpi({EEG.event.type},'E  1'));
        % delete other wrong conditions
        EEG = eeg_eegrej( EEG, [EEG.event(b(6)).latency - 5, EEG.event(e(6)).latency + 5;...
                                EEG.event(b(15)).latency - 5, EEG.event(e(15)).latency + 5;...
                                EEG.event(b(16)).latency - 5, EEG.event(e(16)).latency + 5]);
        EEG = eeg_checkset(EEG,'eventconsistency');
        EEG = eeg_checkset(EEG,'eventconsistency');

    end

  	% find start and end points and check for sanity (should be 15)
    b = find(strcmpi({EEG.event.type},'B  1'))
    e = find(strcmpi({EEG.event.type},'E  1'))
    [b;e]
    e-b
    
    % exclcude data segments that are between blocks
    EEG = eeg_eegrej( EEG, [1 (EEG.event(b(1)).latency-1);...
                            (EEG.event(e(1)).latency+5) (EEG.event(b(2)).latency-5); ...
                            (EEG.event(e(2)).latency+5) (EEG.event(b(3)).latency-5); ...
                            (EEG.event(e(3)).latency+5) (EEG.event(b(4)).latency-5); ...
                            (EEG.event(e(4)).latency+5) (EEG.event(b(5)).latency-5); ...
                            (EEG.event(e(5)).latency+5) (EEG.event(b(6)).latency-5); ...
                            (EEG.event(e(6)).latency+5) (EEG.event(b(7)).latency-5); ...
                            (EEG.event(e(7)).latency+5) (EEG.event(b(8)).latency-5); ...
                            (EEG.event(e(8)).latency+5) (EEG.event(b(9)).latency-5); ...
                            (EEG.event(e(9)).latency+5) (EEG.event(b(10)).latency-5); ...
                            (EEG.event(e(10)).latency+5) (EEG.event(b(11)).latency-5); ...
                            (EEG.event(e(11)).latency+5) (EEG.event(b(12)).latency-5); ...
                            (EEG.event(e(12)).latency+5) (EEG.event(b(13)).latency-5); ...
                            (EEG.event(e(13)).latency+5) (EEG.event(b(14)).latency-5); ...
                            (EEG.event(e(14)).latency+5) (EEG.event(b(15)).latency-5); ...
                            (EEG.event(e(15)).latency+5) length(EEG.data)]);
    
  	% find start and end points and check for sanity (should be 15)
    b = find(strcmpi({EEG.event.type},'B  1'));
    e = find(strcmpi({EEG.event.type},'E  1'));
    [b;e]
    e-b
    
    % sync with dflow data
    syncedEEG = syncDataSwitch(EEG,SUBJ,PATH_DFLOW);
	syncedEEG = eeg_checkset(syncedEEG,'eventconsistency');

    % get steps
    startpoint = b(4);
    detectrange = [syncedEEG.event(startpoint).latency:syncedEEG.event(end).latency];
    stepsyncEEG = get_stepsgrail(syncedEEG,{'LeftHeel.PosY'},detectrange,'98',-1,[1],1);
    saveas(gcf,[PATH_PLOT '/stepdetect/leftfoot_subj' SUBJ '.png']);
    rightfoot = get_stepsgrail(syncedEEG,{'RightHeel.PosY'},detectrange,'99',-1,[1],1);
    saveas(gcf,[PATH_PLOT '/stepdetect/rightfoot_subj' SUBJ '.png']);
    
    % copy right toe events to stepsyncEEG (until now only left toe events)
    rightfootidx = find(strcmpi({rightfoot.event.type},'99'));
    oldevents = length(stepsyncEEG.event);
    eventcount = 1;
    for i = 1:length(rightfootidx)
        stepsyncEEG.event(oldevents+eventcount) = rightfoot.event(rightfootidx(i));
        eventcount = eventcount + 1;
    end
    stepsyncEEG = eeg_checkset(stepsyncEEG,'eventconsistency');
    
    % save the prepped set for preprocessing
    stepsyncEEG = pop_saveset(stepsyncEEG, 'filename', ['gs_vp' SUBJ '_prep.set'], 'filepath', [PATH '/data/2_raw-EEGLAB/gs_vp' num2str(subjects(subj)) '/']);

end % subj