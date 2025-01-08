%% Preparation for preprocessing
clear all

if ispc == 1
    PATH = 'D:/Data/grailswitch';
    PATH_LAB = ['D:\ownCloud\projects\functions\eeglab2021.0'];
elseif ismac == 1
    PATH = '/Users/julianreiser/owncloud/grailswitch/';
elseif isunix == 1
    PATH = '/mnt/data/reiserOn8TB/grailswitch';
    PATH_LAB = '/mnt/data/reiserOn8TB/projects/functions/eeglab2021.0';
end

% EEG
PATH_RAW = [PATH '/data/RAW/'];
PATH_ICA = [PATH '/data/processed/ICA'];
PATH_EEG = [PATH '/data/processed/EEG'];
PATH_PRUNED =[PATH '/data/processed/PRUNED'];
PATH_STUDY = [PATH '/data/processed/STUDY'];
PATH_PLOT = [PATH '/plots/'];

% GRAIL FILES
PATH_DFLOW = [PATH '/unprocessed/dflow']

% LISTS / STATS
PATH_LIST = [PATH '/lists/'];
PATH_COND = [PATH_LIST '/conditions'];
PATH_STAT = [PATH_LIST '/stats/'];
subjlist = {'02','03','04','05','06','07','09','10','11','12','13','15','16','17','18','19','20','21','22','23','24','25','27','28','29'};
allsubs = {'02','03','04','05','06','07','09','10','11','12','13','15','16','17','18','19','20','21','22','23','24','25','27','28','29'};
% 08 no markers...
% 14 with break in one condition


for subj = 1:length(subjlist)
    eeglab;
    
    SUBJ = subjlist{subj};
    EEG = pop_loadset('filename',['gs_vp' SUBJ '_prep.set'],'filepath',PATH_RAW);

    % save channel locations
    EEG = pop_chanedit(EEG, 'lookup', [PATH_LAB '/plugins/dipfit3.7/standard_BESA/standard-10-5-cap385.elp']);
    channel_locations = EEG.chanlocs;
    OLD = EEG;
    EEG = pop_select(EEG,'nochannel',{EEG.chanlocs(65:end).labels});

    %% define basics
    
    %% use zapline
    % noise frequencies are at 50 Hz (line), 21 Hz (probably PhaseSpace motion capture system), 89.35 Hz (VR refresh rate)
    zapline_config = struct('linefreqs',[]);
    zapline_config.initialsigma = 2.8; 
    EEG_clean = bemobil_clean_data_with_zapline(EEG,zapline_config);
    saveas(gcf,[PATH_PLOT '/line/linenoise_subj' SUBJ '.png']);

    % plot cleaned
    pop_eegplot( EEG_clean, 1, 1, 1);
    linescoreplot = gcf;
    saveas(gcf,[PATH_PLOT '/line/sub' SUBJ '.png']);
    %% channel correction / rejection
    chancorr_crit = 0.8;
    chan_max_broken_time = 0.3;
    chan_detect_num_iter = 10;
    chan_detected_fraction_threshold = 0.5;
    flatline_crit = 'off';
    line_noise_crit = 'off';

    chans_to_interp = bemobil_detect_bad_channels(EEG_clean, ALLEEG, CURRENTSET,chancorr_crit,chan_max_broken_time,...
        chan_detect_num_iter,chan_detected_fraction_threshold,flatline_crit,line_noise_crit);
    
    save([PATH_STAT '/analysis/' SUBJ '_chans2interp.mat'],'chans_to_interp');


    % do the actual interpolation and full rank average referencing (no rank reduction)
    [ALLEEG, EEG_preprocessed, CURRENTSET] = bemobil_interp_avref( EEG_clean , ALLEEG, CURRENTSET, chans_to_interp);
    
    %% filter for AMICA
    % See Klug & Gramann (2020) for an investigation of filter effect on AMICA -> 1.25 Hz should be a good compromise if you
    % don't know how much movement exists, otherwise even higher may be good, up to 2Hz, and you need to subtract 0.25 to
    % obtain the correct cutoff value for a filter order of 1650
    filter_lowCutoffFreqAMICA = 1.5; % 1.75 is 1.5Hz cutoff!
    filter_AMICA_highPassOrder = 1650; % was used by Klug & Gramann (2020)
    filter_highCutoffFreqAMICA = []; % not used
    filter_AMICA_lowPassOrder = []; 

    out_filename = [];
    out_filepath= [];

    [ ALLEEG, EEG_filtered, CURRENTSET ] = bemobil_filter(ALLEEG, EEG_preprocessed, CURRENTSET, filter_lowCutoffFreqAMICA,...
        filter_highCutoffFreqAMICA, out_filename, out_filepath, filter_AMICA_highPassOrder, filter_AMICA_lowPassOrder);
    
    %% compute AMICA
    % data rank is the number of channels that were not interpolated
    data_rank = EEG_filtered.nbchan - length(EEG_filtered.etc.interpolated_channels);

    % additional AMICA settings
    amica = true;
    numb_models = 1; % default 1
    AMICA_autoreject = 1; % uses automatic rejection method of AMICA. no time-cleaning (manual or automatic) is needed then!
    AMICA_n_rej = 10; % number of iterations during automatic rejection, default 10
    AMICA_reject_sigma_threshold = 3; % threshold for rejection, default 3

    % reduced number of iterations for demo purpose only!
    AMICA_max_iter = 2000; % maximum number of AMICA iterations, default 2000

    % 4 threads are most effective for single subject speed, more threads don't really shorten the calculation time much.
    % best efficiency is using just 1 thread and have as many matlab instances open as possible (limited by the CPU usage).
    % Remember your RAM limit in this case.
    max_threads = 4;

    other_algorithm = [];
    out_filename = [];
    out_filepath = [];

    % compute AMICA
    [ALLEEG, EEG_amica, CURRENTSET] = bemobil_signal_decomposition(ALLEEG, EEG_filtered, CURRENTSET,...
        amica, numb_models, max_threads, data_rank, other_algorithm, out_filename, out_filepath, AMICA_autoreject,...
        AMICA_n_rej, AMICA_reject_sigma_threshold,AMICA_max_iter);
    
    %% Plot AMICA results
    % plot topographies
    pop_topoplot(EEG_amica, 0, [1:49] ,'final_processed',[7 7] ,0,'electrodes','off');

    % plot autorejection
    data2plot = EEG_amica.data(1:round(EEG_amica.nbchan/10):EEG_amica.nbchan,:)';
    figure;
    set(gcf,'color','w','Position', get(0,'screensize'));
    plot(data2plot,'g');
    data2plot(~EEG_amica.etc.bad_samples,:) = NaN;
    hold on
    plot(data2plot,'r');
    xlim([-10000 EEG_amica.pnts+10000])
    ylim([-1000 1000])
    title(['AMICA autorejection, removed ' num2str(round(EEG_amica.etc.bad_samples_percent,2)) '% of the samples'])
    xlabel('Samples')
    ylabel('\muV')
    saveas(gcf,[PATH_PLOT '/line/amicaresults_subj' SUBJ '.png']);

    
    %% fit dipoles
    RV_threshold = 100; % no ICs are rejected
    remove_outside_head = 'off'; % no ICs are rejected
    number_of_dipoles = 1; % it is possible to fit dual dipoles but usually not recommended

    [ALLEEG, EEG_dipfit, CURRENTSET] = bemobil_dipfit( EEG_amica , ALLEEG, CURRENTSET, [], RV_threshold,...
        remove_outside_head, number_of_dipoles);
    
    %% Copy the spatial filter data into the unfiltered data set for further processing
    % e.g. ICLabel should be done on unfiltered data (or less than 1Hz highpass), and ERPs should have a lower filter than
    % the one used for ICA usually
    disp('Copying all information into full length dataset for single subject processing...');
    [ALLEEG, EEG_AMICA_copied, CURRENTSET] = bemobil_copy_spatial_filter(EEG_preprocessed, ALLEEG, CURRENTSET,...
        EEG_dipfit);
    
    EEG_AMICA_copied = pop_iclabel(EEG_AMICA_copied, 'default');
    EEG_AMICA_copied = DO_Blink_detect_WO(EEG_AMICA_copied,.9,.75);
    EEG_AMICA_copied = DO_CodeBlinksWO_gs(EEG_AMICA_copied);
    EEG_AMICA_copied = eeg_checkset(EEG_AMICA_copied,'eventconsistency')
    %% clean with IClabel

    % 'default' classifier is slightly better for brain ICs but did not lead to good classification of muscles (see Klug &
    % Gramann (2020)), 'lite' was better overall.
    iclabel_classifier = 'default';
    % 'Brain', 'Muscle', 'Eye', 'Heart', 'Line Noise', 'Channel Noise', 'Other'
    iclabel_classes = 1; 
    % if the threshold is set to -1, the popularity classifier is used (i.e. every IC gets the class with the highest
    % probability), if it is set to a value, the summed score of the iclabel_classes must be higher than this threshold to
    % keep an IC. Must be in the [0 1] range!
    iclabel_threshold = 0.5;

    [ALLEEG, EEG_preprocessed_and_ICA, CURRENTSET, ICs_keep, ICs_throw] = bemobil_clean_with_iclabel( EEG_AMICA_copied ,...
        ALLEEG, CURRENTSET, iclabel_classifier,...
        iclabel_classes, iclabel_threshold);
    
    % plot final ICA cleaned
    saveas(gcf,[PATH_PLOT '/line/iclabel_subj' SUBJ '.png']);

    % pop_eegplot( EEG_preprocessed_and_ICA, 1, 1, 1);

	EEG_AMICA_copied = pop_saveset(EEG_AMICA_copied, 'filename', ['gs_vp' SUBJ '_amicacopied.set'], 'filepath', [PATH_ICA]);
    EEG_preprocessed_and_ICA = pop_saveset(EEG_preprocessed_and_ICA, 'filename', ['gs_vp' SUBJ '_pruned.set'], 'filepath', [PATH_PRUNED]);

end % subj
