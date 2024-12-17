%% epoching, epoch rejection, response times and ERP calculation
clear all
close all

%% PATHS, names and variables
% EEG
mac_switch = 1; % 1: work, 2: home

if ispc
   PATH = 'D:\DATA\grailswitch\';
elseif ismac
    if mac_switch == 1
        PATH = '/Volumes/Work2TB/Seafile/grailswitch';
    elseif mac_switch == 2
        PATH = '/Users/julianreiser/Seafile/grailswitch';
    end % mac_switch
end

PATH_RAW = [PATH '/data/RAW/'];
PATH_ICA = [PATH '/data/processed/ICA'];
PATH_EEG = [PATH '/data/processed/EEG'];
PATH_PRUNED =[PATH '/data/processed/PRUNED'];

% GRAIL FILES
PATH_DFLOW = [PATH '/unprocessed/dflow']
PATH_LAB = ['D:\ownCloud\projects\functions\eeglab2021.0'];
% LISTS / STATS
PATH_LIST = [PATH '/lists'];
PATH_COND = [PATH_LIST '/conditions'];
PATH_STAT = [PATH '/stats_rework'];
PATH_PLOT = [PATH '/plots'];
subjlist = {'02','03','04','05','06','07','09','10','11','12','13','15','16','18','19','20','21','22','23','24','25','27'};
% no file: 17
% no '98' step events: 28
%% initialize necessary toolboxes and variables
eeglab;
all_ufresult = [];
all_badICs = cell(length(subjlist),1);
all_winrej = cell(length(subjlist),1);
all_eventlist = cell(length(subjlist),1);

% load ERPs
load([PATH_STAT '/erp_data/allerps_cue.mat'],'allerps_cue');
load([PATH_STAT '/erp_data/allerps_tar.mat'],'allerps_tar');

% load vars
load([PATH_STAT '/ersp_data/erptime.mat'],'erp_time');

%% load subject data and combine to a large matrix
for subj = 1:length(subjlist)
    clear erp_cue erp_tar subjEventList baseidx trimat_v trimat_a

    % load eeglab file, unfold file
    VP = subjlist{subj};
    if subj == 1
        EEG = pop_loadset('filename',['gs_vp' VP '_pruned.set'],'filepath',PATH_PRUNED);
    end

    % calculate responses
    load([PATH_STAT '/beha_data/' VP '_trimat_v.mat'],'trimat_v')
    load([PATH_STAT '/beha_data/' VP '_trimat_a.mat'],'trimat_a')
    
    % create list of trials per condition
    subcond = dlmread([PATH_COND '/vp' VP '_conditions.txt'], ' ', 1, 0);
    condTrials(subj,:) = [subj,trimat_a(1,1)+trimat_a(2,1),trimat_a(5,1)+trimat_a(6,1),trimat_a(7,1)+trimat_a(8,1),trimat_v(1,1)+trimat_v(2,1),trimat_v(5,1)+trimat_v(6,1),trimat_v(7,1)+trimat_v(8,1),...
                               trimat_a(1,3)+trimat_a(2,3),trimat_a(5,3)+trimat_a(6,3),trimat_a(7,3)+trimat_a(8,3),trimat_v(1,3)+trimat_v(2,3),trimat_v(5,3)+trimat_v(6,3),trimat_v(7,3)+trimat_v(8,3),...
                               trimat_a(1,5)+trimat_a(2,5),trimat_a(5,5)+trimat_a(6,5),trimat_a(7,5)+trimat_a(8,5),trimat_v(1,5)+trimat_v(2,5),trimat_v(5,5)+trimat_v(6,5),trimat_v(7,5)+trimat_v(8,5)];

    % extract subject data
    erp_cue = squeeze(allerps_cue(subj,:,:,:));
    erp_tar = squeeze(allerps_tar(subj,:,:,:));

    % baseline idx
    baseidx = [-800,-600];

    %% extract unfold betas and baseline
    for chans = 1:size(erp_cue,2)
        
        erp_cueresult_audi(chans,:,1,subj) = squeeze(erp_cue(1,chans,:));
        erp_cueresult_audi(chans,:,2,subj) = squeeze(erp_cue(2,chans,:));
        erp_cueresult_audi(chans,:,3,subj) = squeeze(erp_cue(3,chans,:));
        erp_cueresult_audi(chans,:,4,subj) = squeeze(erp_cue(7,chans,:));
        erp_cueresult_audi(chans,:,5,subj) = squeeze(erp_cue(8,chans,:));
        erp_cueresult_audi(chans,:,6,subj) = squeeze(erp_cue(9,chans,:));
        erp_cueresult_audi(chans,:,7,subj) = squeeze(erp_cue(13,chans,:));
        erp_cueresult_audi(chans,:,8,subj) = squeeze(erp_cue(14,chans,:));
        erp_cueresult_audi(chans,:,9,subj) = squeeze(erp_cue(15,chans,:));

        erp_cueresult_visu(chans,:,1,subj) = squeeze(erp_cue(4,chans,:));
        erp_cueresult_visu(chans,:,2,subj) = squeeze(erp_cue(5,chans,:));
        erp_cueresult_visu(chans,:,3,subj) = squeeze(erp_cue(6,chans,:));
        erp_cueresult_visu(chans,:,4,subj) = squeeze(erp_cue(10,chans,:));
        erp_cueresult_visu(chans,:,5,subj) = squeeze(erp_cue(11,chans,:));
        erp_cueresult_visu(chans,:,6,subj) = squeeze(erp_cue(12,chans,:));
        erp_cueresult_visu(chans,:,7,subj) = squeeze(erp_cue(16,chans,:));
        erp_cueresult_visu(chans,:,8,subj) = squeeze(erp_cue(17,chans,:));
        erp_cueresult_visu(chans,:,9,subj) = squeeze(erp_cue(18,chans,:));

        erp_tarresult_audi(chans,:,1,subj) = squeeze(erp_tar(1,chans,:));
        erp_tarresult_audi(chans,:,2,subj) = squeeze(erp_tar(2,chans,:));
        erp_tarresult_audi(chans,:,3,subj) = squeeze(erp_tar(3,chans,:));
        erp_tarresult_audi(chans,:,4,subj) = squeeze(erp_tar(7,chans,:));
        erp_tarresult_audi(chans,:,5,subj) = squeeze(erp_tar(8,chans,:));
        erp_tarresult_audi(chans,:,6,subj) = squeeze(erp_tar(9,chans,:));
        erp_tarresult_audi(chans,:,7,subj) = squeeze(erp_tar(13,chans,:));
        erp_tarresult_audi(chans,:,8,subj) = squeeze(erp_tar(14,chans,:));
        erp_tarresult_audi(chans,:,9,subj) = squeeze(erp_tar(15,chans,:));

        erp_tarresult_visu(chans,:,1,subj) = squeeze(erp_tar(4,chans,:));
        erp_tarresult_visu(chans,:,2,subj) = squeeze(erp_tar(5,chans,:));
        erp_tarresult_visu(chans,:,3,subj) = squeeze(erp_tar(6,chans,:));
        erp_tarresult_visu(chans,:,4,subj) = squeeze(erp_tar(10,chans,:));
        erp_tarresult_visu(chans,:,5,subj) = squeeze(erp_tar(11,chans,:));
        erp_tarresult_visu(chans,:,6,subj) = squeeze(erp_tar(12,chans,:));
        erp_tarresult_visu(chans,:,7,subj) = squeeze(erp_tar(16,chans,:));
        erp_tarresult_visu(chans,:,8,subj) = squeeze(erp_tar(17,chans,:));
        erp_tarresult_visu(chans,:,9,subj) = squeeze(erp_tar(18,chans,:));

    end % for chans
    %% extract number of trials for each bin for later GND construction
    % bin_st_a_singl
    trialsPerBin_audi(subj,1) = condTrials(subj,2);
    % bin_st_a_repeat
    trialsPerBin_audi(subj,2) = condTrials(subj,3);
    % bin_st_a_mixed
    trialsPerBin_audi(subj,3) = condTrials(subj,4);
    % bin_wa_a_singl
    trialsPerBin_audi(subj,4) = condTrials(subj,8);
    % bin_wa_a_repeat
    trialsPerBin_audi(subj,5) = condTrials(subj,9);
    % bin_wa_a_mixed
    trialsPerBin_audi(subj,6) = condTrials(subj,10);
    % bin_pe_a_singl
    trialsPerBin_audi(subj,7) = condTrials(subj,14);
    % bin_pe_a_repeat
    trialsPerBin_audi(subj,8) = condTrials(subj,15);
    % bin_pe_a_mixed
    trialsPerBin_audi(subj,9) = condTrials(subj,16);
    
    % bin_st_v_singl
    trialsPerBin_visu(subj,1) = condTrials(subj,5);
    % bin_st_v_repeat
    trialsPerBin_visu(subj,2) = condTrials(subj,6);
    % bin_st_v_mixed
    trialsPerBin_visu(subj,3) = condTrials(subj,7);
    % bin_wa_v_singl
    trialsPerBin_visu(subj,4) = condTrials(subj,11);
    % bin_wa_v_repeat
    trialsPerBin_visu(subj,5) = condTrials(subj,12);
    % bin_wa_v_mixed
    trialsPerBin_visu(subj,6) = condTrials(subj,13);
    % bin_pe_v_singl
    trialsPerBin_visu(subj,7) = condTrials(subj,17);
    % bin_pe_v_repeat
    trialsPerBin_visu(subj,8) = condTrials(subj,18);
    % bin_pe_v_mixed
    trialsPerBin_visu(subj,9) = condTrials(subj,19);

end % for subj

save([PATH_STAT '/erp_data/condTrials_cue.mat'],'condTrials');

%% calculate GAcue mean
% calculate GAcues
GAcue_mean_audi = [];
GAcue_sem_audi = [];
GAcue_mean_audi = squeeze(mean(erp_cueresult_audi,4));

for bins = 1:size(erp_cueresult_audi,3)
    for chans = 1:size(erp_cueresult_audi,1)
        GAcue_sem_audi(chans,:,bins) = std(erp_cueresult_audi(chans,:,bins,:),0,4)/sqrt(size(erp_cueresult_audi(chans,:,bins,:),4));
    end % for chans
end % for bins

GAcue_t_audi = GAcue_mean_audi ./ GAcue_sem_audi;

GAcue_mean_visu = [];
GAcue_sem_visu = [];
GAcue_mean_visu = squeeze(mean(erp_cueresult_visu,4));

for bins = 1:size(erp_cueresult_visu,3)
    for chans = 1:size(erp_cueresult_visu,1)
        GAcue_sem_visu(chans,:,bins) = std(erp_cueresult_visu(chans,:,bins,:),0,4)/sqrt(size(erp_cueresult_visu(chans,:,bins,:),4));
    end % for chans
end % for bins

GAcue_t_visu = GAcue_mean_visu ./ GAcue_sem_visu;

%% build GND struct for further FMUT processing
% auditory
binlistGND = {'st_singl','st_repeat','st_switch',...
              'wa_singl','wa_repeat','wa_switch',...
              'pe_singl','pe_repeat','pe_switch'};

GND = [];
GND.exp_desc = ['Cued Switch Task in the GRAIL - auditory - cue locked'];
GND.filename = ['grailswitch_GND_cue_audi'];
GND.filepath = [PATH_STAT '/FMUT/'];
GND.saved = 'no';
GND.grands = [GAcue_mean_audi];
GND.grands_stder = [GAcue_sem_audi];
GND.grands_t = [GAcue_t_audi];
GND.sub_ct = [size(erp_cueresult_audi,4),size(erp_cueresult_audi,4),size(erp_cueresult_audi,4),size(erp_cueresult_audi,4),...
              size(erp_cueresult_audi,4),size(erp_cueresult_audi,4),size(erp_cueresult_audi,4),size(erp_cueresult_audi,4),...
              size(erp_cueresult_audi,4),size(erp_cueresult_audi,4),size(erp_cueresult_audi,4),size(erp_cueresult_audi,4)];
GND.chanlocs = EEG.chanlocs;
for bini = 1:length(binlistGND)
    GND.bin_info(bini).bindesc = binlistGND{bini};
    GND.bin_info(bini).condcode = 1;
end
GND.condesc = binlistGND;
GND.time_pts = erp_time;
GND.bsln_wind = [-800,-600];
GND.odelay = [];
GND.srate = [EEG.srate];
GND.indiv_fnames = [];
GND.indiv_subnames = subjlist;
GND.indiv_traits = [];
GND.indiv_bin_ct = trialsPerBin_audi;
GND.indiv_bin_raw_ct = [];
GND.indiv_erps = [erp_cueresult_audi];
GND.indiv_art_ics = [all_badICs];
GND.cals = [];
GND.history = [];
GND.t_tests = [];

save([PATH_STAT '/FMUT/GND_grailswitch_cue_audi.mat'],'GND')

% visual
GND = [];
GND.exp_desc = ['Cued Switch Task in the GRAIL - visual - cue locked'];
GND.filename = ['grailswitch_GND_cue_visu'];
GND.filepath = [PATH_STAT '/FMUT/'];
GND.saved = 'no';
GND.grands = [GAcue_mean_visu];
GND.grands_stder = [GAcue_sem_visu];
GND.grands_t = [GAcue_t_visu];
GND.sub_ct = [size(erp_cueresult_visu,4),size(erp_cueresult_visu,4),size(erp_cueresult_visu,4),size(erp_cueresult_visu,4),...
              size(erp_cueresult_visu,4),size(erp_cueresult_visu,4),size(erp_cueresult_visu,4),size(erp_cueresult_visu,4),...
              size(erp_cueresult_visu,4),size(erp_cueresult_visu,4),size(erp_cueresult_visu,4),size(erp_cueresult_visu,4)];
GND.chanlocs = EEG.chanlocs;
for bini = 1:length(binlistGND)
    GND.bin_info(bini).bindesc = binlistGND{bini};
    GND.bin_info(bini).condcode = 1;
end
GND.condesc = binlistGND;
GND.time_pts = erp_time;
GND.bsln_wind = [-800,-600];
GND.odelay = [];
GND.srate = [EEG.srate];
GND.indiv_fnames = [];
GND.indiv_subnames = subjlist;
GND.indiv_traits = [];
GND.indiv_bin_ct = trialsPerBin_visu;
GND.indiv_bin_raw_ct = [];
GND.indiv_erps = [erp_cueresult_visu];
GND.indiv_art_ics = [all_badICs];
GND.cals = [];
GND.history = [];
GND.t_tests = [];

save([PATH_STAT '/FMUT/GND_grailswitch_cue_visu.mat'],'GND') 

%% calculate GAtar mean
% calculate GAtars
GAtar_mean_audi = [];
GAtar_sem_audi = [];
GAtar_mean_audi = squeeze(mean(erp_tarresult_audi,4));

for bins = 1:size(erp_tarresult_audi,3)
    for chans = 1:size(erp_tarresult_audi,1)
        GAtar_sem_audi(chans,:,bins) = std(erp_tarresult_audi(chans,:,bins,:),0,4)/sqrt(size(erp_tarresult_audi(chans,:,bins,:),4));
    end % for chans
end % for bins

GAtar_t_audi = GAtar_mean_audi ./ GAtar_sem_audi;

GAtar_mean_visu = [];
GAtar_sem_visu = [];
GAtar_mean_visu = squeeze(mean(erp_tarresult_visu,4));

for bins = 1:size(erp_tarresult_visu,3)
    for chans = 1:size(erp_tarresult_visu,1)
        GAtar_sem_visu(chans,:,bins) = std(erp_tarresult_visu(chans,:,bins,:),0,4)/sqrt(size(erp_tarresult_visu(chans,:,bins,:),4));
    end % for chans
end % for bins

GAtar_t_visu = GAtar_mean_visu ./ GAtar_sem_visu;

%% build GND struct for further FMUT processing
% auditory
binlistGND = {'st_singl','st_repeat','st_switch',...
              'wa_singl','wa_repeat','wa_switch',...
              'pe_singl','pe_repeat','pe_switch'};

GND = [];
GND.exp_desc = ['tard Switch Task in the GRAIL - auditory - tar locked'];
GND.filename = ['grailswitch_GND_tar_audi'];
GND.filepath = [PATH_STAT '/FMUT/'];
GND.saved = 'no';
GND.grands = [GAtar_mean_audi];
GND.grands_stder = [GAtar_sem_audi];
GND.grands_t = [GAtar_t_audi];
GND.sub_ct = [size(erp_tarresult_audi,4),size(erp_tarresult_audi,4),size(erp_tarresult_audi,4),size(erp_tarresult_audi,4),...
              size(erp_tarresult_audi,4),size(erp_tarresult_audi,4),size(erp_tarresult_audi,4),size(erp_tarresult_audi,4),...
              size(erp_tarresult_audi,4),size(erp_tarresult_audi,4),size(erp_tarresult_audi,4),size(erp_tarresult_audi,4)];
GND.chanlocs = EEG.chanlocs;
for bini = 1:length(binlistGND)
    GND.bin_info(bini).bindesc = binlistGND{bini};
    GND.bin_info(bini).condcode = 1;
end
GND.condesc = binlistGND;
GND.time_pts = erp_time;
GND.bsln_wind = [-200,0];
GND.odelay = [];
GND.srate = [EEG.srate];
GND.indiv_fnames = [];
GND.indiv_subnames = subjlist;
GND.indiv_traits = [];
GND.indiv_bin_ct = trialsPerBin_audi;
GND.indiv_bin_raw_ct = [];
GND.indiv_erps = [erp_tarresult_audi];
GND.indiv_art_ics = [all_badICs];
GND.cals = [];
GND.history = [];
GND.t_tests = [];

save([PATH_STAT '/FMUT/GND_grailswitch_tar_audi.mat'],'GND')

% visual
GND = [];
GND.exp_desc = ['tard Switch Task in the GRAIL - visual - tar locked'];
GND.filename = ['grailswitch_GND_tar_visu'];
GND.filepath = [PATH_STAT '/FMUT/'];
GND.saved = 'no';
GND.grands = [GAtar_mean_visu];
GND.grands_stder = [GAtar_sem_visu];
GND.grands_t = [GAtar_t_visu];
GND.sub_ct = [size(erp_tarresult_visu,4),size(erp_tarresult_visu,4),size(erp_tarresult_visu,4),size(erp_tarresult_visu,4),...
              size(erp_tarresult_visu,4),size(erp_tarresult_visu,4),size(erp_tarresult_visu,4),size(erp_tarresult_visu,4),...
              size(erp_tarresult_visu,4),size(erp_tarresult_visu,4),size(erp_tarresult_visu,4),size(erp_tarresult_visu,4)];
GND.chanlocs = EEG.chanlocs;
for bini = 1:length(binlistGND)
    GND.bin_info(bini).bindesc = binlistGND{bini};
    GND.bin_info(bini).condcode = 1;
end
GND.condesc = binlistGND;
GND.time_pts = erp_time;
GND.bsln_wind = [-200,0];
GND.odelay = [];
GND.srate = [EEG.srate];
GND.indiv_fnames = [];
GND.indiv_subnames = subjlist;
GND.indiv_traits = [];
GND.indiv_bin_ct = trialsPerBin_visu;
GND.indiv_bin_raw_ct = [];
GND.indiv_erps = [erp_tarresult_visu];
GND.indiv_art_ics = [all_badICs];
GND.cals = [];
GND.history = [];
GND.t_tests = [];

save([PATH_STAT '/FMUT/GND_grailswitch_tar_visu.mat'],'GND')