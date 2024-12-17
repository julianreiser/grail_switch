 %% epoching, epoch rejection, response times and ERP calculation
clear all

%% PATHS, names and variables
mac_switch = 1; % 1: work, 2: home
% EEG
if ispc
    PATH = 'D:\data\grailswitch\';
elseif ismac
    if mac_switch == 1
        PATH = '/Users/julianreiser/ownCloud/grailswitch/';
    else
    end % if mac_switch
end % if ispc


PATH_RAW = [PATH '/data/RAW/'];
PATH_ICA = [PATH '/data/processed/ICA'];
PATH_EEG = [PATH '/data/processed/EEG'];
PATH_PRUNED =[PATH '/data/processed/PRUNED'];

% GRAIL FILES
PATH_DFLOW = [PATH '/unprocessed/dflow']
PATH_LAB = ['D:\ownCloud\projects\functions\eeglab2021.0'];
% LISTS / STATS
PATH_LIST = [PATH '/lists/'];
PATH_COND = [PATH_LIST '/conditions'];
PATH_STAT = ['/Volumes/Work2TB/Seafile/grailswitch/stats_rework'];
subjlist = {'27'};
allsubjlist = {'02','03','04','05','06','07','09','10','11','12','13','15','16','18','19','20','21','22','23','24','25','27'};

for subj = 1:length(subjlist)

    eeglab;
    % load subject eeg data
    SUBJ = subjlist{subj};
    EEG = pop_loadset('filename',['gs_vp' SUBJ '_pruned.set'],'filepath',PATH_PRUNED);
    subcond = dlmread([PATH_COND '/vp' SUBJ '_conditions.txt'], ' ', 1, 0);

    % filter data
    EEG = pop_basicfilter(EEG, 1:EEG.nbchan , 'Boundary', 'boundary', 'Cutoff', [ 0.1 35], 'Design', 'butter', 'Filter', 'bandpass', 'Order',  4, 'RemoveDC', 'on' );

    % calculate responses
    EEG = pop_epoch(EEG, {'S  1'}, [-1.1  1.5], 'epochinfo', 'yes');
	EEG = pop_rmbase(EEG, [-800 -600]);
    [EEG rejidx] = pop_autorej(EEG, 'nogui','on','threshold',500,'eegplot','off', 'startprob', 5, 'maxrej', 10);
    [EEG, ~, reacmat_a, reacmat_v, accumat_a, accumat_v ,trimat_a, trimat_v, indexstruct] = gs_reaccalc(EEG,subcond);

    % single task
    Idx_st_a_single = [indexstruct.st_a_s_c,indexstruct.st_a_s_i];
    Idx_st_v_single = [indexstruct.st_v_s_c,indexstruct.st_v_s_i];
    Idx_wa_a_single = [indexstruct.wa_a_s_c,indexstruct.wa_a_s_i];
    Idx_wa_v_single = [indexstruct.wa_v_s_c,indexstruct.wa_v_s_i];
    Idx_pe_a_single = [indexstruct.pe_a_s_c,indexstruct.pe_a_s_i];
    Idx_pe_v_single = [indexstruct.pe_v_s_c,indexstruct.pe_v_s_i];

    % dual task
    Idx_st_a_repeat = [indexstruct.st_a_d_c_r,indexstruct.st_a_d_i_r];
    Idx_st_a_switch =  [indexstruct.st_a_d_c_s,indexstruct.st_a_d_i_s];
    Idx_st_v_repeat = [indexstruct.st_v_d_c_r,indexstruct.st_v_d_i_r];
    Idx_st_v_switch =  [indexstruct.st_v_d_c_s,indexstruct.st_v_d_i_s];
    Idx_wa_a_repeat = [indexstruct.wa_a_d_c_r,indexstruct.wa_a_d_i_r];
    Idx_wa_a_switch =  [indexstruct.wa_a_d_c_s,indexstruct.wa_a_d_i_s];
    Idx_wa_v_repeat = [indexstruct.wa_v_d_c_r,indexstruct.wa_v_d_i_r];
    Idx_wa_v_switch =  [indexstruct.wa_v_d_c_s,indexstruct.wa_v_d_i_s];
    Idx_pe_a_repeat = [indexstruct.pe_a_d_c_r,indexstruct.pe_a_d_i_r];
    Idx_pe_a_switch =  [indexstruct.pe_a_d_c_s,indexstruct.pe_a_d_i_s];
    Idx_pe_v_repeat = [indexstruct.pe_v_d_c_r,indexstruct.pe_v_d_i_r];
    Idx_pe_v_switch =  [indexstruct.pe_v_d_c_s,indexstruct.pe_v_d_i_s];

     save([PATH_STAT '/beha_data/' SUBJ '_accumat_a.mat'],'accumat_a');
     save([PATH_STAT '/beha_data/' SUBJ '_accumat_v.mat'],'accumat_v');
     save([PATH_STAT '/beha_data/' SUBJ '_reacmat_a.mat'],'reacmat_a');
     save([PATH_STAT '/beha_data/' SUBJ '_reacmat_v.mat'],'reacmat_v');
     save([PATH_STAT '/beha_data/' SUBJ '_trimat_a.mat'],'trimat_a');
     save([PATH_STAT '/beha_data/' SUBJ '_trimat_v.mat'],'trimat_v');
     save([PATH_STAT '/analysis/' SUBJ '_rejidx.mat'],'rejidx');

    % calculate erps and stuff
    [erp_cue, erp_time, ersp_cue, erspdb_cue, zersp_cue, itpc_cue, tf_time, tf_frqs] = gs_geterpnerspnitpc_baselineadded(EEG,{Idx_st_a_single,Idx_st_a_repeat,Idx_st_a_switch,...
                                                                                                                              Idx_st_v_single,Idx_st_v_repeat,Idx_st_v_switch,...
                                                                                                                              Idx_wa_a_single,Idx_wa_a_repeat,Idx_wa_a_switch,...
                                                                                                                              Idx_wa_v_single,Idx_wa_v_repeat,Idx_wa_v_switch,...
                                                                                                                              Idx_pe_a_single,Idx_pe_a_repeat,Idx_pe_a_switch,...
                                                                                                                              Idx_pe_v_single,Idx_pe_v_repeat,Idx_pe_v_switch},...
                                                                                                                              'blerp',[-800,-600],'blersp',[-800,-600],'BL_switch',2);

    [erp_tar, erp_time, ersp_tar, erspdb_tar, zersp_tar, itpc_tar, tf_time, tf_frqs] = gs_geterpnerspnitpc_baselineadded(EEG,{Idx_st_a_single,Idx_st_a_repeat,Idx_st_a_switch,...
                                                                                                                              Idx_st_v_single,Idx_st_v_repeat,Idx_st_v_switch,...
                                                                                                                              Idx_wa_a_single,Idx_wa_a_repeat,Idx_wa_a_switch,...
                                                                                                                              Idx_wa_v_single,Idx_wa_v_repeat,Idx_wa_v_switch,...
                                                                                                                              Idx_pe_a_single,Idx_pe_a_repeat,Idx_pe_a_switch,...
                                                                                                                              Idx_pe_v_single,Idx_pe_v_repeat,Idx_pe_v_switch},...
                                                                                                                              'blerp',[-200,0],'blersp',[-200,-0],'BL_switch',2);

    % save the data - cue
    save([PATH_STAT '/erp_data/' SUBJ '_erp_cue.mat'],'erp_cue');
    save([PATH_STAT '/ersp_data/' SUBJ '_ersp_cue.mat'],'ersp_tar');
    save([PATH_STAT '/ersp_data/' SUBJ '_erspdb_cue.mat'],'erspdb_tar');
    
    % save the data - target
    save([PATH_STAT '/erp_data/' SUBJ '_erp_tar.mat'],'erp_tar');
    save([PATH_STAT '/ersp_data/' SUBJ '_ersp_tar.mat'],'ersp_tar');
    save([PATH_STAT '/ersp_data/' SUBJ '_erspdb_tar.mat'],'erspdb_tar');

    % save overall information
    save([PATH_STAT '/ersp_data/erptime.mat'],'erp_time');
    save([PATH_STAT '/ersp_data/ersptime.mat'],'tf_time');

end % subj

%% now merge
for subi = 1:length(allsubjlist)
    
    % clear variables
    clear erp ersp erspdb SUBJ

    % create naming variable
    SUBJ = allsubjlist{subi};
    % load data
    load([PATH_STAT '/erp_data/' SUBJ '_erp_cue.mat'],'erp_cue');
    load([PATH_STAT '/erp_data/' SUBJ '_erp_tar.mat'],'erp_tar');

    % merge data
    allerps_cue(subi,:,:,:) = erp_cue;
    allerps_tar(subi,:,:,:) = erp_tar;

end %subi

save([PATH_STAT '/erp_data/allerps_cue.mat'],'allerps_cue')
save([PATH_STAT '/erp_data/allerps_tar.mat'],'allerps_tar')