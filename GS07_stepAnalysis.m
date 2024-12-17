%% GrailSwitch: Analyze step data
% path information
clear all; clc;

% define working machine
macID = 2; % 1: work, 2: home

if ispc
    PATH = ['D:\Data\grailswitch'];
elseif ismac
    if macID == 1
        PATH = ['/Volumes/Work2TB/Seafile/grailswitch/'];
    elseif macID == 2
        PATH = ['/Users/julianreiser/Seafile/grailswitch/'];
    end
elseif isunix
    PATH = '/mnt/data/reiserOn8TB/grailswitch';
end

% define subjects
subjlist = {'02','03','04','05','06','07','09','10','11','12','13','15','16','18','19','20','21','22','23','24','25','27'};

% define factors
movelist = {'ST','WA','PE'};
modalist = {'V','A'};
tasklist = {'single','mixed'};

%% run for all subjects
for subi = 1:length(subjlist)
    %% load data
    % define subject
    SUB = subjlist{subi};
    % load data for left and right leg step cycles
    load([PATH 'data/walk/gs_vp' SUB '_stepmatl.mat'],'stepmatl')
    load([PATH 'data/walk/gs_vp' SUB '_stepmatr.mat'],'stepmatr')

    %% create larger matrix
    % define factor names
    condlist = {'st_v_single','st_v_mixed','st_a_single','st_a_mixed',...
                'wa_v_single','wa_v_mixed','wa_a_single','wa_a_mixed',...
                'pe_v_single','pe_v_mixed','pe_a_single','pe_a_mixed'};
    
    % full model
    stepnuml(:,subi) = stepmatl(:,8);
    steptimemeanl(:,subi) = stepmatl(:,9);
    steptimevarl(:,subi) = stepmatl(:,10);
    
    stepnumr(:,subi) = stepmatr(:,8);
    steptimemeanr(:,subi) = stepmatr(:,9);
    steptimevarr(:,subi) = stepmatr(:,10);

    % switch-blocks only model
    stepnuml_sbo(:,subi) = stepmatl(2:2:end,8);
    steptimemeanl_sbo(:,subi) = stepmatl(2:2:end,9);
    steptimevarl_sbo(:,subi) = stepmatl(2:2:end,10);
    
    stepnumr_sbo(:,subi) = stepmatr(2:2:end,8);
    steptimemeanr_sbo(:,subi) = stepmatr(2:2:end,9);
    steptimevarr_sbo(:,subi) = stepmatr(2:2:end,10);

end % for subi

%% create anova vars for full model
multcomptype = 'lsd';
full_length = (length(movelist)-1) * length(modalist) * length(tasklist);

% Create a table reflecting the within subject factors
MOVE = cell(full_length,1); % lat ring conditions: 1, 5
MODA = cell(full_length,1); % lat move conditions: stand, walk, pert
TASK = cell(full_length,1); % lat ecce conditions: 0, 18, 40

% Assiging the values to the parameters based on the data sorting
c1 = cell(1,1); c1{1} = movelist{2}; c1 = repmat(c1,full_length/2,1); MOVE(1:4,1) = c1;
c1 = cell(1,1); c1{1} = movelist{3}; c1 = repmat(c1,full_length/2,1); MOVE(5:8,1) = c1;

c1 = cell(1,1); c1{1} = modalist{1}; c1 = repmat(c1,full_length/2,1); MODA([1:2,5:6],1) = c1;
c1 = cell(1,1); c1{1} = modalist{2}; c1 = repmat(c1,full_length/2,1); MODA([3:4,7:8],1) = c1;

c1 = cell(1,1); c1{1} = tasklist{1}; c1 = repmat(c1,full_length/2,1); TASK([1:2:end],1) = c1;
c1 = cell(1,1); c1{1} = tasklist{2}; c1 = repmat(c1,full_length/2,1); TASK([2:2:end],1) = c1;

% Create the within table
factorNames = {'MOVE','MODA','TASK'};
within = table(MOVE,MODA,TASK, 'VariableNames', factorNames);

% now create ANOVA table
varNames = cell(full_length,1);
for i = 1 : full_length
    v = strcat('V',num2str(i));
    varNames{i,1} = v;
end% Create a table storing the respones

% create table from data and headings for ANOVA
stepnumlTable = array2table(stepnuml(5:end,:)', 'VariableNames',varNames);
steptimemeanlTable = array2table(steptimemeanl(5:end,:)', 'VariableNames',varNames);
steptimevarlTable = array2table(steptimevarl(5:end,:)', 'VariableNames',varNames);

stepnumrTable = array2table(stepnumr(5:end,:)', 'VariableNames',varNames);
steptimemeanrTable = array2table(steptimemeanr(5:end,:)', 'VariableNames',varNames);
steptimevarrTable = array2table(steptimevarr(5:end,:)', 'VariableNames',varNames);


%% run ANOVA left
% number of steps - fit the repeated measures model
stepnuml_rm = fitrm(stepnumlTable,'V1-V8~1','WithinDesign',within);
[rmStat.stepnuml.rmanova] = ranova(stepnuml_rm, 'WithinModel','MOVE*MODA*TASK');
%[FCzN2_amp_move FCzN2_amp_ecce FCzN2_amp_int rmStat.FCzN2.ampMOVE rmStat.FCzN2.ampECCE rmStat.FCzN2.ampINT FCzN2_amp_movestat FCzN2_amp_eccestat FCzN2_amp_intstat] = lw_phtest(n1amps1,within,2);
lw_qq(stepnuml(5:end,:)',4,2,stepnuml_rm.mauchly)
rmStat.stepnuml.rmanova = lw_fdr(rmStat.stepnuml.rmanova,'corrMethod',1);
rmStat.stepnuml.rmanova = lw_eta(rmStat.stepnuml.rmanova);
rmStat.stepnuml.mauchly = stepnuml_rm.mauchly;

% inter-step time mean - fit the repeated measures model
steptimemeanl_rm = fitrm(steptimemeanlTable,'V1-V8~1','WithinDesign',within);
[rmStat.steptimemeanl.rmanova] = ranova(steptimemeanl_rm, 'WithinModel','MOVE*MODA*TASK');
%[FCzN2_amp_move FCzN2_amp_ecce FCzN2_amp_int rmStat.FCzN2.ampMOVE rmStat.FCzN2.ampECCE rmStat.FCzN2.ampINT FCzN2_amp_movestat FCzN2_amp_eccestat FCzN2_amp_intstat] = lw_phtest(n1amps1,within,2);
lw_qq(steptimemeanl(5:end,:)',4,2,steptimemeanl_rm.mauchly)
rmStat.steptimemeanl.rmanova = lw_fdr(rmStat.steptimemeanl.rmanova,'corrMethod',1);
rmStat.steptimemeanl.rmanova = lw_eta(rmStat.steptimemeanl.rmanova);
rmStat.steptimemeanl.mauchly = steptimemeanl_rm.mauchly;

% inter-step time variability - fit the repeated measures model
steptimevarl_rm = fitrm(steptimevarlTable,'V1-V8~1','WithinDesign',within);
[rmStat.steptimevarl.rmanova] = ranova(steptimevarl_rm, 'WithinModel','MOVE*MODA*TASK');
%[FCzN2_amp_move FCzN2_amp_ecce FCzN2_amp_int rmStat.FCzN2.ampMOVE rmStat.FCzN2.ampECCE rmStat.FCzN2.ampINT FCzN2_amp_movestat FCzN2_amp_eccestat FCzN2_amp_intstat] = lw_phtest(n1amps1,within,2);
lw_qq(steptimevarl(5:end,:)',4,2,steptimevarl_rm.mauchly)
rmStat.steptimevarl.rmanova = lw_fdr(rmStat.steptimevarl.rmanova,'corrMethod',1);
rmStat.steptimevarl.rmanova = lw_eta(rmStat.steptimevarl.rmanova);
rmStat.steptimevarl.mauchly = steptimevarl_rm.mauchly;
close all
%% run ANOVA right
% number of steps - fit the repeated measures model
stepnumr_rm = fitrm(stepnumrTable,'V1-V8~1','WithinDesign',within);
[rmStat.stepnumr.rmanova] = ranova(stepnumr_rm, 'WithinModel','MOVE*MODA*TASK');
%[FCzN2_amp_move FCzN2_amp_ecce FCzN2_amp_int rmStat.FCzN2.ampMOVE rmStat.FCzN2.ampECCE rmStat.FCzN2.ampINT FCzN2_amp_movestat FCzN2_amp_eccestat FCzN2_amp_intstat] = lw_phtest(n1amps1,within,2);
lw_qq(stepnumr(5:end,:)',4,2,stepnumr_rm.mauchly)
rmStat.stepnumr.rmanova = lw_fdr(rmStat.stepnumr.rmanova,'corrMethod',1);
rmStat.stepnumr.rmanova = lw_eta(rmStat.stepnumr.rmanova);
rmStat.stepnumr.mauchly = stepnumr_rm.mauchly;

% inter-step time mean - fit the repeated measures model
stepteimemeanr_rm = fitrm(steptimemeanrTable,'V1-V8~1','WithinDesign',within);
[rmStat.stepteimemeanr.rmanova] = ranova(stepteimemeanr_rm, 'WithinModel','MOVE*MODA*TASK');
%[FCzN2_amp_move FCzN2_amp_ecce FCzN2_amp_int rmStat.FCzN2.ampMOVE rmStat.FCzN2.ampECCE rmStat.FCzN2.ampINT FCzN2_amp_movestat FCzN2_amp_eccestat FCzN2_amp_intstat] = lw_phtest(n1amps1,within,2);
lw_qq(steptimemeanr(5:end,:)',4,2,stepteimemeanr_rm.mauchly)
rmStat.stepteimemeanr.rmanova = lw_fdr(rmStat.stepteimemeanr.rmanova,'corrMethod',1);
rmStat.stepteimemeanr.rmanova = lw_eta(rmStat.stepteimemeanr.rmanova);
rmStat.stepteimemeanr.mauchly = stepteimemeanr_rm.mauchly;

% inter-step time variability - fit the repeated measures model
stepteimevarr_rm = fitrm(steptimevarrTable,'V1-V8~1','WithinDesign',within);
[rmStat.stepteimevarr.rmanova] = ranova(stepteimevarr_rm, 'WithinModel','MOVE*MODA*TASK');
%[FCzN2_amp_move FCzN2_amp_ecce FCzN2_amp_int rmStat.FCzN2.ampMOVE rmStat.FCzN2.ampECCE rmStat.FCzN2.ampINT FCzN2_amp_movestat FCzN2_amp_eccestat FCzN2_amp_intstat] = lw_phtest(n1amps1,within,2);
lw_qq(steptimevarr(5:end,:)',4,2,stepteimevarr_rm.mauchly)
rmStat.stepteimevarr.rmanova = lw_fdr(rmStat.stepteimevarr.rmanova,'corrMethod',1);
rmStat.stepteimevarr.rmanova = lw_eta(rmStat.stepteimevarr.rmanova);
rmStat.stepteimevarr.mauchly = stepteimevarr_rm.mauchly;
close all

%% create anova vars for switch block only model
multcomptype = 'lsd';
full_length = (length(movelist)-1) * length(modalist);

% Create a table reflecting the within subject factors
MOVE = cell(full_length,1); % lat ring conditions: 1, 5
MODA = cell(full_length,1); % lat move conditions: stand, walk, pert

% Assiging the values to the parameters based on the data sorting
c1 = cell(1,1); c1{1} = movelist{2}; c1 = repmat(c1,full_length/2,1); MOVE(1:2,1) = c1;
c1 = cell(1,1); c1{1} = movelist{3}; c1 = repmat(c1,full_length/2,1); MOVE(3:4,1) = c1;

c1 = cell(1,1); c1{1} = modalist{1}; c1 = repmat(c1,full_length/2,1); MODA(1:2:end,1) = c1;
c1 = cell(1,1); c1{1} = modalist{2}; c1 = repmat(c1,full_length/2,1); MODA(2:2:end,1) = c1;

% Create the within table
factorNames = {'MOVE','MODA'};
within = table(MOVE,MODA, 'VariableNames', factorNames);

% now create ANOVA table
varNames = cell(full_length,1);
for i = 1 : full_length
    v = strcat('V',num2str(i));
    varNames{i,1} = v;
end% Create a table storing the respones

% create table from data and headings for ANOVA
stepnumlTable_sbo = array2table(stepnuml_sbo(3:end,:)', 'VariableNames',varNames);
steptimemeanlTable_sbo = array2table(steptimemeanl_sbo(3:end,:)', 'VariableNames',varNames);
steptimevarlTable_sbo = array2table(steptimevarl_sbo(3:end,:)', 'VariableNames',varNames);

stepnumrTable_sbo = array2table(stepnumr_sbo(3:end,:)', 'VariableNames',varNames);
steptimemeanrTable_sbo = array2table(steptimemeanr_sbo(3:end,:)', 'VariableNames',varNames);
steptimevarrTable_sbo = array2table(steptimevarr_sbo(3:end,:)', 'VariableNames',varNames);


%% run ANOVA left - switch block only
% number of steps - fit the repeated measures model
stepnuml_rm_sbo = fitrm(stepnumlTable_sbo,'V1-V4~1','WithinDesign',within);
[rmStat.stepnuml_sbo.rmanova] = ranova(stepnuml_rm_sbo, 'WithinModel','MOVE*MODA');
%[FCzN2_amp_move FCzN2_amp_ecce FCzN2_amp_int rmStat.FCzN2.ampMOVE rmStat.FCzN2.ampECCE rmStat.FCzN2.ampINT FCzN2_amp_movestat FCzN2_amp_eccestat FCzN2_amp_intstat] = lw_phtest(n1amps1,within,2);
lw_qq(stepnuml_sbo(3:end,:)',2,2,stepnuml_rm_sbo.mauchly)
rmStat.stepnuml_sbo.rmanova = lw_fdr(rmStat.stepnuml_sbo.rmanova,'corrMethod',1);
rmStat.stepnuml_sbo.rmanova = lw_eta(rmStat.stepnuml_sbo.rmanova);
rmStat.stepnuml_sbo.mauchly = stepnuml_rm_sbo.mauchly;

% inter-step time mean - fit the repeated measures model
stepteimemeanl_rm_sbo = fitrm(steptimemeanlTable_sbo,'V1-V4~1','WithinDesign',within);
[rmStat.stepteimemeanl_sbo.rmanova] = ranova(stepteimemeanl_rm_sbo, 'WithinModel','MOVE*MODA');
%[FCzN2_amp_move FCzN2_amp_ecce FCzN2_amp_int rmStat.FCzN2.ampMOVE rmStat.FCzN2.ampECCE rmStat.FCzN2.ampINT FCzN2_amp_movestat FCzN2_amp_eccestat FCzN2_amp_intstat] = lw_phtest(n1amps1,within,2);
lw_qq(steptimemeanl_sbo(3:end,:)',2,2,stepteimemeanl_rm_sbo.mauchly)
rmStat.stepteimemeanl_sbo.rmanova = lw_fdr(rmStat.stepteimemeanl_sbo.rmanova,'corrMethod',1);
rmStat.stepteimemeanl_sbo.rmanova = lw_eta(rmStat.stepteimemeanl_sbo.rmanova);
rmStat.stepteimemeanl_sbo.mauchly = stepteimemeanl_rm_sbo.mauchly;

% inter-step time variability - fit the repeated measures model
stepteimevarl_rm_sbo = fitrm(steptimevarlTable_sbo,'V1-V4~1','WithinDesign',within);
[rmStat.stepteimevarl_sbo.rmanova] = ranova(stepteimevarl_rm_sbo, 'WithinModel','MOVE*MODA');
%[FCzN2_amp_move FCzN2_amp_ecce FCzN2_amp_int rmStat.FCzN2.ampMOVE rmStat.FCzN2.ampECCE rmStat.FCzN2.ampINT FCzN2_amp_movestat FCzN2_amp_eccestat FCzN2_amp_intstat] = lw_phtest(n1amps1,within,2);
lw_qq(steptimevarl_sbo(3:end,:)',2,2,stepteimevarl_rm_sbo.mauchly)
rmStat.stepteimevarl_sbo.rmanova = lw_fdr(rmStat.stepteimevarl_sbo.rmanova,'corrMethod',1);
rmStat.stepteimevarl_sbo.rmanova = lw_eta(rmStat.stepteimevarl_sbo.rmanova);
rmStat.stepteimevarl_sbo.mauchly = stepteimevarl_rm_sbo.mauchly;
close all

%% run ANOVA right - switch block only
% number of steps - fit the repeated measures model
stepnumr_rm_sbo = fitrm(stepnumrTable_sbo,'V1-V4~1','WithinDesign',within);
[rmStat.stepnumr_sbo.rmanova] = ranova(stepnumr_rm, 'WithinModel','MOVE*MODA');
%[FCzN2_amp_move FCzN2_amp_ecce FCzN2_amp_int rmStat.FCzN2.ampMOVE rmStat.FCzN2.ampECCE rmStat.FCzN2.ampINT FCzN2_amp_movestat FCzN2_amp_eccestat FCzN2_amp_intstat] = lw_phtest(n1amps1,within,2);
lw_qq(stepnumr_sbo(3:end,:)',2,2,stepnumr_rm_sbo.mauchly)
rmStat.stepnumr_sbo.rmanova = lw_fdr(rmStat.stepnumr_sbo.rmanova,'corrMethod',1);
rmStat.stepnumr_sbo.rmanova = lw_eta(rmStat.stepnumr_sbo.rmanova);
rmStat.stepnumr_sbo.mauchly = stepnumr_rm_sbo.mauchly;

% inter-step time mean - fit the repeated measures model
stepteimemeanr_rm_sbo = fitrm(steptimemeanrTable_sbo,'V1-V4~1','WithinDesign',within);
[rmStat.stepteimemeanr_sbo.rmanova] = ranova(stepteimemeanr_rm, 'WithinModel','MOVE*MODA');
%[FCzN2_amp_move FCzN2_amp_ecce FCzN2_amp_int rmStat.FCzN2.ampMOVE rmStat.FCzN2.ampECCE rmStat.FCzN2.ampINT FCzN2_amp_movestat FCzN2_amp_eccestat FCzN2_amp_intstat] = lw_phtest(n1amps1,within,2);
lw_qq(steptimemeanr_sbo(3:end,:)',2,2,stepteimemeanr_rm_sbo.mauchly)
rmStat.stepteimemeanr_sbo.rmanova = lw_fdr(rmStat.stepteimemeanr_sbo.rmanova,'corrMethod',1);
rmStat.stepteimemeanr_sbo.rmanova = lw_eta(rmStat.stepteimemeanr_sbo.rmanova);
rmStat.stepteimemeanr_sbo.mauchly = stepteimemeanr_rm_sbo.mauchly;

% inter-step time variability - fit the repeated measures model
stepteimevarr_rm_sbo = fitrm(steptimevarrTable_sbo,'V1-V4~1','WithinDesign',within);
[rmStat.stepteimevarr_sbo.rmanova] = ranova(stepteimevarr_rm, 'WithinModel','MOVE*MODA');
%[FCzN2_amp_move FCzN2_amp_ecce FCzN2_amp_int rmStat.FCzN2.ampMOVE rmStat.FCzN2.ampECCE rmStat.FCzN2.ampINT FCzN2_amp_movestat FCzN2_amp_eccestat FCzN2_amp_intstat] = lw_phtest(n1amps1,within,2);
lw_qq(steptimevarr_sbo(3:end,:)',2,2,stepteimevarr_rm_sbo.mauchly)
rmStat.stepteimevarr_sbo.rmanova = lw_fdr(rmStat.stepteimevarr_sbo.rmanova,'corrMethod',1);
rmStat.stepteimevarr_sbo.rmanova = lw_eta(rmStat.stepteimevarr_sbo.rmanova);
rmStat.stepteimevarr_sbo.mauchly = stepteimevarr_rm_sbo.mauchly;
close all