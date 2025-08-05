%% calculate response related statistics
clear all

%% load paths
if ispc
    PATH = 'D:\Data\grailswitch\';
elseif ismac
    PATH = '/Volumes/Work4TB/Seafile/grailswitch/';
end

PATH_STAT = [PATH '/stats_rework'];
PATH_PLOT = [PATH '/plots/beha']

% delare subjects
subjlist = {'02','03','04','05','06','07','09','10','11','12','13','15','16','18','19','20','21','22','23','24','25','27'};

%% load data
for subj = 1:length(subjlist)
    SUBJ = subjlist{subj};

    tmpdata_a = [];
    tmpdata_v = [];
    tmpdata_a = load([PATH_STAT '/beha_data/' SUBJ '_accumat_a.mat']);
    tmpdata_v = load([PATH_STAT '/beha_data/' SUBJ '_accumat_v.mat']);
    tmpreac_a = load([PATH_STAT '/beha_data/' SUBJ '_reacmat_a.mat']);
    tmpreac_v = load([PATH_STAT '/beha_data/' SUBJ '_reacmat_v.mat']);
    tmptri_a = load([PATH_STAT '/beha_data/' SUBJ '_trimat_a.mat']);
    tmptri_v = load([PATH_STAT '/beha_data/' SUBJ '_trimat_v.mat']);
    tmpmat_a(subj,:,:) = tmpdata_a.accumat_a;
    tmpmat_v(subj,:,:) = tmpdata_v.accumat_v;
    tmpmatr_a(subj,:,:) = tmpreac_a.reacmat_a;
    tmpmatr_v(subj,:,:) = tmpreac_v.reacmat_v;
    tmpmatt_a(subj,:,:) = tmptri_a.trimat_a;
    tmpmatt_v(subj,:,:) = tmptri_v.trimat_v;

end % subj

% create means
accumean_a = squeeze(mean(tmpmat_a,1,'omitnan'));
accumean_v = squeeze(mean(tmpmat_v,1,'omitnan'));
accustd_a = squeeze(std(tmpmat_a,1,'omitnan'));
accustd_v = squeeze(std(tmpmat_v,1,'omitnan'));
reacmean_a = squeeze(mean(tmpmatr_a,1,'omitnan'));
reacmean_v = squeeze(mean(tmpmatr_v,1,'omitnan'));
reacstd_a = squeeze(std(tmpmatr_a,1,'omitnan'));
reacstd_v = squeeze(std(tmpmatr_v,1,'omitnan'));
trimean_a = squeeze(mean(tmpmatt_a,1,'omitnan'));
trimean_v = squeeze(mean(tmpmatt_v,1,'omitnan'));
tristd_a = squeeze(std(tmpmatt_a,1,'omitnan'));
tristd_v = squeeze(std(tmpmatt_v,1,'omitnan'));

% description for reacmat & accumat cols (rows always standt,walk,pert):
%  1 - single task congruent
%  2 - single task incongruent
%  3 - dual task congruent
%  4 - dual task incongruent
%  5 - dual task congruent repeat
%  6 - dual task incongruent repeat
%  7 - dual task congruent switch
%  8 - dual task incongruent switch
%
% ----- 9-16 only for reacmat -----
%
%  9 - single task congruent incorrect
% 10 - single task incongruent incorrect
% 11 - dual task congruent incorrect
% 12 - dual task incongruent incorrect
% 13 - dual task congruent repeat incorrect
% 14 - dual task incongruent repeat incorrect
% 15 - dual task congruent switch incorrect
% 16 - dual task incongruent switch incorrect

%% create anovamats
% 1:6 - stand auditory single, stand auditory repeat, stand auditory switch, stand visual single, stand visual repeat, stand visual switch
% 7:12 - walk auditory single, walk auditory repeat, walk auditory switch, walk visual single, walk visual repeat, walk visual switch
% 13:18 - pert auditory single, pert auditory repeat, pert auditory switch, pert visual single, pert visual repeat, pert visual switch
clear anovaReacMat

% RT
anovaReacMat(:,1) = mean(squeeze(tmpmatr_a(:,1:2,1)),2,'omitnan');
anovaReacMat(:,2) = mean(squeeze(tmpmatr_a(:,5:6,1)),2,'omitnan');
anovaReacMat(:,3) = mean(squeeze(tmpmatr_a(:,7:8,1)),2,'omitnan');
anovaReacMat(:,4) = mean(squeeze(tmpmatr_v(:,1:2,1)),2,'omitnan');
anovaReacMat(:,5) = mean(squeeze(tmpmatr_v(:,5:6,1)),2,'omitnan');
anovaReacMat(:,6) = mean(squeeze(tmpmatr_v(:,7:8,1)),2,'omitnan');

anovaReacMat(:,7) = mean(squeeze(tmpmatr_a(:,1:2,2)),2,'omitnan');
anovaReacMat(:,8) = mean(squeeze(tmpmatr_a(:,5:6,2)),2,'omitnan');
anovaReacMat(:,9) = mean(squeeze(tmpmatr_a(:,7:8,2)),2,'omitnan');
anovaReacMat(:,10) = mean(squeeze(tmpmatr_v(:,1:2,2)),2,'omitnan');
anovaReacMat(:,11) = mean(squeeze(tmpmatr_v(:,5:6,2)),2,'omitnan');
anovaReacMat(:,12) = mean(squeeze(tmpmatr_v(:,7:8,2)),2,'omitnan');

anovaReacMat(:,13) = mean(squeeze(tmpmatr_a(:,1:2,3)),2,'omitnan');
anovaReacMat(:,14) = mean(squeeze(tmpmatr_a(:,5:6,3)),2,'omitnan');
anovaReacMat(:,15) = mean(squeeze(tmpmatr_a(:,7:8,3)),2,'omitnan');
anovaReacMat(:,16) = mean(squeeze(tmpmatr_v(:,1:2,3)),2,'omitnan');
anovaReacMat(:,17) = mean(squeeze(tmpmatr_v(:,5:6,3)),2,'omitnan');
anovaReacMat(:,18) = mean(squeeze(tmpmatr_v(:,7:8,3)),2,'omitnan');

% accuracy
anovaAccuMat(:,1) = mean(squeeze(tmpmat_a(:,1:2,1)),2,'omitnan');
anovaAccuMat(:,2) = mean(squeeze(tmpmat_a(:,5:6,1)),2,'omitnan');
anovaAccuMat(:,3) = mean(squeeze(tmpmat_a(:,7:8,1)),2,'omitnan');
anovaAccuMat(:,4) = mean(squeeze(tmpmat_v(:,1:2,1)),2,'omitnan');
anovaAccuMat(:,5) = mean(squeeze(tmpmat_v(:,5:6,1)),2,'omitnan');
anovaAccuMat(:,6) = mean(squeeze(tmpmat_v(:,7:8,1)),2,'omitnan');

anovaAccuMat(:,7) = mean(squeeze(tmpmat_a(:,1:2,2)),2,'omitnan');
anovaAccuMat(:,8) = mean(squeeze(tmpmat_a(:,5:6,2)),2,'omitnan');
anovaAccuMat(:,9) = mean(squeeze(tmpmat_a(:,7:8,2)),2,'omitnan');
anovaAccuMat(:,10) = mean(squeeze(tmpmat_v(:,1:2,2)),2,'omitnan');
anovaAccuMat(:,11) = mean(squeeze(tmpmat_v(:,5:6,2)),2,'omitnan');
anovaAccuMat(:,12) = mean(squeeze(tmpmat_v(:,7:8,2)),2,'omitnan');

anovaAccuMat(:,13) = mean(squeeze(tmpmat_a(:,1:2,3)),2,'omitnan');
anovaAccuMat(:,14) = mean(squeeze(tmpmat_a(:,5:6,3)),2,'omitnan');
anovaAccuMat(:,15) = mean(squeeze(tmpmat_a(:,7:8,3)),2,'omitnan');
anovaAccuMat(:,16) = mean(squeeze(tmpmat_v(:,1:2,3)),2,'omitnan');
anovaAccuMat(:,17) = mean(squeeze(tmpmat_v(:,5:6,3)),2,'omitnan');
anovaAccuMat(:,18) = mean(squeeze(tmpmat_v(:,7:8,3)),2,'omitnan');

% prepare anova description
anovaTitles = {'st_a_single','st_a_repeat','st_a_switch',...
               'st_v_single','st_v_repeat','st_v_switch',...
               'wa_a_single','wa_a_repeat','wa_a_switch',...
               'wa_v_single','wa_v_repeat','wa_v_switch',...
               'pa_a_single','pa_a_repeat','pa_a_switch',...
               'pa_v_single','pa_v_repeat','pa_v_switch'};

% plot response times
anovaReacMeans = mean(anovaReacMat,1);
anovaReacSEMs = std(anovaReacMat,1) / sqrt(size(anovaReacMat,1));

anovaAccuMeans = mean(anovaAccuMat,1)*100;
anovaAccuSEMs = std(anovaAccuMat,1) / sqrt(size(anovaAccuMat,1))*100;

figure;
subplot(221)
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 40, 20], 'PaperUnits', 'centimeters', 'PaperSize', [29.7,21])
lb_niceInteractionPlot([anovaReacMeans(1:6:end);anovaReacMeans(2:6:end);anovaReacMeans(3:6:end)], ...
                           [anovaReacSEMs(1:6:end);anovaReacSEMs(2:6:end);anovaReacSEMs(3:6:end)], ...
                            'Response times - Auditory','Time in ms',{'repeat-only','repeat','switch'},[400,700],'label_legend',{'stand','walk','perturbation'},'locs_legend','best')
subplot(222)
lb_niceInteractionPlot([anovaReacMeans(4:6:end);anovaReacMeans(5:6:end);anovaReacMeans(6:6:end)], ...
                           [anovaReacSEMs(4:6:end);anovaReacSEMs(5:6:end);anovaReacSEMs(6:6:end)], ...
                            'Response times - Visual','Time in ms',{'repeat-only','repeat','switch'},[400,700])

subplot(223)
lb_niceInteractionPlot([anovaAccuMeans(1:6:end);anovaAccuMeans(2:6:end);anovaAccuMeans(3:6:end)], ...
                           [anovaAccuSEMs(1:6:end);anovaAccuSEMs(2:6:end);anovaAccuSEMs(3:6:end)], ...
                            'Accuracy - Auditory','% correct',{'repeat-only','repeat','switch'},[80,100])
subplot(224)
lb_niceInteractionPlot([anovaAccuMeans(4:6:end);anovaAccuMeans(5:6:end);anovaAccuMeans(6:6:end)], ...
                           [anovaAccuSEMs(4:6:end);anovaAccuSEMs(5:6:end);anovaAccuSEMs(6:6:end)], ...
                            'Accuracy - Visual','% correct',{'repeat-only','repeat','switch'},[80,100])

saveas(gcf,[PATH_PLOT '/accu_means.png']);
saveas(gcf,[PATH_PLOT '/accu_means.fig']); close gcf;
%% ---------------------- now for the rmANOVAs ---------------------------
anova_walklist = 1:3;
anova_modalist = 1:2;
anova_tasklist = 1:3;

multcomptype = 'lsd';
full_length = length(anova_walklist) * length(anova_modalist) * length(anova_tasklist);

% Create a table reflecting the within subject factors
MOVE = cell(full_length,1); % lat ring conditions: 1, 5
MODA = cell(full_length,1); % lat move conditions: stand, walk, pert
TASK = cell(full_length,1); % lat ecce conditions: 0, 18, 40

% Assiging the values to the parameters based on the data sorting
c1 = cell(1,1); c1{1} = 'st'; c1 = repmat(c1,full_length/3,1); MOVE(1:6,1) = c1;
c1 = cell(1,1); c1{1} = 'wa'; c1 = repmat(c1,full_length/3,1); MOVE(7:12,1) = c1;
c1 = cell(1,1); c1{1} = 'pe'; c1 = repmat(c1,full_length/3,1); MOVE(13:18,1) = c1;

c1 = cell(1,1); c1{1} = 'a'; c1 = repmat(c1,full_length/2,1); MODA([1:3,7:9,13:15],1) = c1;
c1 = cell(1,1); c1{1} = 'v'; c1 = repmat(c1,full_length/2,1); MODA([4:6,10:12,16:18],1) = c1;

c1 = cell(1,1); c1{1} = 'single'; c1 = repmat(c1,full_length/3,1); TASK([1:3:end],1) = c1;
c1 = cell(1,1); c1{1} = 'repead'; c1 = repmat(c1,full_length/3,1); TASK([2:3:end],1) = c1;
c1 = cell(1,1); c1{1} = 'switch'; c1 = repmat(c1,full_length/3,1); TASK([3:3:end],1) = c1;

% Create the within table
factorNames = {'MOVE','MODA','TASK'};
within = table(MOVE,MODA,TASK, 'VariableNames', factorNames);

% now create ANOVA table
varNames = cell(full_length,1);
for i = 1 : full_length
    v = strcat('V',num2str(i));
    varNames{i,1} = v;
end% Create a table storing the respones
anovaReacTable = array2table(anovaReacMat, 'VariableNames',varNames);
anovaAccuTable = array2table(anovaAccuMat, 'VariableNames',varNames);


%% compute the rmANOVA
clear rmStat

% Response times - fit the repeated measures model
reac_rm = fitrm(anovaReacTable,'V1-V18~1','WithinDesign',within);
[rmStat.reac.rmanova] = ranova(reac_rm, 'WithinModel','MOVE*MODA*TASK');
lw_qq(anovaReacMat,6,3,reac_rm.mauchly)
rmStat.reac.rmanova = lw_fdr(rmStat.reac.rmanova,'corrMethod',1);
rmStat.reac.rmanova = lw_eta(rmStat.reac.rmanova);
rmStat.reac.mauchly = reac_rm.mauchly;
rmStat.reac.phTask = multcompare(reac_rm,'TASK');
rmStat.reac.phModa = multcompare(reac_rm,'MODA');
rmStat.reac.phMoveByTask = multcompare(reac_rm,'MOVE','By','TASK');
rmStat.reac.phModaByTask = multcompare(reac_rm,'MODA','By','TASK');


% REAC
MainTask = [mean(anovaReacMat(:,[1:3:end]),[1,2]),std(anovaReacMat(:,[1:3:end]),0,[1,2]);... % single
            mean(anovaReacMat(:,[2:3:end]),[1,2]),std(anovaReacMat(:,[2:3:end]),0,[1,2]);... % repeat
            mean(anovaReacMat(:,[3:3:end]),[1,2]),std(anovaReacMat(:,[3:3:end]),0,[1,2])]; % switch

MainModa = [mean(anovaReacMat(:,[1:3,7:9,13:15]),[1,2]),std(anovaReacMat(:,[1:3,7:9,13:15]),0,[1,2]); % auditory
            mean(anovaReacMat(:,[4:6,10:12,16:18]),[1,2]),std(anovaReacMat(:,[4:6,10:12,16:18]),0,[1,2])]; % visual

ModaTask = [mean(anovaReacMat(:,[1,7,13]),[1,2]),std(anovaReacMat(:,[1,7,13]),0,[1,2]); % auditory single
            mean(anovaReacMat(:,[2,8,14]),[1,2]),std(anovaReacMat(:,[2,8,14]),0,[1,2]); % auditory repeat
            mean(anovaReacMat(:,[3,9,15]),[1,2]),std(anovaReacMat(:,[3,9,15]),0,[1,2]); % auditory switch
            mean(anovaReacMat(:,[4,10,16]),[1,2]),std(anovaReacMat(:,[4,10,16]),0,[1,2]); % visual single
            mean(anovaReacMat(:,[5,11,17]),[1,2]),std(anovaReacMat(:,[5,11,17]),0,[1,2]); % visual repeat
            mean(anovaReacMat(:,[6,12,18]),[1,2]),std(anovaReacMat(:,[6,12,18]),0,[1,2])]; % visual switch
            
MoveTask = [mean(anovaReacMat(:,[1,4]),[1,2]),std(anovaReacMat(:,[1,4]),0,[1,2]); % stand single
            mean(anovaReacMat(:,[2,5]),[1,2]),std(anovaReacMat(:,[2,5]),0,[1,2]); % stand repeat
            mean(anovaReacMat(:,[3,6]),[1,2]),std(anovaReacMat(:,[3,6]),0,[1,2]); % stand switch
            mean(anovaReacMat(:,[7,10]),[1,2]),std(anovaReacMat(:,[7,10]),0,[1,2]); % walk single
            mean(anovaReacMat(:,[8,11]),[1,2]),std(anovaReacMat(:,[8,11]),0,[1,2]); % walk repeat
            mean(anovaReacMat(:,[9,12]),[1,2]),std(anovaReacMat(:,[9,12]),0,[1,2]); % walk switch
            mean(anovaReacMat(:,[13,16]),[1,2]),std(anovaReacMat(:,[13,16]),0,[1,2]); % pert single
            mean(anovaReacMat(:,[14,17]),[1,2]),std(anovaReacMat(:,[14,17]),0,[1,2]); % pert repeat
            mean(anovaReacMat(:,[15,18]),[1,2]),std(anovaReacMat(:,[15,18]),0,[1,2])]; % pert switch


% Accuracy - fit the repeated measures model
accu_rm = fitrm(anovaAccuTable,'V1-V18~1','WithinDesign',within);
[rmStat.accu.rmanova] = ranova(accu_rm, 'WithinModel','MOVE*MODA*TASK');
%[FCzN2_amp_move FCzN2_amp_ecce FCzN2_amp_int rmStat.FCzN2.ampMOVE rmStat.FCzN2.ampECCE rmStat.FCzN2.ampINT FCzN2_amp_movestat FCzN2_amp_eccestat FCzN2_amp_intstat] = lw_phtest(n1amps1,within,2);
lw_qq(anovaAccuMat,6,3,accu_rm.mauchly)
rmStat.accu.rmanova = lw_fdr(rmStat.accu.rmanova,'corrMethod',1);
rmStat.accu.rmanova = lw_eta(rmStat.accu.rmanova);
rmStat.accu.mauchly = accu_rm.mauchly;
rmStat.accu.phTask = multcompare(accu_rm,'TASK');

% ACCU
MainDiff = [mean(anovaAccuMat(:,[1:3:end]),[1,2]),std(anovaAccuMat(:,[1:3:end]),0,[1,2]);... % single
            mean(anovaAccuMat(:,[2:3:end]),[1,2]),std(anovaAccuMat(:,[2:3:end]),0,[1,2]);... % repeat
            mean(anovaAccuMat(:,[3:3:end]),[1,2]),std(anovaAccuMat(:,[3:3:end]),0,[1,2])]; % switch

close all
%% now save
save([PATH_STAT '/beha_data/anovaReacMat_sbo.mat'],'anovaReacMat')
save([PATH_STAT '/beha_data/anovaAccuMat_sbo.mat'],'anovaAccuMat')


%% Generate comprehensive RT/Accuracy descriptive statistics
fprintf('\n=== GENERATING RT/ACCURACY DESCRIPTIVE STATISTICS ===\n');

% Call the RT/Accuracy descriptive statistics function
rt_accu_descriptive_table = generate_rt_accu_descriptive_stats(...
    anovaReacMat, anovaAccuMat, PATH_STAT);

%% Create publication-ready RT/Accuracy tables

% Separate by measure (RT vs Accuracy)
measures = unique(rt_accu_descriptive_table.Measure);

for meas = 1:length(measures)
    measure_name = measures{meas};
    measure_data = rt_accu_descriptive_table(strcmp(rt_accu_descriptive_table.Measure, measure_name), :);
    
    % Save individual measure tables
    filename = [PATH_STAT '/beha_data/' lower(measure_name) '_descriptive.csv'];
    writetable(measure_data, filename);
    fprintf('%s descriptive statistics saved to: %s\n', measure_name, filename);
    
    % Display interaction matrix for each measure
    interaction_data = measure_data(strcmp(measure_data.Effect_Type, 'Interaction'), :);
    if height(interaction_data) >= 18
        fprintf('\n--- %s: Interaction Means ---\n', measure_name);
        fprintf('%-15s %-10s %-10s %-10s %-10s %-10s %-10s\n', 'Movement\\Task', 'A-Single', 'A-Repeat', 'A-Switch', 'V-Single', 'V-Repeat', 'V-Switch');
        fprintf('%-15s %-10s %-10s %-10s %-10s %-10s %-10s\n', '-------------', '--------', '--------', '--------', '--------', '--------', '--------');
        
        % Reshape into 3x6 matrix (movement x modality/task combination)
        means_matrix = reshape([interaction_data.Mean], 6, 3)'; % 3 movements x 6 combinations
        
        for move = 1:3
            if strcmp(measure_name, 'RT')
                fprintf('%-15s %-10.1f %-10.1f %-10.1f %-10.1f %-10.1f %-10.1f\n', move_labels{move}, ...
                        means_matrix(move, 1), means_matrix(move, 2), means_matrix(move, 3), ...
                        means_matrix(move, 4), means_matrix(move, 5), means_matrix(move, 6));
            else % Accuracy
                fprintf('%-15s %-10.2f %-10.2f %-10.2f %-10.2f %-10.2f %-10.2f\n', move_labels{move}, ...
                        means_matrix(move, 1), means_matrix(move, 2), means_matrix(move, 3), ...
                        means_matrix(move, 4), means_matrix(move, 5), means_matrix(move, 6));
            end
        end
    end
end

% Create summary table with main effects only
main_effects_rt_accu = rt_accu_descriptive_table(~strcmp(rt_accu_descriptive_table.Effect_Type, 'Interaction'), :);
writetable(main_effects_rt_accu, [PATH_STAT '/beha_data/rt_accu_main_effects_summary.csv']);
fprintf('\nRT/Accuracy main effects summary saved to: %s\n', [PATH_STAT '/beha_data/rt_accu_main_effects_summary.csv']);

%% Create combined publication table
publication_rt_accu = table();
for meas = 1:length(measures)
    measure_data = rt_accu_descriptive_table(strcmp(rt_accu_descriptive_table.Measure, measures{meas}), :);
    
    % Add Mean_SD column for publication
    if strcmp(measures{meas}, 'RT')
        measure_data.Mean_SD = arrayfun(@(m,s) sprintf('%.1f (%.1f)', m, s), measure_data.Mean, measure_data.SD, 'UniformOutput', false);
    else
        measure_data.Mean_SD = arrayfun(@(m,s) sprintf('%.2f (%.2f)', m, s), measure_data.Mean, measure_data.SD, 'UniformOutput', false);
    end
    
    publication_rt_accu = [publication_rt_accu; measure_data];
end

writetable(publication_rt_accu, [PATH_STAT '/beha_data/rt_accu_publication_table.csv']);
fprintf('RT/Accuracy publication table saved to: %s\n', [PATH_STAT '/beha_data/rt_accu_publication_table.csv']);

fprintf('\nRT/Accuracy descriptive statistics completed.\n');


%% Function to generate comprehensive descriptive statistics for RT/Accuracy data
function [rt_accu_descriptive_table] = generate_rt_accu_descriptive_stats(anovaReacMat, anovaAccuMat, PATH_STAT)
% Generate comprehensive descriptive statistics table for RT and Accuracy data
% 
% Inputs:
%   anovaReacMat, anovaAccuMat - RT/Accuracy data matrices (subjects × 18 conditions)
%   PATH_STAT - path to save output files
%
% Output:
%   rt_accu_descriptive_table - comprehensive table with means/SDs for all conditions

fprintf('Generating RT/Accuracy descriptive statistics...\n');

%% Define experimental conditions for RT/Accuracy data
% 18 conditions: 3 MOVE × 2 MODA × 3 TASK
move_labels = {'stand', 'walk', 'pert'};
moda_labels = {'auditory', 'visual'};
task_labels = {'repeat-only', 'mixed-repeat', 'mixed-switch'};

% Based on your script comments:
% 1:6 - stand auditory single, stand auditory repeat, stand auditory switch, stand visual single, stand visual repeat, stand visual switch
% 7:12 - walk auditory single, walk auditory repeat, walk auditory switch, walk visual single, walk visual repeat, walk visual switch
% 13:18 - pert auditory single, pert auditory repeat, pert auditory switch, pert visual single, pert visual repeat, pert visual switch

interaction_labels = cell(1,18);
for move = 1:3
    for moda = 1:2
        for task = 1:3
            cond_idx = (move-1)*6 + (moda-1)*3 + task;
            interaction_labels{cond_idx} = sprintf('%s-%s-%s', ...
                move_labels{move}, moda_labels{moda}, task_labels{task});
        end
    end
end

%% RT and Accuracy data
rt_accu_data = struct();
rt_accu_data.RT = anovaReacMat;
rt_accu_data.Accuracy = anovaAccuMat * 100; % Convert to percentage

measures = fieldnames(rt_accu_data);
all_results = [];

%% Process each measure (RT, Accuracy)
for meas = 1:length(measures)
    measure_name = measures{meas};
    data_matrix = rt_accu_data.(measure_name);
    
    fprintf('Processing %s...\n', measure_name);
    
    %% 1. Individual condition statistics (18 conditions)
    for cond = 1:18
        if cond <= size(data_matrix, 2)
            result_row = struct();
            result_row.Measure = measure_name;
            result_row.Effect_Type = 'Interaction';
            result_row.Condition = interaction_labels{cond};
            
            condition_data = data_matrix(:, cond);
            valid_data = condition_data(~isnan(condition_data));
            
            if ~isempty(valid_data)
                result_row.Mean = mean(valid_data);
                result_row.SD = std(valid_data);
                result_row.N = length(valid_data);
                result_row.SEM = result_row.SD / sqrt(result_row.N);
                
                % Calculate factor levels
                move_idx = ceil(cond/6);
                within_move = mod(cond-1, 6) + 1;
                moda_idx = ceil(within_move/3);
                task_idx = mod(within_move-1, 3) + 1;
                
                result_row.Motor_Level = move_labels{move_idx};
                result_row.Modality_Level = moda_labels{moda_idx};
                result_row.Cognitive_Level = task_labels{task_idx};
            else
                result_row.Mean = NaN;
                result_row.SD = NaN;
                result_row.N = 0;
                result_row.SEM = NaN;
                result_row.Motor_Level = '';
                result_row.Modality_Level = '';
                result_row.Cognitive_Level = '';
            end
            
            all_results = [all_results; result_row];
        end
    end
    
    %% 2. Movement main effect (3 levels)
    for move = 1:3
        move_conditions = (move-1)*6 + (1:6); % 6 conditions per movement
        move_data = data_matrix(:, move_conditions);
        valid_data = move_data(~isnan(move_data));
        
        result_row = struct();
        result_row.Measure = measure_name;
        result_row.Effect_Type = 'Motor_Main_Effect';
        result_row.Condition = move_labels{move};
        result_row.Mean = mean(valid_data);
        result_row.SD = std(valid_data);
        result_row.N = length(valid_data);
        result_row.SEM = result_row.SD / sqrt(result_row.N);
        result_row.Motor_Level = move_labels{move};
        result_row.Modality_Level = 'All';
        result_row.Cognitive_Level = 'All';
        
        all_results = [all_results; result_row];
    end
    
    %% 3. Modality main effect (2 levels)
    % Auditory: 1:3, 7:9, 13:15
    % Visual: 4:6, 10:12, 16:18
    auditory_cols = [1:3, 7:9, 13:15];
    visual_cols = [4:6, 10:12, 16:18];
    
    for moda = 1:2
        if moda == 1
            moda_data = data_matrix(:, auditory_cols);
        else
            moda_data = data_matrix(:, visual_cols);
        end
        valid_data = moda_data(~isnan(moda_data));
        
        result_row = struct();
        result_row.Measure = measure_name;
        result_row.Effect_Type = 'Modality_Main_Effect';
        result_row.Condition = moda_labels{moda};
        result_row.Mean = mean(valid_data);
        result_row.SD = std(valid_data);
        result_row.N = length(valid_data);
        result_row.SEM = result_row.SD / sqrt(result_row.N);
        result_row.Motor_Level = 'All';
        result_row.Modality_Level = moda_labels{moda};
        result_row.Cognitive_Level = 'All';
        
        all_results = [all_results; result_row];
    end
    
    %% 4. Task difficulty main effect (3 levels)
    % Single: 1,4,7,10,13,16
    % Repeat: 2,5,8,11,14,17
    % Switch: 3,6,9,12,15,18
    single_cols = [1,4,7,10,13,16];
    repeat_cols = [2,5,8,11,14,17];
    switch_cols = [3,6,9,12,15,18];
    
    task_cols = {single_cols, repeat_cols, switch_cols};
    
    for task = 1:3
        task_data = data_matrix(:, task_cols{task});
        valid_data = task_data(~isnan(task_data));
        
        result_row = struct();
        result_row.Measure = measure_name;
        result_row.Effect_Type = 'Cognitive_Main_Effect';
        result_row.Condition = task_labels{task};
        result_row.Mean = mean(valid_data);
        result_row.SD = std(valid_data);
        result_row.N = length(valid_data);
        result_row.SEM = result_row.SD / sqrt(result_row.N);
        result_row.Motor_Level = 'All';
        result_row.Modality_Level = 'All';
        result_row.Cognitive_Level = task_labels{task};
        
        all_results = [all_results; result_row];
    end
end

%% Convert to table and save
rt_accu_descriptive_table = struct2table(all_results);
rt_accu_descriptive_table = rt_accu_descriptive_table(:, {'Measure', 'Effect_Type', 'Condition', ...
                                                          'Motor_Level', 'Modality_Level', 'Cognitive_Level', ...
                                                          'Mean', 'SD', 'SEM', 'N'});

% Save comprehensive table
writetable(rt_accu_descriptive_table, [PATH_STAT '/beha_data/rt_accu_descriptive_statistics.csv']);
fprintf('RT/Accuracy descriptive statistics saved to: %s\n', [PATH_STAT '/beha_data/rt_accu_descriptive_statistics.csv']);

% Display summary
fprintf('\n=== RT/ACCURACY DESCRIPTIVE STATISTICS SUMMARY ===\n');
disp(rt_accu_descriptive_table);

end

%%
function create_overall_behavioral_summary(PATH_STAT)
% Create an overall summary combining all behavioral measures

fprintf('\n=== CREATING OVERALL BEHAVIORAL SUMMARY ===\n');

% Load all descriptive tables if they exist
summary_data = struct();

try
    nasa_table = readtable([PATH_STAT '/subj/nasa_descriptive_statistics.csv']);
    summary_data.NASA = nasa_table;
    fprintf('Loaded NASA descriptive statistics\n');
catch
    fprintf('NASA descriptive statistics not found\n');
end

try
    rt_accu_table = readtable([PATH_STAT '/beha_data/rt_accu_descriptive_statistics.csv']);
    summary_data.RT_Accuracy = rt_accu_table;
    fprintf('Loaded RT/Accuracy descriptive statistics\n');
catch
    fprintf('RT/Accuracy descriptive statistics not found\n');
end

try
    step_table = readtable([PATH_STAT '/walk/step_descriptive_statistics.csv']);
    summary_data.Step = step_table;
    fprintf('Loaded Step descriptive statistics\n');
catch
    fprintf('Step descriptive statistics not found\n');
end

% Create combined main effects summary
if ~isempty(fieldnames(summary_data))
    combined_main_effects = table();
    
    fields = fieldnames(summary_data);
    for f = 1:length(fields)
        field_name = fields{f};
        data_table = summary_data.(field_name);
        
        % Extract main effects only
        main_effects = data_table(~strcmp(data_table.Effect_Type, 'Interaction'), :);
        
        % Add domain column
        if strcmp(field_name, 'NASA')
            main_effects.Domain = repmat({'Subjective'}, height(main_effects), 1);
        elseif strcmp(field_name, 'RT_Accuracy')
            main_effects.Domain = repmat({'Behavioral'}, height(main_effects), 1);
        elseif strcmp(field_name, 'Step')
            main_effects.Domain = repmat({'Motor'}, height(main_effects), 1);
        end
        
        combined_main_effects = [combined_main_effects; main_effects];
    end
    
    % Save combined summary
    writetable(combined_main_effects, [PATH_STAT '/combined_behavioral_main_effects_summary.csv']);
    fprintf('Combined behavioral main effects summary saved to: %s\n', [PATH_STAT '/combined_behavioral_main_effects_summary.csv']);
    
    % Display summary by domain and effect type
    fprintf('\n=== BEHAVIORAL MAIN EFFECTS SUMMARY ===\n');
    domains = unique(combined_main_effects.Domain);
    
    for d = 1:length(domains)
        domain_data = combined_main_effects(strcmp(combined_main_effects.Domain, domains{d}), :);
        fprintf('\n--- %s Domain ---\n', domains{d});
        
        effect_types = unique(domain_data.Effect_Type);
        for e = 1:length(effect_types)
            effect_data = domain_data(strcmp(domain_data.Effect_Type, effect_types{e}), :);
            fprintf('\n%s:\n', effect_types{e});
            
            if strcmp(domains{d}, 'Subjective')
                measures = unique(effect_data.Subscale);
                for m = 1:length(measures)
                    measure_data = effect_data(strcmp(effect_data.Subscale, measures{m}), :);
                    fprintf('  %s: %d conditions\n', measures{m}, height(measure_data));
                end
            elseif strcmp(domains{d}, 'Behavioral')
                measures = unique(effect_data.Measure);
                for m = 1:length(measures)
                    measure_data = effect_data(strcmp(effect_data.Measure, measures{m}), :);
                    fprintf('  %s: %d conditions\n', measures{m}, height(measure_data));
                end
            elseif strcmp(domains{d}, 'Motor')
                measures = unique(effect_data.Measure);
                for m = 1:length(measures)
                    measure_data = effect_data(strcmp(effect_data.Measure, measures{m}), :);
                    fprintf('  %s: %d conditions\n', measures{m}, height(measure_data));
                end
            end
        end
    end
end

fprintf('\nOverall behavioral summary completed.\n');

end

