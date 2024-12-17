%% calculate response related statistics
clear all

%% load paths
mac_switch = 2; % 1: work, 2: home

if ispc
    PATH = 'D:\Data\grailswitch\';
elseif ismac
    if mac_switch == 1
        PATH = '/Volumes/Work2TB/Seafile/grailswitch/';
    elseif mac_switch == 2
        PATH = '/Users/julianreiser/Seafile/grailswitch/';
    end
end

PATH_STAT = [PATH '/stats_rework'];
PATH_PLOT = [PATH '/plots/beha'];

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
audiReacMat(:,1) = mean(squeeze(tmpmatr_a(:,1:2,1)),2,'omitnan');
audiReacMat(:,2) = mean(squeeze(tmpmatr_a(:,5:6,1)),2,'omitnan');
audiReacMat(:,3) = mean(squeeze(tmpmatr_a(:,7:8,1)),2,'omitnan');
audiReacMat(:,4) = mean(squeeze(tmpmatr_a(:,1:2,2)),2,'omitnan');
audiReacMat(:,5) = mean(squeeze(tmpmatr_a(:,5:6,2)),2,'omitnan');
audiReacMat(:,6) = mean(squeeze(tmpmatr_a(:,7:8,2)),2,'omitnan');
audiReacMat(:,7) = mean(squeeze(tmpmatr_a(:,1:2,3)),2,'omitnan');
audiReacMat(:,8) = mean(squeeze(tmpmatr_a(:,5:6,3)),2,'omitnan');
audiReacMat(:,9) = mean(squeeze(tmpmatr_a(:,7:8,3)),2,'omitnan');

visuReacMat(:,1) = mean(squeeze(tmpmatr_v(:,1:2,1)),2,'omitnan');
visuReacMat(:,2) = mean(squeeze(tmpmatr_v(:,5:6,1)),2,'omitnan');
visuReacMat(:,3) = mean(squeeze(tmpmatr_v(:,7:8,1)),2,'omitnan');
visuReacMat(:,4) = mean(squeeze(tmpmatr_v(:,1:2,2)),2,'omitnan');
visuReacMat(:,5) = mean(squeeze(tmpmatr_v(:,5:6,2)),2,'omitnan');
visuReacMat(:,6) = mean(squeeze(tmpmatr_v(:,7:8,2)),2,'omitnan');
visuReacMat(:,7) = mean(squeeze(tmpmatr_v(:,1:2,3)),2,'omitnan');
visuReacMat(:,8) = mean(squeeze(tmpmatr_v(:,5:6,3)),2,'omitnan');
visuReacMat(:,9) = mean(squeeze(tmpmatr_v(:,7:8,3)),2,'omitnan');

% accuracy
audiAccuMat(:,1) = mean(squeeze(tmpmat_a(:,1:2,1)),2,'omitnan');
audiAccuMat(:,2) = mean(squeeze(tmpmat_a(:,5:6,1)),2,'omitnan');
audiAccuMat(:,3) = mean(squeeze(tmpmat_a(:,7:8,1)),2,'omitnan');
audiAccuMat(:,4) = mean(squeeze(tmpmat_a(:,1:2,2)),2,'omitnan');
audiAccuMat(:,5) = mean(squeeze(tmpmat_a(:,5:6,2)),2,'omitnan');
audiAccuMat(:,6) = mean(squeeze(tmpmat_a(:,7:8,2)),2,'omitnan');
audiAccuMat(:,7) = mean(squeeze(tmpmat_a(:,1:2,3)),2,'omitnan');
audiAccuMat(:,8) = mean(squeeze(tmpmat_a(:,5:6,3)),2,'omitnan');
audiAccuMat(:,9) = mean(squeeze(tmpmat_a(:,7:8,3)),2,'omitnan');

visuAccuMat(:,1) = mean(squeeze(tmpmat_v(:,1:2,1)),2,'omitnan');
visuAccuMat(:,2) = mean(squeeze(tmpmat_v(:,5:6,1)),2,'omitnan');
visuAccuMat(:,3) = mean(squeeze(tmpmat_v(:,7:8,1)),2,'omitnan');
visuAccuMat(:,4) = mean(squeeze(tmpmat_v(:,1:2,2)),2,'omitnan');
visuAccuMat(:,5) = mean(squeeze(tmpmat_v(:,5:6,2)),2,'omitnan');
visuAccuMat(:,6) = mean(squeeze(tmpmat_v(:,7:8,2)),2,'omitnan');
visuAccuMat(:,7) = mean(squeeze(tmpmat_v(:,1:2,3)),2,'omitnan');
visuAccuMat(:,8) = mean(squeeze(tmpmat_v(:,5:6,3)),2,'omitnan');
visuAccuMat(:,9) = mean(squeeze(tmpmat_v(:,7:8,3)),2,'omitnan');

% prepare anova description
anovaTitles = {'st_single','st_repeat','st_switch',...
               'wa_single','wa_repeat','wa_switch',...
               'pa_single','pa_repeat','pa_switch'};

% plot response times
audiReacMeans = mean(audiReacMat,1);
audiReacSEMs = std(audiReacMat,1) / sqrt(size(audiReacMat,1));

audiAccuMeans = mean(audiAccuMat,1)*100;
audiAccuSEMs = std(audiAccuMat,1) / sqrt(size(audiAccuMat,1))*100;

visuReacMeans = mean(visuReacMat,1);
visuReacSEMs = std(visuReacMat,1) / sqrt(size(visuReacMat,1));

visuAccuMeans = mean(visuAccuMat,1)*100;
visuAccuSEMs = std(visuAccuMat,1) / sqrt(size(visuAccuMat,1))*100;

%% ---------------------- now for the rmANOVAs ---------------------------
anova_walklist = 1:3;
anova_tasklist = 1:3;

multcomptype = 'lsd';
full_length = length(anova_walklist) * length(anova_tasklist);

% Create a table reflecting the within subject factors
MOVE = cell(full_length,1); % lat ring conditions: 1, 5
TASK = cell(full_length,1); % lat ecce conditions: 0, 18, 40

% Assiging the values to the parameters based on the data sorting
c1 = cell(1,1); c1{1} = 'st'; c1 = repmat(c1,full_length/3,1); MOVE(1:3,1) = c1;
c1 = cell(1,1); c1{1} = 'wa'; c1 = repmat(c1,full_length/3,1); MOVE(4:6,1) = c1;
c1 = cell(1,1); c1{1} = 'pe'; c1 = repmat(c1,full_length/3,1); MOVE(7:9,1) = c1;

c1 = cell(1,1); c1{1} = 'single'; c1 = repmat(c1,full_length/3,1); TASK([1:3:end],1) = c1;
c1 = cell(1,1); c1{1} = 'repead'; c1 = repmat(c1,full_length/3,1); TASK([2:3:end],1) = c1;
c1 = cell(1,1); c1{1} = 'switch'; c1 = repmat(c1,full_length/3,1); TASK([3:3:end],1) = c1;

% Create the within table
factorNames = {'MOVE','TASK'};
within = table(MOVE,TASK, 'VariableNames', factorNames);

% now create ANOVA table
varNames = cell(full_length,1);
for i = 1 : full_length
    v = strcat('V',num2str(i));
    varNames{i,1} = v;
end% Create a table storing the respones
audiReacTable = array2table(audiReacMat, 'VariableNames',varNames);
audiAccuTable = array2table(audiAccuMat, 'VariableNames',varNames);

visuReacTable = array2table(visuReacMat, 'VariableNames',varNames);
visuAccuTable = array2table(visuAccuMat, 'VariableNames',varNames);

%% compute the rmANOVA
clear rmStat

% Response times - fit the repeated measures model
audireac_rm = fitrm(audiReacTable,'V1-V9~1','WithinDesign',within);
[rmStat.audireac.rmanova] = ranova(audireac_rm, 'WithinModel','MOVE*TASK');
lw_qq(audiReacMat,3,3,audireac_rm.mauchly)
rmStat.audireac.rmanova = lw_fdr(rmStat.audireac.rmanova,'corrMethod',1);
rmStat.audireac.rmanova = lw_eta(rmStat.audireac.rmanova);
rmStat.audireac.mauchly = audireac_rm.mauchly;
rmStat.audireac.phTask = multcompare(audireac_rm,'TASK');
rmStat.audireac.phMoveByTask = multcompare(audireac_rm,'MOVE','By','TASK');

% Response times - fit the repeated measures model
visureac_rm = fitrm(visuReacTable,'V1-V9~1','WithinDesign',within);
[rmStat.visureac.rmanova] = ranova(visureac_rm, 'WithinModel','MOVE*TASK');
lw_qq(visuReacMat,3,3,visureac_rm.mauchly)
rmStat.visureac.rmanova = lw_fdr(rmStat.visureac.rmanova,'corrMethod',1);
rmStat.visureac.rmanova = lw_eta(rmStat.visureac.rmanova);
rmStat.visureac.mauchly = visureac_rm.mauchly;
rmStat.visureac.phTask = multcompare(visureac_rm,'TASK');
rmStat.visureac.phMoveByTask = multcompare(visureac_rm,'MOVE','By','TASK');


% Accuracy - fit the repeated measures model
audiaccu_rm = fitrm(audiAccuTable,'V1-V9~1','WithinDesign',within);
[rmStat.audiaccu.rmanova] = ranova(audiaccu_rm, 'WithinModel','MOVE*TASK');
lw_qq(audiAccuMat,3,3,audiaccu_rm.mauchly)
rmStat.audiaccu.rmanova = lw_fdr(rmStat.audiaccu.rmanova,'corrMethod',1);
rmStat.audiaccu.rmanova = lw_eta(rmStat.audiaccu.rmanova);
rmStat.audiaccu.mauchly = audiaccu_rm.mauchly;
rmStat.audiaccu.phTask = multcompare(audiaccu_rm,'TASK');
rmStat.audiaccu.phMoveByTask = multcompare(audiaccu_rm,'MOVE','By','TASK');

% Accuracy - fit the repeated measures model
visuaccu_rm = fitrm(visuAccuTable,'V1-V9~1','WithinDesign',within);
[rmStat.visuaccu.rmanova] = ranova(visuaccu_rm, 'WithinModel','MOVE*TASK');
lw_qq(visuAccuMat,3,3,visuaccu_rm.mauchly)
rmStat.visuaccu.rmanova = lw_fdr(rmStat.visuaccu.rmanova,'corrMethod',1);
rmStat.visuaccu.rmanova = lw_eta(rmStat.visuaccu.rmanova);
rmStat.visuaccu.mauchly = visuaccu_rm.mauchly;
rmStat.visuaccu.phTask = multcompare(visuaccu_rm,'TASK');
rmStat.visuaccu.phMoveByTask = multcompare(visuaccu_rm,'MOVE','By','TASK');
%% now plot
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
