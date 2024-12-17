 clear all
close all

%% Define PATHS and LISTS
mac_switch = 2; % 1: work; 2: home

if ispc
    PATH = 'D:\Data\grailswitch\';
elseif ismac
    if mac_switch == 1
        PATH = '/Users/julianreiser/owncloud/grailswitch';
        PATH_OUT = '/Volumes/Work2TB/Seafile/grailswitch';
    elseif mac_switch == 2
        PATH = '/Users/julianreiser/owncloud/grailswitch';
        PATH_OUT = '/Users/julianreiser/Seafile/grailswitch';
    end % if mac_switch
end

PATH_STAT = [PATH '/lists/stats'];
PATH_NASA = [PATH '/lists/subjective'];
PATH_PLOT = [PATH_OUT '/plots'];
subjlist = {'02','03','04','05','06','07','09','10','11','12','13','15','16','18','19','20','21','22','23','24','26','27'};

%% LOAD DATA
nasarawmat = xlsread([PATH_NASA '/gs_NASATLX.xlsx']);

%% Rearrange data in matrix
for subj = 1:length(subjlist)
    SUBJ = subjlist{subj};
    subjnum = str2num(SUBJ);

    nasasubidx = find(nasarawmat(:,1) == subjnum);
    nasatmpmat = nasarawmat(nasasubidx,26:end);
    condcount = 1;
    % get everything into submatrices
    for cond = 1:6:length(nasatmpmat)

        % calculate unweighted mean value
        nasamat(subj,condcount,:) = mean(nasatmpmat(cond : cond + 5),'omitnan');
        nasaeasymat(subj,condcount) = mean(nasatmpmat(cond : cond + 5),'omitnan');

        % calculate separate subdimension scores
        nasamatcogn(subj,condcount) = mean(nasatmpmat(cond),'omitnan');
        nasamatphys(subj,condcount) = mean(nasatmpmat(cond + 1),'omitnan');
        nasamattime(subj,condcount) = mean(nasatmpmat(cond + 2),'omitnan');
        nasamatperf(subj,condcount) = mean(nasatmpmat(cond + 3),'omitnan');
        nasamateffo(subj,condcount) = mean(nasatmpmat(cond + 4),'omitnan');
        nasamatfrus(subj,condcount) = mean(nasatmpmat(cond + 5),'omitnan');                 
        condcount = condcount + 1;

    end % for cond
end % for subj

% rearrange mean nasa conds
nasamattmp(:,1:2) = nasamat(:,1:2); % st v s, st v m
nasamattmp(:,3:4) = nasamat(:,7:8); % st a s, st a m
nasamattmp(:,5:6) = nasamat(:,3:4); % wa v s, wa v m
nasamattmp(:,7:8) = nasamat(:,9:10); % wa a s, wa a m
nasamattmp(:,9:10) = nasamat(:,5:6); % wa v s, wa v m
nasamattmp(:,11:12) = nasamat(:,11:12); % wa a s, wa a m
nasamat = nasamattmp; % write rearrangement back

% rearrange cogn conds
nasamattmpcogn(:,1:2) = nasamatcogn(:,1:2); % st v s, st v m
nasamattmpcogn(:,3:4) = nasamatcogn(:,7:8); % st a s, st a m
nasamattmpcogn(:,5:6) = nasamatcogn(:,3:4); % wa v s, wa v m
nasamattmpcogn(:,7:8) = nasamatcogn(:,9:10); % wa a s, wa a m
nasamattmpcogn(:,9:10) = nasamatcogn(:,5:6); % wa v s, wa v m
nasamattmpcogn(:,11:12) = nasamatcogn(:,11:12); % wa a s, wa a m
nasamatcogn = nasamattmpcogn; % write rearrangement back

% rearrange phys conds
nasamattmpphys(:,1:2) = nasamatphys(:,1:2); % st v s, st v m
nasamattmpphys(:,3:4) = nasamatphys(:,7:8); % st a s, st a m
nasamattmpphys(:,5:6) = nasamatphys(:,3:4); % wa v s, wa v m
nasamattmpphys(:,7:8) = nasamatphys(:,9:10); % wa a s, wa a m
nasamattmpphys(:,9:10) = nasamatphys(:,5:6); % wa v s, wa v m
nasamattmpphys(:,11:12) = nasamatphys(:,11:12); % wa a s, wa a m
nasamatphys = nasamattmpphys; % write rearrangement back

% rearrange time conds
nasamattmptime(:,1:2) = nasamattime(:,1:2); % st v s, st v m
nasamattmptime(:,3:4) = nasamattime(:,7:8); % st a s, st a m
nasamattmptime(:,5:6) = nasamattime(:,3:4); % wa v s, wa v m
nasamattmptime(:,7:8) = nasamattime(:,9:10); % wa a s, wa a m
nasamattmptime(:,9:10) = nasamattime(:,5:6); % wa v s, wa v m
nasamattmptime(:,11:12) = nasamattime(:,11:12); % wa a s, wa a m
nasamattime = nasamattmptime; % write rearrangement back

% rearrange perf conds
nasamattmpperf(:,1:2) = nasamatperf(:,1:2); % st v s, st v m
nasamattmpperf(:,3:4) = nasamatperf(:,7:8); % st a s, st a m
nasamattmpperf(:,5:6) = nasamatperf(:,3:4); % wa v s, wa v m
nasamattmpperf(:,7:8) = nasamatperf(:,9:10); % wa a s, wa a m
nasamattmpperf(:,9:10) = nasamatperf(:,5:6); % wa v s, wa v m
nasamattmpperf(:,11:12) = nasamatperf(:,11:12); % wa a s, wa a m
nasamatperf = nasamattmpperf; % write rearrangement back

% rearrange effo conds
nasamattmpeffo(:,1:2) = nasamateffo(:,1:2); % st v s, st v m
nasamattmpeffo(:,3:4) = nasamateffo(:,7:8); % st a s, st a m
nasamattmpeffo(:,5:6) = nasamateffo(:,3:4); % wa v s, wa v m
nasamattmpeffo(:,7:8) = nasamateffo(:,9:10); % wa a s, wa a m
nasamattmpeffo(:,9:10) = nasamateffo(:,5:6); % wa v s, wa v m
nasamattmpeffo(:,11:12) = nasamateffo(:,11:12); % wa a s, wa a m
nasamateffo = nasamattmpeffo; % write rearrangement back

% rearrange frus conds
nasamattmpfrus(:,1:2) = nasamatfrus(:,1:2); % st v s, st v m
nasamattmpfrus(:,3:4) = nasamatfrus(:,7:8); % st a s, st a m
nasamattmpfrus(:,5:6) = nasamatfrus(:,3:4); % wa v s, wa v m
nasamattmpfrus(:,7:8) = nasamatfrus(:,9:10); % wa a s, wa a m
nasamattmpfrus(:,9:10) = nasamatfrus(:,5:6); % wa v s, wa v m
nasamattmpfrus(:,11:12) = nasamatfrus(:,11:12); % wa a s, wa a m
nasamatfrus = nasamattmpfrus; % write rearrangement back

%% plot
% plot response times
anovaNASAMeans = mean(nasamat,1);
anovaNASASEMs = std(nasamat,1) / sqrt(size(nasamat,1));

cognNASAMeans = mean(nasamatcogn,1);
cognNASASEMs = std(nasamatcogn,1) / sqrt(size(nasamatcogn,1));

physNASAMeans = mean(nasamatphys,1);
physNASASEMs = std(nasamatphys,1) / sqrt(size(nasamatphys,1));

timeNASAMeans = mean(nasamattime,1);
timeNASASEMs = std(nasamattime,1) / sqrt(size(nasamattime,1));

perfNASAMeans = mean(nasamatperf,1);
perfNASASEMs = std(nasamatperf,1) / sqrt(size(nasamatperf,1));

effoNASAMeans = mean(nasamateffo,1);
effoNASASEMs = std(nasamateffo,1) / sqrt(size(nasamateffo,1));

frusNASAMeans = mean(nasamatfrus,1);
frusNASASEMs = std(nasamatfrus,1) / sqrt(size(nasamatfrus,1));

figure;
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 20, 20], 'PaperUnits', 'centimeters', 'PaperSize', [29.7,21])

lb_niceBARplot([anovaNASAMeans(1:4:end);anovaNASAMeans(3:4:end);anovaNASAMeans(2:4:end);anovaNASAMeans(4:4:end)], ...
               [anovaNASASEMs(1:4:end);anovaNASASEMs(3:4:end);anovaNASASEMs(2:4:end);anovaNASASEMs(4:4:end)], ...
               'All NASA scales','Mean Raw Score',{'visual - single task','auditory - single task','visual - mixed task','auditory - mixed task'},[0,20])

saveas(gcf,[PATH_PLOT '/subjective/nasa_means.png']); close gcf;

figure;
    subplot(321)
    lb_niceBARplot([cognNASAMeans(1:4);cognNASAMeans(5:8);cognNASAMeans(9:12)], ...
                   [cognNASASEMs(1:4);cognNASASEMs(5:8);cognNASASEMs(9:12)], ...
                    'NASA cognitive','ms',{'stand','walk','pert'},[0,20],'label_legend',{'v_single','v_mixed','a_single','a_mixed'})

    subplot(322)
    lb_niceBARplot([physNASAMeans(1:4);physNASAMeans(5:8);physNASAMeans(9:12)], ...
                   [physNASASEMs(1:4);physNASASEMs(5:8);physNASASEMs(9:12)], ...
                   'NASA physical','ms',{'stand','walk','pert'},[0,20])

    subplot(323)
    lb_niceBARplot([timeNASAMeans(1:4);timeNASAMeans(5:8);timeNASAMeans(9:12)], ...
                   [timeNASASEMs(1:4);timeNASASEMs(5:8);timeNASASEMs(9:12)], ...
                   'NASA time','ms',{'stand','walk','pert'},[0,20])

    subplot(324)
    lb_niceBARplot([perfNASAMeans(1:4);perfNASAMeans(5:8);perfNASAMeans(9:12)], ...
                   [perfNASASEMs(1:4);perfNASASEMs(5:8);perfNASASEMs(9:12)], ...
                    'NASA performance','ms',{'stand','walk','pert'},[0,20])

    subplot(325)
    lb_niceBARplot([effoNASAMeans(1:4);effoNASAMeans(5:8);effoNASAMeans(9:12)], ...
                   [effoNASASEMs(1:4);effoNASASEMs(5:8);effoNASASEMs(9:12)], ...
                   'NASA effort','ms',{'stand','walk','pert'},[0,20])
    subplot(326)
    lb_niceBARplot([frusNASAMeans(1:4);frusNASAMeans(5:8);frusNASAMeans(9:12)], ...
                   [frusNASASEMs(1:4);frusNASASEMs(5:8);frusNASASEMs(9:12)], ...
                   'NASA frustration','ms',{'stand','walk','pert'},[0,20])

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 35, 20], 'PaperUnits', 'centimeters', 'PaperSize', [29.7,21])
saveas(gcf,[PATH_PLOT '/subjective/nasa_subscales.png']); close gcf;

%% ---------------------- now for the rmANOVAs ---------------------------
anova_walklist = 1:3;
anova_modalist = 1:2;
anova_mixlist = 1:2;

multcomptype = 'lsd';
full_length = length(anova_walklist) * length(anova_modalist) * length(anova_mixlist);

% Create a table reflecting the within subject factors
MOVE = cell(full_length,1); % lat ring conditions: 1, 5
MODA = cell(full_length,1); % lat move conditions: stand, walk, pert
MIX = cell(full_length,1); % lat ecce conditions: 0, 18, 40

% Assiging the values to the parameters based on the data sorting
c1 = cell(1,1); c1{1} = 'st'; c1 = repmat(c1,full_length/3,1); MOVE(1:4,1) = c1;
c1 = cell(1,1); c1{1} = 'wa'; c1 = repmat(c1,full_length/3,1); MOVE(5:8,1) = c1;
c1 = cell(1,1); c1{1} = 'pe'; c1 = repmat(c1,full_length/3,1); MOVE(9:12,1) = c1;

c1 = cell(1,1); c1{1} = 'v'; c1 = repmat(c1,full_length/2,1); MODA([1:2,5:6,9:10],1) = c1;
c1 = cell(1,1); c1{1} = 'a'; c1 = repmat(c1,full_length/2,1); MODA([3:4,7:8,11:12],1) = c1;

c1 = cell(1,1); c1{1} = 'singl'; c1 = repmat(c1,full_length/2,1); MIX([1:2:end],1) = c1;
c1 = cell(1,1); c1{1} = 'mixed'; c1 = repmat(c1,full_length/2,1); MIX([2:2:end],1) = c1;

% Create the within table
factorNames = {'MOVE','MODA','MIX'};
within = table(MOVE,MODA,MIX, 'VariableNames', factorNames);

% now create ANOVA table
varNames = cell(full_length,1);
for i = 1 : full_length
    v = strcat('V',num2str(i));
    varNames{i,1} = v;
end% Create a table storing the respones
anovaNASATable = array2table(nasamat, 'VariableNames',varNames);
cognNASATable = array2table(nasamatcogn, 'VariableNames',varNames);
physNASATable = array2table(nasamatphys, 'VariableNames',varNames);
timeNASATable = array2table(nasamattime, 'VariableNames',varNames);
perfNASATable = array2table(nasamatperf, 'VariableNames',varNames);
effoNASATable = array2table(nasamateffo, 'VariableNames',varNames);
frusNASATable = array2table(nasamatfrus, 'VariableNames',varNames);

%% compute the rmANOVA
clear rmStat

% complete nasa - fit the repeated measures model
nasa_rm = fitrm(anovaNASATable,'V1-V12~1','WithinDesign',within);
[rmStat.nasa.rmanova] = ranova(nasa_rm, 'WithinModel','MOVE*MODA*MIX');
%[FCzN2_amp_move FCzN2_amp_ecce FCzN2_amp_int rmStat.FCzN2.ampMOVE rmStat.FCzN2.ampECCE rmStat.FCzN2.ampINT FCzN2_amp_movestat FCzN2_amp_eccestat FCzN2_amp_intstat] = lw_phtest(n1amps1,within,2);
lw_qq(nasamat,4,3,nasa_rm.mauchly)
rmStat.nasa.rmanova = lw_fdr(rmStat.nasa.rmanova,'corrMethod',1);
rmStat.nasa.rmanova = lw_eta(rmStat.nasa.rmanova);
rmStat.nasa.mauchly = nasa_rm.mauchly;
rmStat.nasa.phMove = multcompare(nasa_rm,'MOVE');
rmStat.nasa.phMix = multcompare(nasa_rm,'MIX');
rmStat.nasa.phMoveByModa = multcompare(nasa_rm,'MOVE','By','MODA');
rmStat.nasa.phModaByMove = multcompare(nasa_rm,'MODA','By','MOVE');
rmStat.nasa.phMixByModa = multcompare(nasa_rm,'MIX','By','MODA');
rmStat.nasa.phModaByMix = multcompare(nasa_rm,'MODA','By','MIX');

MainMove= [mean(nasamat(:,[1:4]),[1,2]),std(nasamat(:,[1:4]),0,[1,2]);
           mean(nasamat(:,[5:8]),[1,2]),std(nasamat(:,[5:8]),0,[1,2]);
           mean(nasamat(:,[9:12]),[1,2]),std(nasamat(:,[9:12]),0,[1,2])];

MainDiff = [mean(nasamat(:,[1:2:end]),[1,2]),std(nasamat(:,[1:2:end]),0,[1,2]);
            mean(nasamat(:,[2:2:end]),[1,2]),std(nasamat(:,[2:2:end]),0,[1,2])];


% cogn nasa - fit the repeated measures model
cogn_rm = fitrm(cognNASATable,'V1-V12~1','WithinDesign',within);
[rmStat.cogn.rmanova] = ranova(cogn_rm, 'WithinModel','MOVE*MODA*MIX');
%[FCzN2_amp_move FCzN2_amp_ecce FCzN2_amp_int rmStat.FCzN2.ampMOVE rmStat.FCzN2.ampECCE rmStat.FCzN2.ampINT FCzN2_amp_movestat FCzN2_amp_eccestat FCzN2_amp_intstat] = lw_phtest(n1amps1,within,2);
lw_qq(nasamatcogn,4,3,cogn_rm.mauchly)
rmStat.cogn.rmanova = lw_fdr(rmStat.cogn.rmanova,'corrMethod',1);
rmStat.cogn.rmanova = lw_eta(rmStat.cogn.rmanova);
rmStat.cogn.mauchly = cogn_rm.mauchly;
rmStat.cogn.phMoveByModa = multcompare(cogn_rm,'MOVE','By','MODA');
rmStat.cogn.phModaByMove = multcompare(cogn_rm,'MODA','By','MOVE');
rmStat.cogn.phMixByModa = multcompare(cogn_rm,'MIX','By','MODA');

% phys nasa - fit the repeated measures model
phys_rm = fitrm(physNASATable,'V1-V12~1','WithinDesign',within);
[rmStat.phys.rmanova] = ranova(phys_rm, 'WithinModel','MOVE*MODA*MIX');
%[FCzN2_amp_move FCzN2_amp_ecce FCzN2_amp_int rmStat.FCzN2.ampMOVE rmStat.FCzN2.ampECCE rmStat.FCzN2.ampINT FCzN2_amp_movestat FCzN2_amp_eccestat FCzN2_amp_intstat] = lw_phtest(n1amps1,within,2);
lw_qq(nasamatphys,4,3,phys_rm.mauchly)
rmStat.phys.rmanova = lw_fdr(rmStat.phys.rmanova,'corrMethod',1);
rmStat.phys.rmanova = lw_eta(rmStat.phys.rmanova);
rmStat.phys.mauchly = phys_rm.mauchly;
rmStat.phys.phMoveByModa = multcompare(phys_rm,'MOVE','By','MODA');
rmStat.phys.phModaByMove = multcompare(phys_rm,'MODA','By','MOVE');
rmStat.phys.phMixByModa = multcompare(phys_rm,'MIX','By','MODA');

% time nasa - fit the repeated measures model
time_rm = fitrm(timeNASATable,'V1-V12~1','WithinDesign',within);
[rmStat.time.rmanova] = ranova(time_rm, 'WithinModel','MOVE*MODA*MIX');
%[FCzN2_amp_move FCzN2_amp_ecce FCzN2_amp_int rmStat.FCzN2.ampMOVE rmStat.FCzN2.ampECCE rmStat.FCzN2.ampINT FCzN2_amp_movestat FCzN2_amp_eccestat FCzN2_amp_intstat] = lw_phtest(n1amps1,within,2);
lw_qq(nasamattime,4,3,time_rm.mauchly)
rmStat.time.rmanova = lw_fdr(rmStat.time.rmanova,'corrMethod',1);
rmStat.time.rmanova = lw_eta(rmStat.time.rmanova);
rmStat.time.mauchly = time_rm.mauchly;
rmStat.time.phMoveByModa = multcompare(time_rm,'MOVE','By','MODA');
rmStat.time.phModaByMove = multcompare(time_rm,'MODA','By','MOVE');
rmStat.time.phMixByModa = multcompare(time_rm,'MIX','By','MODA');

% perf nasa - fit the repeated measures model
perf_rm = fitrm(perfNASATable,'V1-V12~1','WithinDesign',within);
[rmStat.perf.rmanova] = ranova(perf_rm, 'WithinModel','MOVE*MODA*MIX');
%[FCzN2_amp_move FCzN2_amp_ecce FCzN2_amp_int rmStat.FCzN2.ampMOVE rmStat.FCzN2.ampECCE rmStat.FCzN2.ampINT FCzN2_amp_movestat FCzN2_amp_eccestat FCzN2_amp_intstat] = lw_phtest(n1amps1,within,2);
lw_qq(nasamatperf,4,3,perf_rm.mauchly)
rmStat.perf.rmanova = lw_fdr(rmStat.perf.rmanova,'corrMethod',1);
rmStat.perf.rmanova = lw_eta(rmStat.perf.rmanova);
rmStat.perf.mauchly = perf_rm.mauchly;
rmStat.perf.phMoveByModa = multcompare(perf_rm,'MOVE','By','MODA');
rmStat.perf.phModaByMove = multcompare(perf_rm,'MODA','By','MOVE');
rmStat.perf.phMixByModa = multcompare(perf_rm,'MIX','By','MODA');

% effo nasa - fit the repeated measures model
effo_rm = fitrm(effoNASATable,'V1-V12~1','WithinDesign',within);
[rmStat.effo.rmanova] = ranova(effo_rm, 'WithinModel','MOVE*MODA*MIX');
%[FCzN2_amp_move FCzN2_amp_ecce FCzN2_amp_int rmStat.FCzN2.ampMOVE rmStat.FCzN2.ampECCE rmStat.FCzN2.ampINT FCzN2_amp_movestat FCzN2_amp_eccestat FCzN2_amp_intstat] = lw_phtest(n1amps1,within,2);
lw_qq(nasamateffo,4,3,effo_rm.mauchly)
rmStat.effo.rmanova = lw_fdr(rmStat.effo.rmanova,'corrMethod',1);
rmStat.effo.rmanova = lw_eta(rmStat.effo.rmanova);
rmStat.effo.mauchly = effo_rm.mauchly;
rmStat.effo.phMoveByModa = multcompare(effo_rm,'MOVE','By','MODA');
rmStat.effo.phModaByMove = multcompare(effo_rm,'MODA','By','MOVE');
rmStat.effo.phMixByModa = multcompare(effo_rm,'MIX','By','MODA');

% frus nasa - fit the repeated measures model
frus_rm = fitrm(frusNASATable,'V1-V12~1','WithinDesign',within);
[rmStat.frus.rmanova] = ranova(frus_rm, 'WithinModel','MOVE*MODA*MIX');
%[FCzN2_amp_move FCzN2_amp_ecce FCzN2_amp_int rmStat.FCzN2.ampMOVE rmStat.FCzN2.ampECCE rmStat.FCzN2.ampINT FCzN2_amp_movestat FCzN2_amp_eccestat FCzN2_amp_intstat] = lw_phtest(n1amps1,within,2);
lw_qq(nasamatfrus,4,3,frus_rm.mauchly)
rmStat.frus.rmanova = lw_fdr(rmStat.frus.rmanova,'corrMethod',1);
rmStat.frus.rmanova = lw_eta(rmStat.frus.rmanova);
rmStat.frus.mauchly = frus_rm.mauchly;
rmStat.frus.phMoveByModa = multcompare(frus_rm,'MOVE','By','MODA');
rmStat.frus.phModaByMove = multcompare(frus_rm,'MODA','By','MOVE');
rmStat.frus.phMixByModa = multcompare(frus_rm,'MIX','By','MODA');

close all
%% save this
save([PATH_STAT '/subj/nasa_all.mat'],'nasamat')
save([PATH_STAT '/subj/nasa_cogn.mat'],'nasamatcogn')
save([PATH_STAT '/subj/nasa_phys.mat'],'nasamatphys')
save([PATH_STAT '/subj/nasa_time.mat'],'nasamattime')
save([PATH_STAT '/subj/nasa_perf.mat'],'nasamatperf')
save([PATH_STAT '/subj/nasa_effo.mat'],'nasamateffo')
save([PATH_STAT '/subj/nasa_frus.mat'],'nasamatfrus')