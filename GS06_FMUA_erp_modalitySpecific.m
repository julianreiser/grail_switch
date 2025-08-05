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
        PATH = '/Volumes/Work4TB/Seafile/grailswitch';
    elseif mac_switch == 2
        PATH = '/Users/julianreiser/Seafile/grailswitch';
    end % mac_switch
end

PATH_RAW = [PATH '/data/RAW/'];
PATH_ICA = [PATH '/data/processed/ICA'];
PATH_EEG = [PATH '/data/processed/EEG'];
PATH_PRUNED =[PATH '/data/processed/PRUNED'];

% LISTS / STATS
PATH_LIST = [PATH '/lists'];
PATH_COND = [PATH_LIST '/conditions'];
PATH_STAT = [PATH '/stats_rework'];
PATH_PLOT = [PATH '/plots'];
subjlist = {'02','03','04','05','06','07','09','10','11','12','13','15','16','18','19','20','21','22','23','24','25','27'};

%% Load GND_grailswitch
GND = load([PATH_STAT '/FMUT/GND_grailswitch_cue_audi.mat']);
GND_cue_audi = GND.GND;

GND = load([PATH_STAT '/FMUT/GND_grailswitch_cue_visu.mat']);
GND_cue_visu = GND.GND;

GND = load([PATH_STAT '/FMUT/GND_grailswitch_tar_audi.mat']);
GND_tar_audi = GND.GND;

GND = load([PATH_STAT '/FMUT/GND_grailswitch_tar_visu.mat']);
GND_tar_visu = GND.GND;

% correct wrong times-scaling
if GND_cue_audi.time_pts(end) < 1000
    GND_cue_audi.time_pts = GND_cue_audi.time_pts * 1000;
    GND_cue_visu.time_pts = GND_cue_visu.time_pts * 1000;
    GND_tar_audi.time_pts = GND_tar_audi.time_pts * 1000;
    GND_tar_visu.time_pts = GND_tar_visu.time_pts * 1000;
end

%% set GND parameters
%Define some variables
nPerm = 1e4;
chanHood = 75;

% P1 params
P1StimWindow = [80 160];
P1Chans = {'POz','PO3','PO4','PO7','PO8','PO9','PO10','Oz','O1','O2'};

% CNV params
CNVStimWindow = [-150 50];
CNVPlotWindow = [-800 1000];
CNVChans = {'F1','F2','Fz','FC1','FC2','Cz','C1','C2'};
CNVylim = [-8,4];

% N2 params
N2StimWindow = [250 400];
N2Chans = {'F1','F2','Fz','FC1','FC2','Cz','C1','C2'};

% P3 params
P3StimWindow = [300 600];
P3PlotWindow = [-800 1000];
P3Chans = {'CPz','CP1','CP2','Pz','P1','P2','P3','P4','POz','PO3','PO4'};
P3ylim = [-4,6];

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
        PATH = '/Volumes/Work4TB/Seafile/grailswitch';
    elseif mac_switch == 2
        PATH = '/Users/julianreiser/Seafile/grailswitch';
    end % mac_switch
end

PATH_RAW = [PATH '/data/RAW/'];
PATH_ICA = [PATH '/data/processed/ICA'];
PATH_EEG = [PATH '/data/processed/EEG'];
PATH_PRUNED =[PATH '/data/processed/PRUNED'];

% LISTS / STATS
PATH_LIST = [PATH '/lists'];
PATH_COND = [PATH_LIST '/conditions'];
PATH_STAT = [PATH '/stats_rework'];
PATH_PLOT = [PATH '/plots'];
subjlist = {'02','03','04','05','06','07','09','10','11','12','13','15','16','18','19','20','21','22','23','24','25','27'};

%% Load GND_grailswitch
GND = load([PATH_STAT '/FMUT/GND_grailswitch_cue_audi.mat']);
GND_cue_audi = GND.GND;

GND = load([PATH_STAT '/FMUT/GND_grailswitch_cue_visu.mat']);
GND_cue_visu = GND.GND;

GND = load([PATH_STAT '/FMUT/GND_grailswitch_tar_audi.mat']);
GND_tar_audi = GND.GND;

GND = load([PATH_STAT '/FMUT/GND_grailswitch_tar_visu.mat']);
GND_tar_visu = GND.GND;

% correct wrong times-scaling
if GND_cue_audi.time_pts(end) < 1000
    GND_cue_audi.time_pts = GND_cue_audi.time_pts * 1000;
    GND_cue_visu.time_pts = GND_cue_visu.time_pts * 1000;
    GND_tar_audi.time_pts = GND_tar_audi.time_pts * 1000;
    GND_tar_visu.time_pts = GND_tar_visu.time_pts * 1000;
end

%% set GND parameters
%Define some variables
nPerm = 1e4;
chanHood = 75;

% P1 params
P1StimWindow = [80 160];
P1Chans = {'POz','PO3','PO4','PO7','PO8','PO9','PO10','Oz','O1','O2'};

% CNV params
CNVStimWindow = [-150 50];
CNVChans = {'F1','F2','Fz','FC1','FC2','Cz','C1','C2'};

% N2 params
N2StimWindow = [250 400];
N2Chans = {'F1','F2','Fz','FC1','FC2','Cz','C1','C2'};

% P3 params
P3StimWindow = [300 600];
P3Chans = {'CPz','CP1','CP2','Pz','P1','P2','P3','P4','POz','PO3','PO4'};

%% Function to extract comprehensive statistics from F_tests (handling multiple clusters)
%% Function to extract comprehensive statistics from F_tests (handling multiple clusters)
function [stats_table] = extract_statistics(F_test_struct, effect_names)
    stats_table = table();
    
    for i = 1:length(effect_names)
        effect_name = effect_names{i};
        
        if isfield(F_test_struct.df, effect_name) && isfield(F_test_struct.clust_info, effect_name)
            % Extract degrees of freedom
            if isfield(F_test_struct.df, 'error')
                df_effect = F_test_struct.df.(effect_name);
                df_error = F_test_struct.df.error;
            else
                % Alternative structure - extract from effect field
                df_vals = F_test_struct.df.(effect_name);
                if length(df_vals) >= 2
                    df_effect = df_vals(1);
                    df_error = df_vals(2);
                else
                    df_effect = df_vals(1);
                    df_error = NaN; % Will need to be filled manually
                end
            end
            
            % Check cluster information
            clust_ids = F_test_struct.clust_info.(effect_name).clust_ids;
            p_val_raw = F_test_struct.clust_info.(effect_name).pval;
            
            % Determine number of clusters
            unique_cluster_ids = unique(clust_ids(:));
            unique_cluster_ids = unique_cluster_ids(unique_cluster_ids > 0); % Remove 0 (non-cluster)
            
            if isempty(unique_cluster_ids) % No significant clusters
                % Create row for non-significant effect
                effect_label = effect_name;
                if isempty(p_val_raw)
                    new_row = table({effect_label}, df_effect, df_error, ...
                    NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
                    NaN, 0, 0, {[NaN, NaN]}, {[]}, {'n.s.'}, ...
                    'VariableNames', {'Effect', 'df_effect', 'df_error', ...
                    'F_max', 'F_min', 'F_mean', ...
                    'eta_sq_max', 'eta_sq_min', 'eta_sq_mean', 'eta_sq_cluster_mean', 'eta_sq_cluster_std', ...
                    'p_value', 'cluster_number', 'cluster_size', 'time_range_ms', 'cluster_channels', 'Status'});
                else
                new_row = table({effect_label}, df_effect, df_error, ...
                    NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
                    p_val_raw(1), 0, 0, {[NaN, NaN]}, {[]}, {'n.s.'}, ...
                    'VariableNames', {'Effect', 'df_effect', 'df_error', ...
                    'F_max', 'F_min', 'F_mean', ...
                    'eta_sq_max', 'eta_sq_min', 'eta_sq_mean', 'eta_sq_cluster_mean', 'eta_sq_cluster_std', ...
                    'p_value', 'cluster_number', 'cluster_size', 'time_range_ms', 'cluster_channels', 'Status'});
                end

                if isempty(stats_table)
                    stats_table = new_row;
                else
                    stats_table = [stats_table; new_row];
                end
                
            else % Process each significant cluster separately
                for cluster_num = 1:length(unique_cluster_ids)
                    cluster_id = unique_cluster_ids(cluster_num);
                    
                    % Get p-value for this cluster (if multiple clusters exist)
                    if length(p_val_raw) >= cluster_num
                        p_val = p_val_raw(cluster_num);
                    else
                        p_val = p_val_raw(1); % Fallback to first p-value
                    end
                    
                    % Create effect label for this cluster
                    if length(unique_cluster_ids) > 1
                        effect_label = sprintf('%s_cluster%d', effect_name, cluster_num);
                    else
                        effect_label = effect_name;
                    end
                    
                    % Extract F values only from this specific cluster
                    if isfield(F_test_struct.F_obs, effect_name)
                        F_obs_matrix = F_test_struct.F_obs.(effect_name);
                        cluster_mask = clust_ids == cluster_id;
                        cluster_F_values = F_obs_matrix(cluster_mask);
                        
                        % Calculate statistics for all points in this cluster
                        F_max = max(cluster_F_values);
                        F_min = min(cluster_F_values);
                        F_mean = mean(cluster_F_values);
                        
                        % Calculate partial eta squared for max, min, and mean F values
                        eta_sq_max = (F_max * df_effect) / (F_max * df_effect + df_error);
                        eta_sq_min = (F_min * df_effect) / (F_min * df_effect + df_error);
                        eta_sq_mean = (F_mean * df_effect) / (F_mean * df_effect + df_error);
                        
                        % Calculate eta squared for ALL points in this cluster
                        eta_sq_all_points = (cluster_F_values * df_effect) ./ (cluster_F_values * df_effect + df_error);
                        eta_sq_cluster_mean = mean(eta_sq_all_points);
                        eta_sq_cluster_std = std(eta_sq_all_points);
                        
                        % Count cluster size (number of significant time-electrode points)
                        cluster_size = sum(cluster_mask(:));
                        
                        % Find cluster extent information
                        [cluster_chans, cluster_times] = find(cluster_mask);
                        if ~isempty(cluster_chans)
                            % Get actual time points and channel names
                            if isfield(F_test_struct, 'used_tpt_ids')
                                actual_times = F_test_struct.used_tpt_ids(unique(cluster_times));
                                time_range = [min(actual_times), max(actual_times)];
                            else
                                time_range = [min(cluster_times), max(cluster_times)];
                            end
                            
                            if isfield(F_test_struct, 'used_chan_ids')
                                actual_chans = F_test_struct.used_chan_ids(unique(cluster_chans));
                            else
                                actual_chans = unique(cluster_chans);
                            end
                            
                            % Convert time indices back to actual time values if needed
                            if isfield(F_test_struct, 'time_wind')
                                % Time points are likely in the time_wind range
                                time_start = F_test_struct.time_wind(1) + (min(cluster_times)-1) * 2; % Assuming 2ms resolution
                                time_end = F_test_struct.time_wind(1) + (max(cluster_times)-1) * 2;
                                time_range_ms = [time_start, time_end];
                            else
                                time_range_ms = time_range;
                            end
                        else
                            time_range_ms = [NaN, NaN];
                            actual_chans = [];
                        end
                        
                        % Store comprehensive statistics in table
                        new_row = table({effect_label}, df_effect, df_error, ...
                            F_max, F_min, F_mean, ...
                            eta_sq_max, eta_sq_min, eta_sq_mean, eta_sq_cluster_mean, eta_sq_cluster_std, ...
                            p_val, cluster_num, cluster_size, {time_range_ms}, {actual_chans}, {'significant'}, ...
                            'VariableNames', {'Effect', 'df_effect', 'df_error', ...
                            'F_max', 'F_min', 'F_mean', ...
                            'eta_sq_max', 'eta_sq_min', 'eta_sq_mean', 'eta_sq_cluster_mean', 'eta_sq_cluster_std', ...
                            'p_value', 'cluster_number', 'cluster_size', 'time_range_ms', 'cluster_channels', 'Status'});
                    else
                        % No F_obs field
                        new_row = table({effect_label}, df_effect, df_error, ...
                            NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
                            p_val, cluster_num, 0, {[NaN, NaN]}, {[]}, {'n.s.'}, ...
                            'VariableNames', {'Effect', 'df_effect', 'df_error', ...
                            'F_max', 'F_min', 'F_mean', ...
                            'eta_sq_max', 'eta_sq_min', 'eta_sq_mean', 'eta_sq_cluster_mean', 'eta_sq_cluster_std', ...
                            'p_value', 'cluster_number', 'cluster_size', 'time_range_ms', 'cluster_channels', 'Status'});
                    end
                    
                    % Concatenate tables properly
                    if isempty(stats_table)
                        stats_table = new_row;
                    else
                        stats_table = [stats_table; new_row];
                    end
                end
            end
        end
    end
end

%% Function to create cluster visualization with eta squared values (handles multiple clusters)
function create_cluster_eta_plot(F_test_struct, effect_name, save_path, component_name)
    if ~isfield(F_test_struct.clust_info, effect_name)
        fprintf('No cluster info found for effect: %s\n', effect_name);
        return; % No cluster info to plot
    end
    
    clust_ids = F_test_struct.clust_info.(effect_name).clust_ids;
    unique_cluster_ids = unique(clust_ids(:));
    unique_cluster_ids = unique_cluster_ids(unique_cluster_ids > 0); % Remove 0 (non-cluster)
    
    if isempty(unique_cluster_ids)
        fprintf('No significant clusters found for effect: %s\n', effect_name);
        return; % No significant clusters to plot
    end
    
    fprintf('Creating cluster plots for %s - %s (%d clusters)\n', component_name, effect_name, length(unique_cluster_ids));
    
    % Extract cluster information
    F_obs_matrix = F_test_struct.F_obs.(effect_name);
    if isfield(F_test_struct.df, 'error')
        df_effect = F_test_struct.df.(effect_name);
        df_error = F_test_struct.df.error;
    else
        df_vals = F_test_struct.df.(effect_name);
        if length(df_vals) >= 2
            df_effect = df_vals(1);
            df_error = df_vals(2);
        else
            df_effect = df_vals(1);
            df_error = 42; % Default error df for typical ERP studies
        end
    end
    
    % Calculate eta squared for all points
    eta_sq_matrix = (F_obs_matrix * df_effect) ./ (F_obs_matrix * df_effect + df_error);
    
    % Create plot for each cluster or combined plot
    if length(unique_cluster_ids) == 1
        % Single cluster - create one plot
        cluster_mask = clust_ids == unique_cluster_ids(1);
        eta_sq_cluster = eta_sq_matrix;
        eta_sq_cluster(~cluster_mask) = NaN;
        
        create_single_cluster_plot(eta_sq_cluster, eta_sq_matrix, cluster_mask, ...
            F_test_struct, effect_name, component_name, save_path, 1);
        
    else
        % Multiple clusters - create separate plots for each + combined plot
        
        % Combined plot
        eta_sq_all_clusters = eta_sq_matrix;
        eta_sq_all_clusters(clust_ids == 0) = NaN;
        
        [path_dir, filename, ext] = fileparts(save_path);
        combined_path = fullfile(path_dir, [filename '_all_clusters' ext]);
        
        figure;
        set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 25, 15], 'PaperUnits', 'centimeters', 'PaperSize', [26,16])
        
        imagesc(eta_sq_all_clusters);
        
        % Use unipolar colormap (MATLAB 2025a compatible)
        try
            colormap(viridis); % Modern unipolar colormap
        catch
            try
                colormap(parula); % Default unipolar colormap
            catch
                colormap(hot); % Fallback unipolar colormap
            end
        end
        
        % Create colorbar with eta squared label
        cb = colorbar;
        ylabel(cb, '\eta^2_{partial}', 'FontSize', 14, 'Rotation', 270, 'VerticalAlignment', 'bottom');
        caxis([0, .5]);
        
        % Add contour lines for all cluster boundaries
        hold on;
        colors = {'k', 'w', 'r', 'g', 'b', 'm', 'c'}; % Different colors for different clusters
        for c = 1:length(unique_cluster_ids)
            cluster_id = unique_cluster_ids(c);
            color_idx = mod(c-1, length(colors)) + 1;
            contour(clust_ids == cluster_id, [0.5 0.5], colors{color_idx}, 'LineWidth', 3);
        end
        
        % Simplified title without eta squared
        title([component_name ' - ' effect_name ' - All Clusters'], 'FontSize', 16);
        xlabel('Time Points', 'FontSize', 14);
        ylabel('Channels', 'FontSize', 14);
        
        add_axis_labels(F_test_struct, eta_sq_all_clusters);
        set(gcf,'color','w');
        saveas(gcf, combined_path);
        close(gcf);
        
        % Individual cluster plots
        for c = 1:length(unique_cluster_ids)
            cluster_id = unique_cluster_ids(c);
            cluster_mask = clust_ids == cluster_id;
            eta_sq_cluster = eta_sq_matrix;
            eta_sq_cluster(~cluster_mask) = NaN;
            
            individual_path = fullfile(path_dir, [filename '_cluster' num2str(c) ext]);
            create_single_cluster_plot(eta_sq_cluster, eta_sq_matrix, cluster_mask, ...
                F_test_struct, effect_name, component_name, individual_path, cluster_id);
        end
    end
end

function create_single_cluster_plot(eta_sq_cluster, eta_sq_matrix, cluster_mask, ...
    F_test_struct, effect_name, component_name, save_path, cluster_num)
    
    figure;
    set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 25, 15], 'PaperUnits', 'centimeters', 'PaperSize', [26,16])
    
    % Plot eta squared values with cluster outline
    imagesc(eta_sq_cluster);
    
    % Use unipolar colormap (MATLAB 2025a compatible)
    try
        colormap(viridis); % Modern unipolar colormap
    catch
        try
            colormap(parula); % Default unipolar colormap
        catch
            colormap(hot); % Fallback unipolar colormap
        end
    end
    
    % Create colorbar with eta squared label
    cb = colorbar;
    ylabel(cb, '\eta^2_{partial}', 'FontSize', 14, 'Rotation', 270, 'VerticalAlignment', 'bottom');
    caxis([0, .5]);
    
    % Add contour lines for cluster boundaries
    hold on;
    contour(cluster_mask, [0.5 0.5], 'k-', 'LineWidth', 3);
    
    % Simplified title without eta squared
    if cluster_num == 1 && sum(cluster_mask(:)) == sum(~isnan(eta_sq_cluster(:)))
        title_str = [component_name ' - ' effect_name ' - Significant Cluster'];
    else
        title_str = [component_name ' - ' effect_name ' - Cluster ' num2str(cluster_num)];
    end
    title(title_str, 'FontSize', 16);
    xlabel('Time Points', 'FontSize', 14);
    ylabel('Channels', 'FontSize', 14);
    
    add_axis_labels(F_test_struct, eta_sq_cluster);
    set(gcf,'color','w');
    
    % Save the plot
    saveas(gcf, save_path);
    close(gcf);
end

function add_axis_labels(F_test_struct, plot_data)
    % Add channel labels if available
    if isfield(F_test_struct, 'include_chans')
        yticks(1:length(F_test_struct.include_chans));
        yticklabels(F_test_struct.include_chans);
    end
    
    % Add time information if available
    if isfield(F_test_struct, 'time_wind')
        time_start = F_test_struct.time_wind(1);
        time_end = F_test_struct.time_wind(2);
        n_timepoints = size(plot_data, 2);
        if n_timepoints > 1
            time_ticks = linspace(1, n_timepoints, 5);
            time_labels = linspace(time_start, time_end, 5);
            xticks(time_ticks);
            xticklabels(arrayfun(@(x) sprintf('%.0f', x), time_labels, 'UniformOutput', false));
            xlabel('Time (ms)', 'FontSize', 14);
        end
    end
end

%% Function to extract post-hoc comparisons for significant clusters
function [posthoc_table] = extract_cluster_posthoc(GND_struct, F_test_struct, effect_name, factor_names, factor_levels)
    posthoc_table = table();
    
    if ~isfield(F_test_struct.clust_info, effect_name)
        return;
    end
    
    clust_ids = F_test_struct.clust_info.(effect_name).clust_ids;
    unique_cluster_ids = unique(clust_ids(:));
    unique_cluster_ids = unique_cluster_ids(unique_cluster_ids > 0);
    
    if isempty(unique_cluster_ids)
        return;
    end
    
    % Factor level labels
    if length(factor_names) == 1
        if strcmp(factor_names{1}, 'mixed')
            level_labels = {'repeat-only', 'mixed - repeat', 'mixed - switch'};
            bin_grouping = {[1,4,7], [2,5,8], [3,6,9]}; % single, repeat, switch
        elseif strcmp(factor_names{1}, 'move')
            level_labels = {'stand', 'walk', 'perturbation'};
            bin_grouping = {[1,2,3], [4,5,6], [7,8,9]}; % stand, walk, perturbation
        end
    elseif length(factor_names) == 2
        % Interaction: mixed x move
        level_labels = {'Single-Stand', 'Single-Walk', 'Single-Pert', ...
                       'Repeat-Stand', 'Repeat-Walk', 'Repeat-Pert', ...
                       'Switch-Stand', 'Switch-Walk', 'Switch-Pert'};
        bin_grouping = {1, 2, 3, 4, 5, 6, 7, 8, 9}; % each bin separately
    end
    
    % define channels
    for chani = 1:length(GND_struct.F_tests.include_chans)
        chan_mask(chani) = find(strcmpi({GND_struct.chanlocs.labels},GND_struct.F_tests.include_chans{chani}));
    end

    % define time-window
    time_mask = find(ismember(GND_struct.time_pts,GND_struct.F_tests.time_wind))

    % Process each cluster
    for cluster_num = 1:length(unique_cluster_ids)
        cluster_id = unique_cluster_ids(cluster_num);
        cluster_mask = clust_ids == cluster_id;
        
        if ~any(cluster_mask(:))
            continue;
        end
        
        % Extract cluster data (average across cluster time-channel points)
        cluster_data = [];
        for subj = 1:size(GND_struct.indiv_erps, 4)
            subj_data = squeeze(GND_struct.indiv_erps(chan_mask, time_mask(1):time_mask(2), :, subj));
            cluster_mean = mean(subj_data(cluster_mask), 'all'); % Average across cluster points
            cluster_data(subj, :) = cluster_mean;
        end
        
        % Organize data by factor levels
        level_means = [];
        level_stds = [];
        level_data = {};
        
        for level = 1:length(level_labels)
            if length(factor_names) == 1
                % Main effect: average across bins for this level
                level_bins = bin_grouping{level};
                level_values = [];
                for bin = level_bins
                    if bin <= size(GND_struct.indiv_erps, 3)
                        bin_data = squeeze(GND_struct.indiv_erps(chan_mask, time_mask(1):time_mask(2), bin, :));
                        for subj = 1:size(bin_data, 3)
                            subj_bin_data = bin_data(:, :, subj);
                            level_values(end+1) = mean(subj_bin_data(cluster_mask), 'all');
                        end
                    end
                end
            else
                % Interaction: specific bin
                bin = bin_grouping{level};
                if bin <= size(GND_struct.indiv_erps, 3)
                    bin_data = squeeze(GND_struct.indiv_erps(chan_mask, time_mask(1):time_mask(2), bin, :));
                    level_values = [];
                    for subj = 1:size(bin_data, 3)
                        subj_bin_data = bin_data(:, :, subj);
                        level_values(end+1) = mean(subj_bin_data(cluster_mask), 'all');
                    end
                end
            end
            
            level_data{level} = level_values;
            level_means(level) = mean(level_values);
            level_stds(level) = std(level_values);
        end
        
        % Create pairwise comparisons
        for i = 1:length(level_labels)
            for j = (i+1):length(level_labels)
                mean_diff = level_means(i) - level_means(j);
                
                % Calculate pooled standard error for difference
                n1 = length(level_data{i});
                n2 = length(level_data{j});
                pooled_se = sqrt((level_stds(i)^2/n1) + (level_stds(j)^2/n2));
                
                % Calculate t-statistic and partial eta squared
                if pooled_se > 0
                    t_stat = mean_diff / pooled_se;
                    df = n1 + n2 - 2;
                    partial_eta_sq = (t_stat^2) / (t_stat^2 + df);
                    % two-tailed p-value from Student-t distribution
                    p_val = 2 * (1 - tcdf(abs(t_stat), df));
                else
                    t_stat = 0;
                    partial_eta_sq = 0;
                    p_val = NaN;                       % No variance → undefined p

                end
                
                % Create effect label
                if length(unique_cluster_ids) > 1
                    effect_label = sprintf('%s_cluster%d', effect_name, cluster_num);
                else
                    effect_label = effect_name;
                end
                
                % Add to table
                new_row = table({effect_label}, cluster_num, sum(cluster_mask(:)), ...
                    {level_labels{i}}, level_means(i), level_stds(i), ...
                    {level_labels{j}}, level_means(j), level_stds(j), ...
                    mean_diff, pooled_se, t_stat, df, p_val, partial_eta_sq, ...
                    'VariableNames', {'Effect','Cluster_Number','Cluster_Size', ...
                    'Level_A','Mean_A','SD_A','Level_B','Mean_B','SD_B', ...
                    'Mean_Diff','Pooled_SE','t_stat','df','p_value','Partial_Eta_Squared'});

                
                posthoc_table = [posthoc_table; new_row];
            end
        end
        
        % Also add descriptive statistics table
        fprintf('\n=== CLUSTER %d DESCRIPTIVE STATISTICS ===\n', cluster_num);
        fprintf('Effect: %s, Cluster Size: %d points\n', effect_name, sum(cluster_mask(:)));
        fprintf('%-15s %-10s %-10s %-10s\n', 'Level', 'Mean(µV)', 'SD(µV)', 'N');
        fprintf('%-15s %-10s %-10s %-10s\n', '-----', '--------', '------', '-');
        for level = 1:length(level_labels)
            fprintf('%-15s %-10.3f %-10.3f %-10d\n', level_labels{level}, ...
                level_means(level), level_stds(level), length(level_data{level}));
        end
    end
end
%% plotting
function create_cluster_mean_plots(GND_struct, F_test_struct, effect_name, save_path, component_name, spacing)
% Function to create cluster-based mean plots showing average amplitude 
% within significant clusters for each experimental condition

%% Default spacing if not provided
if nargin < 6 || isempty(spacing)
    spacing = 0.05;
end

%% Check for cluster information
if ~isfield(F_test_struct.clust_info, effect_name)
    fprintf('No cluster info found for effect: %s\n', effect_name);
    return;
end

clust_ids = F_test_struct.clust_info.(effect_name).clust_ids;
unique_cluster_ids = unique(clust_ids(:));
unique_cluster_ids = unique_cluster_ids(unique_cluster_ids > 0);

if isempty(unique_cluster_ids)
    fprintf('No significant clusters found for effect: %s\n', effect_name);
    return;
end

fprintf('Creating cluster mean plots for %s - %s (%d clusters)\n', component_name, effect_name, length(unique_cluster_ids));

%% Determine experimental design
factor_info = setup_experimental_design(effect_name);

%% Create plots for each cluster
for cluster_num = 1:length(unique_cluster_ids)
    cluster_id = unique_cluster_ids(cluster_num);
    cluster_mask = clust_ids == cluster_id;
    
    if ~any(cluster_mask(:))
        continue;
    end
    
    % Extract cluster data for this cluster
    [condition_data, condition_labels, condition_colors] = extract_cluster_condition_data(GND_struct, cluster_mask, factor_info);
    
    % Create individual cluster plot
    if length(unique_cluster_ids) == 1
        plot_save_path = save_path;
        plot_title = sprintf('%s - %s - Cluster Mean', component_name, effect_name);
    else
        [path_dir, filename, ext] = fileparts(save_path);
        plot_save_path = fullfile(path_dir, [filename '_cluster' num2str(cluster_num) ext]);
        plot_title = sprintf('%s - %s - Cluster %d Mean', component_name, effect_name, cluster_num);
    end
    
    create_single_cluster_mean_plot(condition_data, condition_labels, condition_colors, ...
        plot_title, plot_save_path, sum(cluster_mask(:)), spacing, factor_info);
end

%% Create combined plot if multiple clusters
if length(unique_cluster_ids) > 1
    [path_dir, filename, ext] = fileparts(save_path);
    combined_save_path = fullfile(path_dir, [filename '_all_clusters' ext]);
    
    create_multi_cluster_comparison_plot(GND_struct, F_test_struct, effect_name, ...
        component_name, combined_save_path, factor_info, spacing);
end

end

function factor_info = setup_experimental_design(effect_name)
% Setup experimental design information based on effect name

factor_info = struct();

if contains(effect_name, 'X') || contains(effect_name, 'interaction')
    % Interaction effect
    factor_info.type = 'interaction';
    factor_info.labels = {'Single-Stand', 'Single-Walk', 'Single-Pert', ...
                         'Repeat-Stand', 'Repeat-Walk', 'Repeat-Pert', ...
                         'Switch-Stand', 'Switch-Walk', 'Switch-Pert'};
    factor_info.colors = {'#2E2E2E', '#4A4A4A', '#666666', ...  % Single (dark grays)
                         '#8B1A73', '#B8357E', '#D655A0', ...   % Repeat (magentas) 
                         '#1B4F8C', '#235DB8', '#4A7BC8'};      % Switch (blues)
    factor_info.bin_grouping = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    factor_info.x_positions = 1:9;
    factor_info.xlabel = 'Mixed Task Condition';
    
elseif contains(effect_name, 'mixed')
    % Mixed task main effect
    factor_info.type = 'main_effect';
    factor_info.labels = {'repeat-only', 'mixed - repeat', 'mixed - switch'};
    factor_info.colors = {'#272730', '#B8357E', '#235DB8'};
    factor_info.bin_grouping = {[1,4,7], [2,5,8], [3,6,9]};
    factor_info.x_positions = 1:3;
    factor_info.xlabel = 'Mixed Task Condition';
    
elseif contains(effect_name, 'move')
    % Movement main effect
    factor_info.type = 'main_effect';
    factor_info.labels = {'stand', 'walk', 'perturbation'};
    factor_info.colors = {'#2D5016', '#4A7C22', '#73A942'};
    factor_info.bin_grouping = {[1,2,3], [4,5,6], [7,8,9]};
    factor_info.x_positions = 1:3;
    factor_info.xlabel = 'Movement Condition';
    
else
    % Default/unknown effect
    factor_info.type = 'unknown';
    factor_info.labels = {'Condition1', 'Condition2', 'Condition3'};
    factor_info.colors = {'#333333', '#666666', '#999999'};
    factor_info.bin_grouping = {1, 2, 3};
    factor_info.x_positions = 1:3;
    factor_info.xlabel = 'Condition';
end

end

function [condition_data, condition_labels, condition_colors] = extract_cluster_condition_data(GND_struct, cluster_mask, factor_info)
% Extract ERP data averaged within cluster for each experimental condition

condition_data = {};
condition_labels = factor_info.labels;
condition_colors = factor_info.colors;

for cond = 1:length(condition_labels)
    cond_values = [];

    % define channels
    for chani = 1:length(GND_struct.F_tests.include_chans)
        chan_mask(chani) = find(strcmpi({GND_struct.chanlocs.labels},GND_struct.F_tests.include_chans{chani}));
    end

    % define time-window
    time_mask = find(ismember(GND_struct.time_pts,GND_struct.F_tests.time_wind))
    
    if strcmp(factor_info.type, 'main_effect')
        % Main effect: average across multiple bins
        bins = factor_info.bin_grouping{cond};
        
        for bin = bins
            if bin <= size(GND_struct.indiv_erps, 3) && bin > 0
                for subj = 1:size(GND_struct.indiv_erps, 4)
                    subj_data = squeeze(GND_struct.indiv_erps(chan_mask, time_mask(1):time_mask(2), bin, subj));
                    cluster_mean = mean(subj_data(cluster_mask), 'all');
                    if ~isnan(cluster_mean) && ~isinf(cluster_mean)
                        cond_values(end+1) = cluster_mean;
                    end
                end
            end
        end
        
    else
        % Interaction: specific bin only
        bin = factor_info.bin_grouping{cond};
        
        if bin <= size(GND_struct.indiv_erps, 3) && bin > 0
            for subj = 1:size(GND_struct.indiv_erps, 4)
                subj_data = squeeze(GND_struct.indiv_erps(chan_mask, time_mask(1):time_mask(2), bin, subj));
                cluster_mean = mean(subj_data(cluster_mask), 'all');
                if ~isnan(cluster_mean) && ~isinf(cluster_mean)
                    cond_values(end+1) = cluster_mean;
                end
            end
        end
    end
    
    condition_data{cond} = cond_values;
end

end

function create_single_cluster_mean_plot(condition_data, condition_labels, condition_colors, plot_title, save_path, cluster_size, spacing, factor_info)
% Create mean plot for a single cluster

% Calculate means and SEMs
means = zeros(1, length(condition_data));
sems = zeros(1, length(condition_data));

for cond = 1:length(condition_data)
    if ~isempty(condition_data{cond})
        means(cond) = mean(condition_data{cond});
        sems(cond) = std(condition_data{cond}) / sqrt(length(condition_data{cond}));
    else
        means(cond) = NaN;
        sems(cond) = NaN;
    end
end

% Create figure
figure;
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 20, 15], ...
    'PaperUnits', 'centimeters', 'PaperSize', [21, 16]);

hold on;

x_positions = factor_info.x_positions;

if strcmp(factor_info.type, 'interaction')
    % For interaction plots, show mixed conditions on x-axis with lines for movement
    movement_groups = {[1,4,7], [2,5,8], [3,6,9]}; % Stand, Walk, Perturbation groups
    movement_colors = {'#272730', '#B8357E', '#235DB8'}; % Stand (black), Walk (magenta), Pert (blue)
    movement_labels = {'stand', 'walk', 'perturbation'};
    mixed_x_positions = 1:3; % Single, Repeat, Switch
    
    % Plot each movement condition as a connected line across mixed conditions
    for move = 1:3
        move_indices = movement_groups{move};
        move_means = means(move_indices);
        move_sems = sems(move_indices);
        
        % Only plot if we have valid data
        valid_idx = ~isnan(move_means);
        if any(valid_idx)
            errorbar(mixed_x_positions(valid_idx), move_means(valid_idx), move_sems(valid_idx), ...
                '-o', 'LineWidth', 2, 'MarkerSize', 8, ...
                'Color', movement_colors{move}, 'MarkerFaceColor', movement_colors{move}, ...
                'DisplayName', movement_labels{move});
        end
    end
    
    % Set x-axis labels for mixed conditions
    set(gca, 'XTick', mixed_x_positions, 'XTickLabel', {'repeat-only', 'mixed - repeat', 'mixed - switch'}, 'FontSize', 12);
    
    % Add legend
    legend('Location', 'best', 'FontSize', 12);
    
else
    % For main effects, use connected lines
    valid_idx = ~isnan(means);
    if any(valid_idx)
        errorbar(x_positions(valid_idx), means(valid_idx), sems(valid_idx), ...
            '-o', 'LineWidth', 2, 'MarkerSize', 8, ...
            'Color', condition_colors{1}, 'MarkerFaceColor', condition_colors{1});
    end
    
    set(gca, 'XTick', x_positions, 'XTickLabel', condition_labels);
end

hold off;

% Formatting
xlabel(factor_info.xlabel, 'FontSize', 14);
ylabel('Cluster Mean Amplitude (µV)', 'FontSize', 14);
title(sprintf('%s (n=%d points)', plot_title, cluster_size), 'FontSize', 16);
set(gca, 'FontSize', 12);
grid on;
set(gcf, 'Color', 'w');
set(gca, 'TickDir', 'out');

% Auto-scale y-axis with some padding
valid_means = means(~isnan(means));
valid_sems = sems(~isnan(sems));
if ~isempty(valid_means)
    y_range = [min(valid_means - valid_sems), max(valid_means + valid_sems)];
    y_padding = diff(y_range) * 0.1;
    ylim([y_range(1) - y_padding, y_range(2) + y_padding]);
end

% Save and close
saveas(gcf, save_path);
close(gcf);

% Print summary statistics
fprintf('\n=== CLUSTER MEAN SUMMARY ===\n');
fprintf('Cluster Size: %d time-channel points\n', cluster_size);
fprintf('%-20s %-12s %-12s %-8s\n', 'Condition', 'Mean(µV)', 'SEM(µV)', 'N');
fprintf('%-20s %-12s %-12s %-8s\n', '--------', '--------', '-------', '-');
for cond = 1:length(condition_labels)
    if ~isempty(condition_data{cond})
        fprintf('%-20s %-12.3f %-12.3f %-8d\n', condition_labels{cond}, ...
            means(cond), sems(cond), length(condition_data{cond}));
    end
end

end

function create_multi_cluster_comparison_plot(GND_struct, F_test_struct, effect_name, component_name, save_path, factor_info, spacing)
% Create comparison plot showing all clusters side by side

clust_ids = F_test_struct.clust_info.(effect_name).clust_ids;
unique_cluster_ids = unique(clust_ids(:));
unique_cluster_ids = unique_cluster_ids(unique_cluster_ids > 0);

% Create subplot figure
n_clusters = length(unique_cluster_ids);
figure;
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 7*n_clusters, 15], ...
    'PaperUnits', 'centimeters', 'PaperSize', [7*n_clusters+1, 16]);

% Track y-limits for consistent scaling
all_means = [];
all_sems = [];

subplot_data = {};
for cluster_num = 1:n_clusters
    cluster_id = unique_cluster_ids(cluster_num);
    cluster_mask = clust_ids == cluster_id;
    
    [condition_data, condition_labels, condition_colors] = extract_cluster_condition_data(GND_struct, cluster_mask, factor_info);
    
    % Calculate means and SEMs
    means = zeros(1, length(condition_data));
    sems = zeros(1, length(condition_data));
    
    for cond = 1:length(condition_data)
        if ~isempty(condition_data{cond})
            means(cond) = mean(condition_data{cond});
            sems(cond) = std(condition_data{cond}) / sqrt(length(condition_data{cond}));
        else
            means(cond) = NaN;
            sems(cond) = NaN;
        end
    end
    
    subplot_data{cluster_num} = struct('means', means, 'sems', sems, 'cluster_size', sum(cluster_mask(:)));
    valid_means = means(~isnan(means));
    valid_sems = sems(~isnan(sems));
    all_means = [all_means, valid_means];
    all_sems = [all_sems, valid_sems];
end

% Determine common y-limits
if ~isempty(all_means)
    y_range = [min(all_means - all_sems), max(all_means + all_sems)];
    y_padding = diff(y_range) * 0.1;
    common_ylim = [y_range(1) - y_padding, y_range(2) + y_padding];
else
    common_ylim = [-1, 1];
end

% Create subplots
for cluster_num = 1:n_clusters
    subplot(1, n_clusters, cluster_num);
    
    means = subplot_data{cluster_num}.means;
    sems = subplot_data{cluster_num}.sems;
    cluster_size = subplot_data{cluster_num}.cluster_size;
    
    hold on;
    
    x_positions = factor_info.x_positions;
    
    if strcmp(factor_info.type, 'interaction')
        % For interaction plots in subplots, show mixed conditions on x-axis
        movement_groups = {[1,4,7], [2,5,8], [3,6,9]}; % Stand, Walk, Perturbation groups
        movement_colors = {'#272730', '#B8357E', '#235DB8'};
        movement_labels = {'stand', 'walk', 'perturabtion'};
        mixed_x_positions = 1:3; % Single, Repeat, Switch
        
        % Plot each movement condition as a connected line
        for move = 1:3
            move_indices = movement_groups{move};
            move_means = means(move_indices);
            move_sems = sems(move_indices);
            
            valid_idx = ~isnan(move_means);
            if any(valid_idx)
                errorbar(mixed_x_positions(valid_idx), move_means(valid_idx), move_sems(valid_idx), ...
                    '-o', 'LineWidth', 2, 'MarkerSize', 6, ...
                    'Color', movement_colors{move}, 'MarkerFaceColor', movement_colors{move});
            end
        end
        
        % Set x-axis labels for mixed conditions
        set(gca, 'XTick', mixed_x_positions, 'XTickLabel', {'repeat-only', 'mixed - repeat', 'mixed - switch'}, 'FontSize', 10);
        
    else
        % For main effects in subplots
        valid_idx = ~isnan(means);
        if any(valid_idx)
            errorbar(x_positions(valid_idx), means(valid_idx), sems(valid_idx), ...
                '-o', 'LineWidth', 2, 'MarkerSize', 8, ...
                'Color', condition_colors{1}, 'MarkerFaceColor', condition_colors{1});
        end
        set(gca, 'XTick', x_positions, 'XTickLabel', condition_labels, 'FontSize', 10);
    end
    
    hold off;
    
    % Formatting
    if cluster_num == 1
        ylabel('Cluster Mean Amplitude (µV)', 'FontSize', 12);
    end
    xlabel(factor_info.xlabel, 'FontSize', 10);
    title(sprintf('Cluster %d\n(n=%d points)', cluster_num, cluster_size), 'FontSize', 14);
    ylim(common_ylim);
    grid on;
    set(gca, 'TickDir', 'out');
end

% Add overall title
sgtitle(sprintf('%s - %s - All Clusters Comparison', component_name, effect_name), 'FontSize', 16);

% Add legend for interaction effects to help identify the 9 conditions
if strcmp(factor_info.type, 'interaction') && n_clusters > 0
    % Create a single legend for the entire figure
    legend(condition_labels, 'Location', 'eastoutside', 'FontSize', 10);
end

set(gcf, 'Color', 'w');

% Save and close
saveas(gcf, save_path);
close(gcf);

end

%% Full factorial designs
%% Stim-related CNV
GND_cue_SCNV_audi = FclustGND(GND_cue_audi, ...
                    'bins', 1:9, ...
                    'factor_names', {'mixed', 'move'}, ...
                    'factor_levels', [3, 3], ...
                    'time_wind', CNVStimWindow, ...
                    'include_chans',CNVChans, ...
                    'chan_hood', chanHood, ...
                    'n_perm', nPerm, ...
                    'plot_raster','yes', ...
                    'save_GND', 'no');

save([PATH_STAT '/FMUT/GND_GS_StimCNV_fullFactorial_cue_audi.mat'],'GND_cue_SCNV_audi');

%% Extract statistics for CNV auditory
effect_names = {'mixed', 'move', 'mixedXmove'};
CNV_audi_stats = extract_statistics(GND_cue_SCNV_audi.F_tests, effect_names);
disp('CNV Auditory Statistics:');
disp(CNV_audi_stats);

% Create cluster eta squared plots for CNV auditory
for i = 1:length(effect_names)
    create_cluster_eta_plot(GND_cue_SCNV_audi.F_tests, effect_names{i}, ...
        [PATH_PLOT '/FMUT/clusters/CNV_audi_' effect_names{i} '_eta_cluster.fig'], ...
        'CNV Auditory');
end
close all

%% Extract post-hoc comparisons for significant CNV auditory clusters
fprintf('\n=== CNV AUDITORY POST-HOC COMPARISONS ===\n');
CNV_audi_posthoc = table();
for i = 1:length(effect_names)
    if strcmp(effect_names{i}, 'mixed')
        posthoc = extract_cluster_posthoc(GND_cue_SCNV_audi, GND_cue_SCNV_audi.F_tests, effect_names{i}, {'mixed'}, [3]);
    elseif strcmp(effect_names{i}, 'move')
        posthoc = extract_cluster_posthoc(GND_cue_SCNV_audi, GND_cue_SCNV_audi.F_tests, effect_names{i}, {'move'}, [3]);
    else
        posthoc = extract_cluster_posthoc(GND_cue_SCNV_audi, GND_cue_SCNV_audi.F_tests, effect_names{i}, {'mixed', 'move'}, [3, 3]);
    end

    % Add cluster mean plots
    save_path = [PATH_PLOT '/FMUT/cluster_means/CNV_audi_' effect_names{i} '_eta_cluster.fig'];
    create_cluster_mean_plots(GND_cue_SCNV_audi, GND_cue_SCNV_audi.F_tests, ...
        effect_names{i}, save_path, 'CNV auditory', 0.05);

    CNV_audi_posthoc = [CNV_audi_posthoc; posthoc];
end
if ~isempty(CNV_audi_posthoc)
    disp(CNV_audi_posthoc);
end
%%
GND_cue_SCNV_visu = FclustGND(GND_cue_visu, ...
                    'bins', 1:9, ...
                    'factor_names', {'mixed', 'move'}, ...
                    'factor_levels', [3, 3], ...
                    'time_wind', CNVStimWindow, ...
                    'include_chans',CNVChans, ...
                    'chan_hood', chanHood, ...
                    'n_perm', nPerm, ...
                    'plot_raster','yes', ...
                    'save_GND', 'no');

save([PATH_STAT '/FMUT/GND_GS_StimCNV_fullFactorial_cue_visu.mat'],'GND_cue_SCNV_visu');

%% Extract statistics for CNV visual
CNV_visu_stats = extract_statistics(GND_cue_SCNV_visu.F_tests, effect_names);
disp('CNV Visual Statistics:');
disp(CNV_visu_stats);

% Create cluster eta squared plots for CNV visual
for i = 1:length(effect_names)
    create_cluster_eta_plot(GND_cue_SCNV_visu.F_tests, effect_names{i}, ...
        [PATH_PLOT '/FMUT/clusters/CNV_visu_' effect_names{i} '_eta_cluster.fig'], ...
        'CNV Visual');
end
close all

% Extract post-hoc comparisons for significant CNV visual clusters
fprintf('\n=== CNV VISUAL POST-HOC COMPARISONS ===\n');
CNV_visu_posthoc = table();
for i = 1:length(effect_names)
    if strcmp(effect_names{i}, 'mixed')
        posthoc = extract_cluster_posthoc(GND_cue_SCNV_visu, GND_cue_SCNV_visu.F_tests, effect_names{i}, {'mixed'}, [3]);
    elseif strcmp(effect_names{i}, 'move')
        posthoc = extract_cluster_posthoc(GND_cue_SCNV_visu, GND_cue_SCNV_visu.F_tests, effect_names{i}, {'move'}, [3]);
    else
        posthoc = extract_cluster_posthoc(GND_cue_SCNV_visu, GND_cue_SCNV_visu.F_tests, effect_names{i}, {'mixed', 'move'}, [3, 3]);
    end

    % Add cluster mean plots
    save_path = [PATH_PLOT '/FMUT/cluster_means/CNV_visu_' effect_names{i} '_eta_cluster.fig'];
% Add this around your function calls to isolate the issue:
try
    fprintf('Starting cluster mean plots...\n');
    create_cluster_mean_plots(GND_cue_SCNV_visu, GND_cue_SCNV_visu.F_tests, ...
        effect_names{i}, save_path, 'CNV Visual', 0.05);
    fprintf('Cluster mean plots completed.\n');
catch ME
    fprintf('Error in cluster mean plots: %s\n', ME.message);
    rethrow(ME);
end
    CNV_visu_posthoc = [CNV_visu_posthoc; posthoc];
end
if ~isempty(CNV_visu_posthoc)
    disp(CNV_visu_posthoc);
end
%% Stim-related P3
GND_tar_SP3_audi = FclustGND(GND_tar_audi, ...
                'bins', 1:9, ...
                'factor_names', {'mixed', 'move'}, ...
                'factor_levels', [3, 3], ...
                'time_wind', P3StimWindow, ...
                'include_chans',P3Chans, ...
                'chan_hood', chanHood, ...
                'n_perm', nPerm, ...
                'plot_raster','yes', ...
                'save_GND', 'no');
save([PATH_STAT '/FMUT/GND_GS_StimP3_fullFactorial_tar_audi.mat'],'GND_tar_SP3_audi');

%% Extract statistics for P3 auditory
P3_audi_stats = extract_statistics(GND_tar_SP3_audi.F_tests, effect_names);
disp('P3 Auditory Statistics:');
disp(P3_audi_stats);

% Create cluster eta squared plots for P3 auditory
for i = 1:length(effect_names)
    create_cluster_eta_plot(GND_tar_SP3_audi.F_tests, effect_names{i}, ...
        [PATH_PLOT '/FMUT/clusters/P3_audi_' effect_names{i} '_eta_cluster.fig'], ...
        'P3 Auditory');
end
close all

% Extract post-hoc comparisons for significant P3 auditory clusters
fprintf('\n=== P3 AUDITORY POST-HOC COMPARISONS ===\n');
P3_audi_posthoc = table();
for i = 1:length(effect_names)
    if strcmp(effect_names{i}, 'mixed')
        posthoc = extract_cluster_posthoc(GND_tar_SP3_audi, GND_tar_SP3_audi.F_tests, effect_names{i}, {'mixed'}, [3]);
    elseif strcmp(effect_names{i}, 'move')
        posthoc = extract_cluster_posthoc(GND_tar_SP3_audi, GND_tar_SP3_audi.F_tests, effect_names{i}, {'move'}, [3]);
    else
        posthoc = extract_cluster_posthoc(GND_tar_SP3_audi, GND_tar_SP3_audi.F_tests, effect_names{i}, {'mixed', 'move'}, [3, 3]);
    end

    % Add cluster mean plots
    save_path = [PATH_PLOT '/FMUT/cluster_means/P3_audi_' effect_names{i} '_eta_cluster.fig'];
    create_cluster_mean_plots(GND_tar_SP3_audi, GND_tar_SP3_audi.F_tests, ...
        effect_names{i}, save_path, 'P3 Auditory', 0.05);

    P3_audi_posthoc = [P3_audi_posthoc; posthoc];
end
if ~isempty(P3_audi_posthoc)
    disp(P3_audi_posthoc);
end

%%
GND_tar_SP3_visu = FclustGND(GND_tar_visu, ...
                'bins', 1:9, ...
                'factor_names', {'mixed', 'move'}, ...
                'factor_levels', [3, 3], ...
                'time_wind', P3StimWindow, ...
                'include_chans',P3Chans, ...
                'chan_hood', chanHood, ...
                'n_perm', nPerm, ...
                'plot_raster','yes', ...
                'save_GND', 'no');
save([PATH_STAT '/FMUT/GND_GS_StimP3_fullFactorial_tar_visu.mat'],'GND_tar_SP3_visu');

%% Extract statistics for P3 visual
P3_visu_stats = extract_statistics(GND_tar_SP3_visu.F_tests, effect_names);
disp('P3 Visual Statistics:');
disp(P3_visu_stats);

% Create cluster eta squared plots for P3 visual
for i = 1:length(effect_names)
    create_cluster_eta_plot(GND_tar_SP3_visu.F_tests, effect_names{i}, ...
        [PATH_PLOT '/FMUT/clusters/P3_visu_' effect_names{i} '_eta_cluster.fig'], ...
        'P3 Visual');
end
close all;

% Extract post-hoc comparisons for significant P3 visual clusters
fprintf('\n=== P3 VISUAL POST-HOC COMPARISONS ===\n');
P3_visu_posthoc = table();
for i = 1:length(effect_names)
    if strcmp(effect_names{i}, 'mixed')
        posthoc = extract_cluster_posthoc(GND_tar_SP3_visu, GND_tar_SP3_visu.F_tests, effect_names{i}, {'mixed'}, [3]);
    elseif strcmp(effect_names{i}, 'move')
        posthoc = extract_cluster_posthoc(GND_tar_SP3_visu, GND_tar_SP3_visu.F_tests, effect_names{i}, {'move'}, [3]);
    else
        posthoc = extract_cluster_posthoc(GND_tar_SP3_visu, GND_tar_SP3_visu.F_tests, effect_names{i}, {'mixed', 'move'}, [3, 3]);
    end

    % Add cluster mean plots
    save_path = [PATH_PLOT '/FMUT/cluster_means/P3_visu_' effect_names{i} '_eta_cluster.fig'];
    create_cluster_mean_plots(GND_tar_SP3_visu, GND_tar_SP3_visu.F_tests, ...
        effect_names{i}, save_path, 'P3 Visual', 0.05);

    P3_visu_posthoc = [P3_visu_posthoc; posthoc];
end
if ~isempty(P3_visu_posthoc)
    disp(P3_visu_posthoc);
end

%% Save all statistics to file
all_stats = struct();
all_stats.CNV_audi = CNV_audi_stats;
all_stats.CNV_visu = CNV_visu_stats;
all_stats.P3_audi = P3_audi_stats;
all_stats.P3_visu = P3_visu_stats;

% Save post-hoc comparisons
%% Add Component column and concatenate
CNV_audi_posthoc.Component = repmat("CNV_audi", height(CNV_audi_posthoc), 1);
CNV_visu_posthoc.Component = repmat("CNV_visu", height(CNV_visu_posthoc), 1);
P3_audi_posthoc.Component  = repmat("P3_audi",  height(P3_audi_posthoc),  1);
P3_visu_posthoc.Component  = repmat("P3_visu",  height(P3_visu_posthoc),  1);

all_posthoc_table = [
    CNV_audi_posthoc;
    CNV_visu_posthoc;
    P3_audi_posthoc;
    P3_visu_posthoc
];

%% Write out CSV (now with Component as first column)
if ~isempty(all_posthoc_table)
    % Move 'Component' to be the first column in the table
    all_posthoc_table = movevars(all_posthoc_table, 'Component', 'Before', 1);
    
    writetable(all_posthoc_table, fullfile(PATH_STAT,'FMUT','posthoc_comparisons.csv'));
    fprintf('Post-hoc comparisons saved to: %s\n', fullfile(PATH_STAT,'FMUT','posthoc_comparisons.csv'));
end
%% Create a comprehensive summary table
all_effects_table = [CNV_audi_stats; CNV_visu_stats; P3_audi_stats; P3_visu_stats];
% Add component and modality columns for easier reading
n_rows = height(all_effects_table);
component = [repmat({'CNV'}, height(CNV_audi_stats), 1); 
             repmat({'CNV'}, height(CNV_visu_stats), 1);
             repmat({'P3'}, height(P3_audi_stats), 1);
             repmat({'P3'}, height(P3_visu_stats), 1)];
modality = [repmat({'Auditory'}, height(CNV_audi_stats), 1); 
            repmat({'Visual'}, height(CNV_visu_stats), 1);
            repmat({'Auditory'}, height(P3_audi_stats), 1);
            repmat({'Visual'}, height(P3_visu_stats), 1)];

comprehensive_table = [table(component, modality, 'VariableNames', {'Component', 'Modality'}), all_effects_table];

fprintf('\n=== COMPREHENSIVE STATISTICS TABLE ===\n');
disp(comprehensive_table);

% Save comprehensive table
writetable(comprehensive_table, [PATH_STAT '/FMUT/comprehensive_statistics.csv']);
fprintf('\nComprehensive statistics saved to: %s\n', [PATH_STAT '/FMUT/comprehensive_statistics.csv']);

% Create a simplified table for publication (reviewer Table 5 format)
% For multiple clusters, report the most significant (lowest p-value) cluster per effect
simplified_stats = table();
publication_stats_all_clusters = table();

for i = 1:height(comprehensive_table)
    row = comprehensive_table(i,:);
    if strcmp(row.Status{1}, 'significant')
        % Add to comprehensive table (all clusters)
        pub_row_all = table(row.Component, row.Modality, row.Effect, ...
            row.df_effect, row.df_error, row.F_max, row.p_value, row.eta_sq_max, row.cluster_number, ...
            'VariableNames', {'Component', 'Modality', 'Effect', 'df_effect', 'df_error', 'F_value', 'p_value', 'eta_squared_partial', 'cluster_number'});
        publication_stats_all_clusters = [publication_stats_all_clusters; pub_row_all];
    end
end

% For simplified table, keep only the most significant cluster per effect
if ~isempty(publication_stats_all_clusters)
    % Get unique combinations of Component, Modality, and base Effect name (without cluster suffix)
    effect_base_names = publication_stats_all_clusters.Effect;
    
    % Remove cluster suffixes to group effects
    for i = 1:length(effect_base_names)
        effect_name = effect_base_names{i};
        if contains(effect_name, '_cluster')
            effect_base_names{i} = extractBefore(effect_name, '_cluster');
        end
    end
    
    publication_stats_all_clusters.EffectBase = effect_base_names;
    
    % Group by Component, Modality, and EffectBase, keep row with lowest p-value
    [unique_effects, ~, group_idx] = unique(strcat(publication_stats_all_clusters.Component, '_', ...
        publication_stats_all_clusters.Modality, '_', publication_stats_all_clusters.EffectBase), 'stable');
    
    for g = 1:length(unique_effects)
        group_rows = publication_stats_all_clusters(group_idx == g, :);
        [~, min_p_idx] = min(group_rows.p_value);
        best_row = group_rows(min_p_idx, :);
        
        % Create simplified row (remove cluster info for cleaner table)
        simple_row = table(best_row.Component, best_row.Modality, best_row.EffectBase, ...
            best_row.df_effect, best_row.df_error, best_row.F_value, best_row.p_value, best_row.eta_squared_partial, ...
            'VariableNames', {'Component', 'Modality', 'Effect', 'df_effect', 'df_error', 'F_value', 'p_value', 'eta_squared_partial'});
        simplified_stats = [simplified_stats; simple_row];
    end
end

fprintf('\n=== SIMPLIFIED TABLE FOR PUBLICATION (Table 5 format - Most significant cluster per effect) ===\n');
disp(simplified_stats);

fprintf('\n=== ALL CLUSTERS TABLE FOR COMPREHENSIVE REPORTING ===\n');
disp(publication_stats_all_clusters);

% Save both tables for publication
writetable(simplified_stats, [PATH_STAT '/FMUT/publication_table5_simplified.csv']);
writetable(publication_stats_all_clusters, [PATH_STAT '/FMUT/publication_table5_all_clusters.csv']);
fprintf('\nSimplified publication table (best cluster per effect) saved to: %s\n', [PATH_STAT '/FMUT/publication_table5_simplified.csv']);
fprintf('All clusters publication table saved to: %s\n', [PATH_STAT '/FMUT/publication_table5_all_clusters.csv']);

%% Create plot directories if they don't exist
if ~exist([PATH_PLOT '/FMUT/ERP'], 'dir')
    mkdir([PATH_PLOT '/FMUT/ERP']);
end
if ~exist([PATH_PLOT '/FMUT/interactions'], 'dir')
    mkdir([PATH_PLOT '/FMUT/interactions']);
end
if ~exist([PATH_PLOT '/FMUT/effects'], 'dir')
    mkdir([PATH_PLOT '/FMUT/effects']);
end

%% Create plot directories if they don't exist
if ~exist([PATH_PLOT '/FMUT/ERP'], 'dir')
    mkdir([PATH_PLOT '/FMUT/ERP']);
end
if ~exist([PATH_PLOT '/FMUT/interactions'], 'dir')
    mkdir([PATH_PLOT '/FMUT/interactions']);
end
if ~exist([PATH_PLOT '/FMUT/effects'], 'dir')
    mkdir([PATH_PLOT '/FMUT/effects']);
end

%% Create ERP Line Plots using lw_niceLINEplot
% Note: Make sure lw_niceLINEplot function is in your MATLAB path
% If you get input parsing errors, check the parameter validation functions in lw_niceLINEplot

%% prepare dataset for plotting
% cnv audi
cnv_audi_repeatOnly = GND_cue_SCNV_audi.grands(:,:,1:3:end);
cnv_audi_repeatMix = GND_cue_SCNV_audi.grands(:,:,2:3:end);
cnv_audi_switchMix = GND_cue_SCNV_audi.grands(:,:,3:3:end);
% cnv visu
cnv_visu_repeatOnly = GND_cue_SCNV_visu.grands(:,:,1:3:end);
cnv_visu_repeatMix = GND_cue_SCNV_visu.grands(:,:,2:3:end);
cnv_visu_switchMix = GND_cue_SCNV_visu.grands(:,:,3:3:end);

% p3 audi
p3_audi_repeatOnly = GND_tar_SP3_audi.grands(:,:,1:3:end);
p3_audi_repeatMix = GND_tar_SP3_audi.grands(:,:,2:3:end);
p3_audi_switchMix = GND_tar_SP3_audi.grands(:,:,3:3:end);
% p3 visu
p3_visu_repeatOnly = GND_tar_SP3_visu.grands(:,:,1:3:end);
p3_visu_repeatMix = GND_tar_SP3_visu.grands(:,:,2:3:end);
p3_visu_switchMix = GND_tar_SP3_visu.grands(:,:,3:3:end);

% S1 and S2 are each 64×1000×3 struct-arrays
n = size(cnv_audi_repeatOnly,3);            % = 3
tpl = cnv_audi_repeatOnly(1,1,1);           % pick one element as a template

% preallocate 64×1000×6
S = repmat(tpl, size(cnv_audi_repeatOnly,1), size(cnv_audi_repeatOnly,2), 2*n);

% fill interleaved: [A1,B1,A2,B2,A3,B3]
for k = 1:n
    % cnv - repeat-only
    cnv_repeatOnly(:,:,2*k-1) = cnv_audi_repeatOnly(:,:,k);   % A_k
    cnv_repeatOnly(:,:,2*k)   = cnv_visu_repeatOnly(:,:,k);   % B_k
    % cnv - repeat-mix
    cnv_repeatMix(:,:,2*k-1) = cnv_audi_repeatMix(:,:,k);   % A_k
    cnv_repeatMix(:,:,2*k)   = cnv_visu_repeatMix(:,:,k);   % B_k
    % cnv - switch-mix
    cnv_switchMix(:,:,2*k-1) = cnv_audi_switchMix(:,:,k);   % A_k
    cnv_switchMix(:,:,2*k)   = cnv_visu_switchMix(:,:,k);   % B_k

    % p3 - repeat-only
    p3_repeatOnly(:,:,2*k-1) = p3_audi_repeatOnly(:,:,k);   % A_k
    p3_repeatOnly(:,:,2*k)   = p3_visu_repeatOnly(:,:,k);   % B_k
    % p3 - repeat-mix
    p3_repeatMix(:,:,2*k-1) = p3_audi_repeatMix(:,:,k);   % A_k
    p3_repeatMix(:,:,2*k)   = p3_visu_repeatMix(:,:,k);   % B_k
    % p3 - switch-mix
    p3_switchMix(:,:,2*k-1) = p3_audi_switchMix(:,:,k);   % A_k
    p3_switchMix(:,:,2*k)   = p3_visu_switchMix(:,:,k);   % B_k

end

%% CNV Auditory ERP Plot - Mixed Effect
figure;
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 25, 15], 'PaperUnits', 'centimeters', 'PaperSize', [26,16])
subplot(311)
    lw_niceLINEplot(GND_cue_SCNV_audi.chanlocs, CNVChans, cnv_repeatOnly, 2, ...
                    CNVPlotWindow, GND_cue_SCNV_audi.time_pts, ...
                    'ERPs at frontal patch - repeat in repeat-only blocks', 'Time (ms)', 'Amplitude (µV)', ...
                    'label_legend', {'stand - auditory', 'stand - visual', 'walk - auditory', 'walk - visual', 'perturbation - auditory', 'perturbation - visual'}, ...
                    'locs_legend', 'northwest',...
                    'xtick_fontsize',10,'ytick_fontsize',10,...
                    'polarity',-1,...
                    'avg_window',CNVStimWindow, ...
                    'ylim',CNVylim);
subplot(312)
    lw_niceLINEplot(GND_cue_SCNV_audi.chanlocs, CNVChans, cnv_repeatMix, 2, ...
                    CNVPlotWindow, GND_cue_SCNV_audi.time_pts, ...
                    'ERPs at frontal patch - repeat in mixed blocks', 'Time (ms)', 'Amplitude (µV)',...
                    'xtick_fontsize',10,'ytick_fontsize',10,...
                    'polarity',-1,...
                    'avg_window',CNVStimWindow,...
                    'ylim',CNVylim);
subplot(313)
    lw_niceLINEplot(GND_cue_SCNV_audi.chanlocs, CNVChans, cnv_switchMix, 2, ...
                    CNVPlotWindow, GND_cue_SCNV_audi.time_pts, ...
                    'ERPs at frontal patch - switch in mixed blocks', 'Time (ms)', 'Amplitude (µV)',...
                    'xtick_fontsize',10,'ytick_fontsize',10,...
                    'polarity',-1,...
                    'avg_window',CNVStimWindow,...
                    'ylim',CNVylim);
                    
saveas(gcf,[PATH_PLOT '/FMUT/ERP/CNV_audi_ERP.png']); 
saveas(gcf,[PATH_PLOT '/FMUT/ERP/CNV_audi_ERP.fig']); close gcf

%% P3 Auditory ERP Plot - Mixed Effect
figure;
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 25, 15], 'PaperUnits', 'centimeters', 'PaperSize', [26,16])
subplot(311)
    lw_niceLINEplot(GND_tar_SP3_audi.chanlocs, P3Chans, p3_repeatOnly, 2, ...
                    P3PlotWindow, GND_tar_SP3_audi.time_pts, ...
                    'ERPs at posterior patch - repeat in repeat-only blocks', 'Time (ms)', 'Amplitude (µV)', ...
                    'label_legend', {'stand - auditory', 'stand - visual', 'walk - auditory', 'walk - visual', 'perturbation - auditory', 'perturbation - visual'}, ...
                    'locs_legend', 'northwest',...
                    'xtick_fontsize',10,'ytick_fontsize',10,...
                    'polarity',-1,...
                    'avg_window',P3StimWindow,...
                    'ylim',P3ylim);
subplot(312)
    lw_niceLINEplot(GND_tar_SP3_audi.chanlocs, P3Chans, p3_repeatMix, 2, ...
                    P3PlotWindow, GND_tar_SP3_audi.time_pts, ...
                    'ERPs at posterior patch - repeat in mixed blocks', 'Time (ms)', 'Amplitude (µV)',...
                    'xtick_fontsize',10,'ytick_fontsize',10,...
                    'polarity',-1,...
                    'avg_window',P3StimWindow,...
                    'ylim',P3ylim);
subplot(313)
    lw_niceLINEplot(GND_tar_SP3_audi.chanlocs, P3Chans, p3_switchMix, 2, ...
                    P3PlotWindow, GND_tar_SP3_audi.time_pts, ...
                    'ERPs at posterior patch - switch in mixed blocks', 'Time (ms)', 'Amplitude (µV)',...
                    'xtick_fontsize',10,'ytick_fontsize',10,...
                    'polarity',-1,...
                    'avg_window',P3StimWindow,...
                    'ylim',P3ylim);
                    
saveas(gcf,[PATH_PLOT '/FMUT/ERP/P3_audi_ERP.png']); 
saveas(gcf,[PATH_PLOT '/FMUT/ERP/P3_audi_ERP.fig']); close gcf

%%
%% Create interaction plots using component-specific parameters
% CNV Auditory - mixed x move
create_interaction_plot(GND_cue_SCNV_audi, CNVStimWindow, CNVChans, ...
                        'Auditory cue CNV: Motor × Cognitive Difficulty Interaction', ...
                        [PATH_PLOT '/FMUT/interactions/CNV_audi_mixedXmove_interaction.fig'], ...
                        [-2, 1]);
% CNV Visual - mixed x move
create_interaction_plot(GND_cue_SCNV_visu, CNVStimWindow, CNVChans, ...
                        'Visual cue CNV: Motor × Cognitive Difficulty Interaction', ...
                        [PATH_PLOT '/FMUT/interactions/CNV_visu_mixedXmove_interaction.fig'], ...
                        [-2, 1]);

% P3 Auditory - mixed x move
create_interaction_plot(GND_tar_SP3_audi, P3StimWindow, P3Chans, ...
                        'Auditory cue P3: Motor × Cognitive Difficulty Interaction', ...
                        [PATH_PLOT '/FMUT/interactions/P3_audi_mixedXmove_interaction.fig'], ...
                        [0, 5]);
% P3 Visual - mixed x move
create_interaction_plot(GND_tar_SP3_visu, P3StimWindow, P3Chans, ...
                        'Visual cue P3: Motor × Cognitive Difficulty Interaction', ...
                        [PATH_PLOT '/FMUT/interactions/P3_visu_mixedXmove_interaction.fig'], ...
                        [0, 5]);

%% Create Interaction Plots for significant interactions using component-specific parameters
function create_interaction_plot(GND_struct, time_window, channels, title_str, save_path, y_limits, spacing)
% Function to create interaction plots with SEM error bars
% using individual ERP data of size
% channels × timepoints × conditions × subjects
%
% Inputs:
% GND_struct.indiv_erps – 4-D array [chan × time × cond × subj]
% time_window – [t_start, t_end]
% channels – cell array of channel labels to include
% title_str – plot title
% save_path – filename (with extension) to save figure
% y_limits – [ymin, ymax] for the y-axis
% spacing – horizontal offset between lines (default: 0.05)

%% Default spacing if not provided
if nargin < 7 || isempty(spacing)
    spacing = 0.1;
end

%% find channel indices
chan_names = {GND_struct.chanlocs.labels};
chan_indices = [];
for ch = 1:numel(channels)
    idx = find(strcmp(chan_names, channels{ch}));
    if ~isempty(idx)
        chan_indices(end+1) = idx;
    end
end

%% find time indices
t = GND_struct.time_pts;
tmask = (t >= time_window(1)) & (t <= time_window(2));
time_indices = find(tmask);

%% extract and average data over channels & time, but keep cond × subj
% indiv_erps is [chan × time × cond × subj]
erp4d = GND_struct.indiv_erps;
% subset: [channels × times × cond × subj]
sub = erp4d(chan_indices, time_indices, :, :);
% first average over channels (dim 1), then over time (dim 2):
m_chan = mean(sub, 1); % 1 × time × cond × subj
m_time = mean(m_chan, 2); % 1 × 1 × cond × subj
% squeeze to get [cond × subj], then transpose to [subj × cond]
tmp = squeeze(m_time); % cond × subj
tmp = tmp'; % subj × cond

%% compute mean and SEM across subjects
nSub = size(tmp,1);
mean_data = mean(tmp, 1); % 1 × 9
sem_data = std(tmp, 0, 1) ./ sqrt(nSub); % 1 × 9

%% reshape into 3 (mixed) × 3 (move)
mean_mat = reshape(mean_data, [3, 3]);
sem_mat = reshape(sem_data, [3, 3]);

%% plotting
figure;
set(gcf, 'Units','centimeters','Position',[0 0 20 15],...
    'PaperUnits','centimeters','PaperSize',[21 16]);

hold on;
colors = {'#272730','#B8357E','#235DB8'};
move_labels = {'stand','walk','perturbation'};
mixed_labels = {'repeat-only','mixed - repeat','mixed - switch'};

% Base x-coordinates
xvals = 1:3;

% Create offset x-coordinates for each line
offsets = [-spacing, 0, spacing]; % Left, center, right

for i = 1:3
    xvals_offset = xvals + offsets(i);
    errorbar(xvals_offset, mean_mat(i,:), sem_mat(i,:), ...
        '-o', 'LineWidth', 2, 'MarkerSize', 8, ...
        'Color', colors{i}, 'MarkerFaceColor', colors{i});
end

hold off;

%% labels, title, legend, styling
xlabel('Movement Condition','FontSize',16);
ylabel('Amplitude (µV)','FontSize',16);
title(title_str,'FontSize',18);
set(gca, 'XTick', xvals, 'XTickLabel', move_labels, 'FontSize',14);
legend(mixed_labels,'Location','best','FontSize',14);
ylim(y_limits);
grid on;
set(gcf,'Color','w');
set(gca,'TickDir','out');

%% save and close
saveas(gcf, save_path);
close(gcf);
end

%% Create summary statistics table
fprintf('\n=== COMPREHENSIVE STATISTICS ===\n');
fprintf('Statistics extracted from significant clusters (p < 0.05)\n');
fprintf('n.s. = not significant (no significant cluster found)\n');
fprintf('F values: max/min/mean within cluster\n');
fprintf('η²p: partial eta squared (max/min/mean/cluster_mean ± SD)\n\n');

fprintf('CNV Auditory:\n');
disp(CNV_audi_stats);
fprintf('\nCNV Visual:\n');
disp(CNV_visu_stats);
fprintf('\nP3 Auditory:\n');
disp(P3_audi_stats);
fprintf('\nP3 Visual:\n');
disp(P3_visu_stats);

% Save all statistics to file
all_stats = struct();
all_stats.CNV_audi = CNV_audi_stats;
all_stats.CNV_visu = CNV_visu_stats;
all_stats.P3_audi = P3_audi_stats;
all_stats.P3_visu = P3_visu_stats;

save([PATH_STAT '/FMUT/all_effect_statistics.mat'], 'all_stats');

% Create a comprehensive summary table
all_effects_table = [CNV_audi_stats; CNV_visu_stats; P3_audi_stats; P3_visu_stats];
% Add component and modality columns for easier reading
n_rows = height(all_effects_table);
component = [repmat({'CNV'}, height(CNV_audi_stats), 1); 
             repmat({'CNV'}, height(CNV_visu_stats), 1);
             repmat({'P3'}, height(P3_audi_stats), 1);
             repmat({'P3'}, height(P3_visu_stats), 1)];
modality = [repmat({'Auditory'}, height(CNV_audi_stats), 1); 
            repmat({'Visual'}, height(CNV_visu_stats), 1);
            repmat({'Auditory'}, height(P3_audi_stats), 1);
            repmat({'Visual'}, height(P3_visu_stats), 1)];

comprehensive_table = [table(component, modality, 'VariableNames', {'Component', 'Modality'}), all_effects_table];

fprintf('\n=== COMPREHENSIVE STATISTICS TABLE ===\n');
disp(comprehensive_table);

% Save comprehensive table
writetable(comprehensive_table, [PATH_STAT '/FMUT/comprehensive_statistics.csv']);
fprintf('\nComprehensive statistics saved to: %s\n', [PATH_STAT '/FMUT/comprehensive_statistics.csv']);

% Create a simplified table for publication (reviewer Table 5 format)
% For multiple clusters, report the most significant (lowest p-value) cluster per effect
simplified_stats = table();
publication_stats_all_clusters = table();

for i = 1:height(comprehensive_table)
    row = comprehensive_table(i,:);
    if strcmp(row.Status{1}, 'significant')
        % Add to comprehensive table (all clusters)
        pub_row_all = table(row.Component, row.Modality, row.Effect, ...
            row.df_effect, row.df_error, row.F_max, row.p_value, row.eta_sq_max, row.cluster_number, ...
            'VariableNames', {'Component', 'Modality', 'Effect', 'df_effect', 'df_error', 'F_value', 'p_value', 'eta_squared_partial', 'cluster_number'});
        publication_stats_all_clusters = [publication_stats_all_clusters; pub_row_all];
    end
end

% For simplified table, keep only the most significant cluster per effect
if ~isempty(publication_stats_all_clusters)
    % Get unique combinations of Component, Modality, and base Effect name (without cluster suffix)
    effect_base_names = publication_stats_all_clusters.Effect;
    
    % Remove cluster suffixes to group effects
    for i = 1:length(effect_base_names)
        effect_name = effect_base_names{i};
        if contains(effect_name, '_cluster')
            effect_base_names{i} = extractBefore(effect_name, '_cluster');
        end
    end
    
    publication_stats_all_clusters.EffectBase = effect_base_names;
    
    % Group by Component, Modality, and EffectBase, keep row with lowest p-value
    [unique_effects, ~, group_idx] = unique(strcat(publication_stats_all_clusters.Component, '_', ...
        publication_stats_all_clusters.Modality, '_', publication_stats_all_clusters.EffectBase), 'stable');
    
    for g = 1:length(unique_effects)
        group_rows = publication_stats_all_clusters(group_idx == g, :);
        [~, min_p_idx] = min(group_rows.p_value);
        best_row = group_rows(min_p_idx, :);
        
        % Create simplified row (remove cluster info for cleaner table)
        simple_row = table(best_row.Component, best_row.Modality, best_row.EffectBase, ...
            best_row.df_effect, best_row.df_error, best_row.F_value, best_row.p_value, best_row.eta_squared_partial, ...
            'VariableNames', {'Component', 'Modality', 'Effect', 'df_effect', 'df_error', 'F_value', 'p_value', 'eta_squared_partial'});
        simplified_stats = [simplified_stats; simple_row];
    end
end

fprintf('\n=== SIMPLIFIED TABLE FOR PUBLICATION (Table 5 format - Most significant cluster per effect) ===\n');
disp(simplified_stats);

fprintf('\n=== ALL CLUSTERS TABLE FOR COMPREHENSIVE REPORTING ===\n');
disp(publication_stats_all_clusters);

% Save both tables for publication
writetable(simplified_stats, [PATH_STAT '/FMUT/publication_table5_simplified.csv']);
writetable(publication_stats_all_clusters, [PATH_STAT '/FMUT/publication_table5_all_clusters.csv']);
fprintf('\nSimplified publication table (best cluster per effect) saved to: %s\n', [PATH_STAT '/FMUT/publication_table5_simplified.csv']);
fprintf('All clusters publication table saved to: %s\n', [PATH_STAT '/FMUT/publication_table5_all_clusters.csv']);

%% Generate comprehensive descriptive statistics table
% Add this section to your existing script after all the F-tests are completed

fprintf('\n=== GENERATING COMPREHENSIVE DESCRIPTIVE STATISTICS ===\n');

% Call the comprehensive descriptive statistics function
descriptive_table = generate_comprehensive_descriptive_stats(...
    GND_cue_audi, GND_cue_visu, GND_tar_audi, GND_tar_visu, ...
    CNVStimWindow, CNVChans, P3StimWindow, P3Chans, PATH_STAT);

%% Create formatted publication-ready tables

% Create a more compact table for publication
fprintf('\n=== CREATING PUBLICATION-READY TABLES ===\n');

% Interaction table with proper formatting
interaction_data = descriptive_table(strcmp(descriptive_table.Effect_Type, 'Interaction'), :);

% Reshape interaction data into a more readable format
components = unique(interaction_data.Component);
modalities = unique(interaction_data.Modality);

% Pre-allocate arrays for consistent data types
all_components = {};
all_modalities = {};
all_conditions = {};
all_mean_sd = {};
all_mean_uv = [];
all_sd_uv = [];
all_n = [];

% Display interaction means in matrix format for each component-modality
for comp = 1:length(components)
    for modal = 1:length(modalities)
        fprintf('\n--- %s %s: Interaction Means (µV) ---\n', components{comp}, modalities{modal});
        fprintf('%-15s %-10s %-10s %-10s\n', 'Mixed\\Move', 'Stand', 'Walk', 'Perturb');
        fprintf('%-15s %-10s %-10s %-10s\n', '----------', '-----', '----', '-------');
        
        subset = interaction_data(strcmp(interaction_data.Component, components{comp}) & ...
                                 strcmp(interaction_data.Modality, modalities{modal}), :);
        
        if ~isempty(subset)
            % Collect data for this component-modality combination
            for row = 1:height(subset)
                all_components{end+1} = subset.Component{row};
                all_modalities{end+1} = subset.Modality{row};
                all_conditions{end+1} = subset.Condition{row};
                all_mean_sd{end+1} = sprintf('%.3f (%.3f)', subset.Mean_uV(row), subset.SD_uV(row));
                all_mean_uv(end+1) = subset.Mean_uV(row);
                all_sd_uv(end+1) = subset.SD_uV(row);
                all_n(end+1) = subset.N(row);
            end
        end
    end
end

% Create table all at once with consistent data types
publication_table = table(all_components', all_modalities', all_conditions', all_mean_sd', ...
                         all_mean_uv', all_sd_uv', all_n', ...
                         'VariableNames', {'Component', 'Modality', 'Condition', 'Mean_SD', 'Mean_uV', 'SD_uV', 'N'});

% Save publication table
writetable(publication_table, [PATH_STAT '/FMUT/publication_descriptive_table.csv']);
fprintf('Publication-ready descriptive table saved to: %s\n', [PATH_STAT '/FMUT/publication_descriptive_table.csv']);

%% Create main effects summary tables
% Process main effects
main_effects_data = descriptive_table(~strcmp(descriptive_table.Effect_Type, 'Interaction'), :);

% Pre-allocate arrays for consistent data types
main_comp_mod = {};
main_effect_types = {};
main_conditions = {};
main_mean_sd = {};
main_mean_uv = [];
main_sd_uv = [];
main_n = [];

for i = 1:height(main_effects_data)
    row = main_effects_data(i, :);
    main_comp_mod{end+1} = sprintf('%s_%s', row.Component{1}, row.Modality{1});
    main_effect_types{end+1} = row.Effect_Type{1};
    main_conditions{end+1} = row.Condition{1};
    main_mean_sd{end+1} = sprintf('%.3f (%.3f)', row.Mean_uV, row.SD_uV);
    main_mean_uv(end+1) = row.Mean_uV;
    main_sd_uv(end+1) = row.SD_uV;
    main_n(end+1) = row.N;
end

% Create table all at once
main_effects_summary = table(main_comp_mod', main_effect_types', main_conditions', main_mean_sd', ...
                            main_mean_uv', main_sd_uv', main_n', ...
                            'VariableNames', {'Component_Modality', 'Effect_Type', 'Condition', 'Mean_SD', 'Mean_uV', 'SD_uV', 'N'});

writetable(main_effects_summary, [PATH_STAT '/FMUT/main_effects_summary.csv']);
fprintf('Main effects summary table saved to: %s\n', [PATH_STAT '/FMUT/main_effects_summary.csv']);

%% Display key summary statistics
fprintf('\n=== KEY DESCRIPTIVE STATISTICS SUMMARY ===\n');

% Display interaction means in matrix format for each component-modality
for comp = 1:length(components)
    for mod = 1:length(modalities)
        fprintf('\n--- %s %s: Interaction Means (µV) ---\n', components{comp}, modalities{modal});
        fprintf('%-15s %-10s %-10s %-10s\n', 'Mixed\\Move', 'Stand', 'Walk', 'Perturb');
        fprintf('%-15s %-10s %-10s %-10s\n', '----------', '-----', '----', '-------');
        
        subset = interaction_data(strcmp(interaction_data.Component, components{comp}) & ...
                                 strcmp(interaction_data.Modality, modalities{modal}), :);
        
        if height(subset) >= 9
            % Reshape into 3x3 matrix (mixed x movement)
            means_matrix = reshape([subset.Mean_uV], 3, 3);
            mixed_labels = {'repeat-only', 'mixed-repeat', 'mixed-switch'};
            
            for mixed = 1:3
                fprintf('%-15s %-10.3f %-10.3f %-10.3f\n', mixed_labels{mixed}, ...
                        means_matrix(mixed, 1), means_matrix(mixed, 2), means_matrix(mixed, 3));
            end
        end
    end
end

fprintf('\nAll descriptive statistics have been generated and saved to the FMUT folder.\n');

%% Function to generate comprehensive descriptive statistics table
function [descriptive_table] = generate_comprehensive_descriptive_stats(GND_cue_audi, GND_cue_visu, GND_tar_audi, GND_tar_visu, CNVStimWindow, CNVChans, P3StimWindow, P3Chans, PATH_STAT)
% Generate comprehensive descriptive statistics table for all conditions
% 
% Inputs:
%   GND_* structures from the main analysis
%   CNVStimWindow, P3StimWindow: time windows for analysis
%   CNVChans, P3Chans: channel groups for analysis
%   PATH_STAT: path to save output files
%
% Output:
%   descriptive_table: comprehensive table with means/SDs for all conditions

fprintf('Generating comprehensive descriptive statistics...\n');

%% Define experimental conditions
mixed_labels = {'repeat-only', 'mixed-repeat', 'mixed-switch'};
move_labels = {'stand', 'walk', 'perturbation'};
interaction_labels = {'Single-Stand', 'Single-Walk', 'Single-Pert', ...
                     'Repeat-Stand', 'Repeat-Walk', 'Repeat-Pert', ...
                     'Switch-Stand', 'Switch-Walk', 'Switch-Pert'};

%% Initialize results storage
all_results = [];

%% Analysis configurations
analyses = {
    struct('GND', GND_cue_audi, 'component', 'CNV', 'modality', 'Auditory', 'time_window', CNVStimWindow, 'channels', {CNVChans}), ...
    struct('GND', GND_cue_visu, 'component', 'CNV', 'modality', 'Visual', 'time_window', CNVStimWindow, 'channels', {CNVChans}), ...
    struct('GND', GND_tar_audi, 'component', 'P3', 'modality', 'Auditory', 'time_window', P3StimWindow, 'channels', {P3Chans}), ...
    struct('GND', GND_tar_visu, 'component', 'P3', 'modality', 'Visual', 'time_window', P3StimWindow, 'channels', {P3Chans})
};

%% Process each analysis
for a = 1:length(analyses)
    analysis = analyses{a};
    fprintf('Processing %s %s...\n', analysis.component, analysis.modality);
    
    % Extract amplitude data for this analysis
    [condition_data, mixed_main_data, move_main_data] = extract_amplitude_data(analysis.GND, analysis.time_window, analysis.channels);
    
    %% 1. Individual condition means and SDs (9 conditions)
    for cond = 1:9
        if ~isempty(condition_data{cond})
            result_row = struct();
            result_row.Component = analysis.component;
            result_row.Modality = analysis.modality;
            result_row.Effect_Type = 'Interaction';
            result_row.Condition = interaction_labels{cond};
            result_row.Mean_uV = mean(condition_data{cond});
            result_row.SD_uV = std(condition_data{cond});
            result_row.N = length(condition_data{cond});
            result_row.SEM_uV = result_row.SD_uV / sqrt(result_row.N);
            
            % Add factor level information
            mixed_level = ceil(cond/3); % 1, 2, or 3
            move_level = mod(cond-1, 3) + 1; % 1, 2, or 3
            result_row.Mixed_Level = mixed_labels{mixed_level};
            result_row.Move_Level = move_labels{move_level};
            
            all_results = [all_results; result_row];
        end
    end
    
    %% 2. Mixed task main effect (3 levels)
    for mixed = 1:3
        if ~isempty(mixed_main_data{mixed})
            result_row = struct();
            result_row.Component = analysis.component;
            result_row.Modality = analysis.modality;
            result_row.Effect_Type = 'Mixed_Main_Effect';
            result_row.Condition = mixed_labels{mixed};
            result_row.Mean_uV = mean(mixed_main_data{mixed});
            result_row.SD_uV = std(mixed_main_data{mixed});
            result_row.N = length(mixed_main_data{mixed});
            result_row.SEM_uV = result_row.SD_uV / sqrt(result_row.N);
            result_row.Mixed_Level = mixed_labels{mixed};
            result_row.Move_Level = 'All'; % Collapsed across movement
            
            all_results = [all_results; result_row];
        end
    end
    
    %% 3. Movement main effect (3 levels)
    for move = 1:3
        if ~isempty(move_main_data{move})
            result_row = struct();
            result_row.Component = analysis.component;
            result_row.Modality = analysis.modality;
            result_row.Effect_Type = 'Move_Main_Effect';
            result_row.Condition = move_labels{move};
            result_row.Mean_uV = mean(move_main_data{move});
            result_row.SD_uV = std(move_main_data{move});
            result_row.N = length(move_main_data{move});
            result_row.SEM_uV = result_row.SD_uV / sqrt(result_row.N);
            result_row.Mixed_Level = 'All'; % Collapsed across mixed task
            result_row.Move_Level = move_labels{move};
            
            all_results = [all_results; result_row];
        end
    end
end

%% Convert to table
descriptive_table = struct2table(all_results);

%% Reorder columns for better readability
descriptive_table = descriptive_table(:, {'Component', 'Modality', 'Effect_Type', 'Condition', ...
                                         'Mixed_Level', 'Move_Level', 'Mean_uV', 'SD_uV', 'SEM_uV', 'N'});

%% Display results
fprintf('\n=== COMPREHENSIVE DESCRIPTIVE STATISTICS ===\n');
disp(descriptive_table);

%% Save to files
% Save complete table
writetable(descriptive_table, [PATH_STAT '/FMUT/comprehensive_descriptive_statistics.csv']);
fprintf('\nComprehensive descriptive statistics saved to: %s\n', [PATH_STAT '/FMUT/comprehensive_descriptive_statistics.csv']);

% Create separate tables for different effect types
interaction_table = descriptive_table(strcmp(descriptive_table.Effect_Type, 'Interaction'), :);
mixed_main_table = descriptive_table(strcmp(descriptive_table.Effect_Type, 'Mixed_Main_Effect'), :);
move_main_table = descriptive_table(strcmp(descriptive_table.Effect_Type, 'Move_Main_Effect'), :);

writetable(interaction_table, [PATH_STAT '/FMUT/interaction_descriptive_stats.csv']);
writetable(mixed_main_table, [PATH_STAT '/FMUT/mixed_main_effect_descriptive_stats.csv']);
writetable(move_main_table, [PATH_STAT '/FMUT/move_main_effect_descriptive_stats.csv']);

fprintf('Separate effect tables saved to FMUT folder\n');

%% Create summary tables by component and modality
fprintf('\n=== SUMMARY BY COMPONENT AND MODALITY ===\n');

% Create pivot-style summaries
components = unique(descriptive_table.Component);
modalities = unique(descriptive_table.Modality);

for comp = 1:length(components)
    for modal = 1:length(modalities)
        subset = descriptive_table(strcmp(descriptive_table.Component, components{comp}) & ...
                                  strcmp(descriptive_table.Modality, modalities{modal}), :);
        
        if ~isempty(subset)
            fprintf('\n--- %s %s ---\n', components{comp}, modalities{modal});
            
            % Interaction effects
            int_subset = subset(strcmp(subset.Effect_Type, 'Interaction'), :);
            if ~isempty(int_subset)
                fprintf('Interaction Effects (all 9 conditions):\n');
                fprintf('%-20s %-12s %-12s %-8s\n', 'Condition', 'Mean(µV)', 'SD(µV)', 'N');
                fprintf('%-20s %-12s %-12s %-8s\n', '--------', '--------', '------', '-');
                for i = 1:height(int_subset)
                    fprintf('%-20s %-12.3f %-12.3f %-8d\n', int_subset.Condition{i}, ...
                            int_subset.Mean_uV(i), int_subset.SD_uV(i), int_subset.N(i));
                end
            end
            
            % Main effects
            main_subset = subset(~strcmp(subset.Effect_Type, 'Interaction'), :);
            if ~isempty(main_subset)
                fprintf('\nMain Effects:\n');
                fprintf('%-25s %-20s %-12s %-12s %-8s\n', 'Effect Type', 'Condition', 'Mean(µV)', 'SD(µV)', 'N');
                fprintf('%-25s %-20s %-12s %-12s %-8s\n', '----------', '--------', '--------', '------', '-');
                for i = 1:height(main_subset)
                    fprintf('%-25s %-20s %-12.3f %-12.3f %-8d\n', main_subset.Effect_Type{i}, ...
                            main_subset.Condition{i}, main_subset.Mean_uV(i), main_subset.SD_uV(i), main_subset.N(i));
                end
            end
        end
    end
end

end

%% Helper function to extract amplitude data
function [condition_data, mixed_main_data, move_main_data] = extract_amplitude_data(GND_struct, time_window, channels)
% Extract amplitude data for specified time window and channels
%
% Outputs:
%   condition_data: cell array with data for each of 9 conditions
%   mixed_main_data: cell array with data for 3 mixed task levels (collapsed across movement)
%   move_main_data: cell array with data for 3 movement levels (collapsed across mixed task)

% Find channel indices
chan_names = {GND_struct.chanlocs.labels};
chan_indices = [];
for ch = 1:length(channels)
    idx = find(strcmp(chan_names, channels{ch}));
    if ~isempty(idx)
        chan_indices(end+1) = idx;
    end
end

% Find time indices
time_mask = GND_struct.time_pts >= time_window(1) & GND_struct.time_pts <= time_window(2);
time_indices = find(time_mask);

% Initialize output
condition_data = cell(1, 9);
mixed_main_data = cell(1, 3);
move_main_data = cell(1, 3);

% Extract data for each individual condition (9 conditions)
for cond = 1:9
    if cond <= size(GND_struct.indiv_erps, 3)
        cond_values = [];
        for subj = 1:size(GND_struct.indiv_erps, 4)
            % Extract data: [channels × timepoints × condition × subject]
            subj_data = squeeze(GND_struct.indiv_erps(chan_indices, time_indices, cond, subj));
            % Average across channels and time
            cond_values(end+1) = mean(subj_data(:));
        end
        condition_data{cond} = cond_values;
    end
end

% Mixed task main effect (collapse across movement conditions)
% Bins 1,4,7 = repeat-only; 2,5,8 = mixed-repeat; 3,6,9 = mixed-switch
mixed_bins = {[1,4,7], [2,5,8], [3,6,9]};
for mixed = 1:3
    mixed_values = [];
    for bin = mixed_bins{mixed}
        if bin <= size(GND_struct.indiv_erps, 3)
            for subj = 1:size(GND_struct.indiv_erps, 4)
                subj_data = squeeze(GND_struct.indiv_erps(chan_indices, time_indices, bin, subj));
                mixed_values(end+1) = mean(subj_data(:));
            end
        end
    end
    mixed_main_data{mixed} = mixed_values;
end

% Movement main effect (collapse across mixed task conditions)
% Bins 1,2,3 = stand; 4,5,6 = walk; 7,8,9 = perturbation
move_bins = {[1,2,3], [4,5,6], [7,8,9]};
for move = 1:3
    move_values = [];
    for bin = move_bins{move}
        if bin <= size(GND_struct.indiv_erps, 3)
            for subj = 1:size(GND_struct.indiv_erps, 4)
                subj_data = squeeze(GND_struct.indiv_erps(chan_indices, time_indices, bin, subj));
                move_values(end+1) = mean(subj_data(:));
            end
        end
    end
    move_main_data{move} = move_values;
end

end
