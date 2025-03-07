function [stat, p, p_ttest, avg, stat_all, p_all, mean_y, err_y, groups,x,y] = plot_in_situ_quantification(data,...
    field, x_order_label, c_order, norm_true, control_cond, triangle_tf, y_limits, control_norm)
%PLOT_IN_SITU_QUANTIFICATION Plot in situ quantification
% 
%   Inputs
%       data: data structure from quantify_in_situ
%       field: names of the field for averaging and plotting
%
%   Outputs
%       None
% 
%   Overview
%       This function plots the data from quantify_in_situ
    
    if ~iscell(data)
        data = {data};
    end

    if ~iscell(x_order_label)
        x_order_label = {x_order_label};
    end

    if ~iscell(c_order)
        c_order = {c_order};
    end

    if ~iscell(control_cond)
        control_cond = {control_cond};
    end

    avg = cell(size(data,2), 1);
    stat = cell(size(data,2), 1);
    p = cell(size(data,2), 1);
    p_ttest = cell(size(data,2), 1);

    for i = 1:size(data,2)
        % Find the unique conditions in data
        conditions = unique({data{i}.condition}', 'stable');
        
        % Calulate the mean
        avg{i} = calc_means(data{i}, conditions);
        
        if size(avg{i},2) > 1
            [stat{i}, p{i}, p_ttest{i}] = statistical_analysis(avg{i}, field, false);
        else
            stat{i} = [];
            p{i} = [];
            p_ttest{i} = [];
        end
        
        if ~norm_true
            % Plot the mean for the given field
            [stat_all, p_all, mean_y, err_y, groups] =  plot_mean_intensity(avg{i}, field, x_order_label{i}, c_order{i},...
                norm_true, control_cond{i}, triangle_tf, y_limits, control_norm);
        end
    end

    if norm_true
        % Plot the mean for the given field
        [stat_all, p_all, mean_y, err_y, groups, x, y] = plot_mean_intensity(avg, field, x_order_label{:}, c_order{:},...
            norm_true, control_cond, triangle_tf, y_limits, control_norm);
    end

%     stat_all = [];
%     p_all = [];
%     mean_y = [];
%     err_y = [];
%     groups = [];
end

function avg = calc_means(data, condition)
%CALC_MEANS Calculate the mean intensity
% 
%   Inputs
%       data: data structure from main function from quantify_in_situ
%       condition: the condition of treatment
% 
%   Outputs
%       avg: a structure containg the average data, the individual data,
%           and the error bar data
% 
%   Overview
%       Calculates the mean and error for data
            
    % Make structure for storing averaged data
    avg = struct('condition', [],...
                 'norm_I_avg', cell(3,size(condition, 1)),...
                 'sum_I_avg', [],...
                 'norm_sum_I_avg', [],...
                 'norm_A_avg', [],...
                 'norm_Maj_avg', [],...
                 'norm_Min_avg', [],...
                 'Maj_em_avg', [],...
                 'Min_em_avg', [],...
                 'A_em_avg', [],...
                 'I_em_avg', [],...
                 'Maj_im_avg', [],...
                 'Min_im_avg', [],...
                 'A_im_avg', [],...
                 'I_im_avg', [],...
                 'A_bg_avg', [],...
                 'I_bg_avg', [],...
                 'varI_im_avg', [],...
                 'disp_ind_im_avg', [],...
                 'norm_MaxFeret_avg', [],...
                 'norm_MinFeret_avg', [],...
                 'MaxFeret_im_avg', [],...
                 'MinFeret_im_avg', [],...
                 'MaxFeret_em_avg', [],...
                 'MinFeret_em_avg', []);
             
    % Save field names
    data_fields = fieldnames(avg);
    data_fields = data_fields(2:end);
               
    % For each entered condition
    for i = 1:size(condition, 1)
        % Find indices of data that match entered condition (logical
        % indexing)
        j = strcmp({data.condition}, condition{i,1});
        
        % For each field
        for k = 1:size(data_fields, 1)
            % Save the individual data
            avg(1,i).(data_fields{k}) = cat(1, data(j).(data_fields{k}));
            
            % Calculate the mean
            avg(2,i).(data_fields{k}) = mean(cat(1,...
                                        data(j).(data_fields{k})));
            
            % Calculate error for error bars
            avg(3,i).(data_fields{k}) = calc_error(cat(1, ...
                                        data(j).(data_fields{k})), 'SD');%SD%SEM
        end
        
        % Save the condition
        avg(1,i).condition = repmat(condition(i,1), size(avg(1,i).norm_I_avg));
        avg(2,i).condition = condition{i,1};
        avg(3,i).condition = condition{i,1};
    end
end

function err_d = calc_error(d, select_error)
%CALC_ERROR Calculates the error for data d
% 
%   Input
%       d: data points 
%       select_error: a string, either SEM, CI, or SD
% 
%   Output
%       err_d: the error for the data
% 
%   Overview
%       This function calculates the error for determining error bars. It
%       takes data d and the choice for calculating the error, either
%       standard error of the mean (SEM), 95% confidence intervals (CI), or
%       standard deviation (SD)
    
    % standard deviation for data in d
%     STD_d = nanstd(d, [], 1);
    STD_d = std(d, 0, 1, 'omitnan');

    % standard error of the mean for data in d
    SEM_d = STD_d ./ sqrt(sum(~isnan(d), 1));

    % confidence interval for data in d
    ts_d = tinv(0.975, sum(~isnan(d), 1) - 1);
    CI_d = ts_d .* SEM_d;

    % If user inputed SD
    if ~isempty(select_error) && isequal(select_error, 'SD')
        % Error is standard deviation
        err_d = STD_d;
    % Else if user inputed CI
    elseif ~isempty(select_error) && isequal(select_error, 'CI')
        % Error is confidence intervals
        err_d = CI_d;
    % Else if user inputed SEM or anything else
    else
        % Error is standard error of the mean
        err_d = SEM_d;
    end
end

function [stat, p, mean_y, err_y, unordered_groups, x, y] = plot_mean_intensity(avg, field, x_label_order, c_order,...
    norm_true, control_cond, change_marker, y_limits, control_norm)
%PLOT_MEAN_INTENSITY Plot the individual intensity, mean, and error bars
% 
%   Inputs
%       avg: the avg structure from calc_means
%       field: the field to be plotted
% 
%   Outputs
%       None
% 
%   Overview
%       A plot is made of the individual data, the mean, and the error bars
%       for the data specified by field
            
    % Array of colors
    colors = [   0.5,    0.5,    0.5;   % 1) gray
              0.4660, 0.6740, 0.1880;   % 2) green
              0.3010, 0.7450, 0.9330;   % 3) cyan
                   0, 0.4470, 0.7410;   % 4) blue
                   0, 0.2235, 0.3705;   % 5) dark blue
              0.4940, 0.1840, 0.5560;   % 6) purple
              0.8500, 0.3250, 0.0980;   % 7) orange
              0.9290, 0.6940, 0.1250;   % 8) yellow/gold
              0.5882, 0.2941,      0;   % 9) brown
              0.9500, 0.0800, 0.0980;   % 10) red
              0.9525, 0.1170, 0.5520;   % 11) pink-magenta
              0.6350, 0.0780, 0.1840;   % 12) dark red
                   0,      0,      0];  % 13) black
              
               
    % Concatenate the data (individual, mean, error, and condition) for the
    % specified field
    if norm_true
        y = cell(size(avg, 1),1);
        mean_y = cell(size(avg, 1),1);
        err_y = cell(size(avg, 1),1);
        g = cell(size(avg, 1),1);
        control_y = cell(size(avg, 1),1);
        control_ind = cell(size(avg, 1),1);

        for k = 1:size(avg, 1)
            control_y{k} = avg{k}(1,strcmp({avg{k}(2,:).condition}, control_cond{k})).(field)...
                / avg{k}(2,strcmp({avg{k}(2,:).condition}, control_cond{k})).(field);

            y{k} = cat(1, avg{k}(1,:).(field)) / avg{k}(2,strcmp({avg{k}(2,:).condition}, control_cond{k})).(field);
%             mean_y{k} = cat(1, avg{k}(2,:).(field)) / avg{k}(2,strcmp({avg{k}(2,:).condition}, control_cond{k})).(field);
%             err_y{k} = cat(1, avg{k}(3,:).(field)) / avg{k}(2,strcmp({avg{k}(2,:).condition}, control_cond{k})).(field);
            g{k} = cat(1, avg{k}(1,:).condition);
            control_ind{k} = strcmp({avg{k}(2,:).condition}, control_cond{k});
        end

%         control_y = cat(1, control_y{:});
        y = cat(1, y{:});% + mean(control_y);
%         mean_y = cat(1, mean_y{:});% + mean(control_y);
%         err_y = cat(1, err_y{:});
        g = cat(1, g{:});

%         control_ind = cat(2, control_ind{:});

%         ind_loc = find(control_ind);
%         err_y(ind_loc(1)) = calc_error(control_y, 'SD');

%         for m = size(ind_loc,2):-1:2
%             mean_y(ind_loc(m)) = [];
%             err_y(ind_loc(m)) = [];
%         end
    else
        if iscell(avg)
            y = cell(size(avg, 1),1);
            g = cell(size(avg, 1),1);
    
            for k = 1:size(avg, 1)
                y{k} = cat(1, avg{k}(1,:).(field));
                g{k} = cat(1, avg{k}(1,:).condition);
            end

            y = cat(1, y{:});
            g = cat(1, g{:});
        else
            y = cat(1, avg(1,:).(field));
            mean_y = cat(1, avg(2,:).(field));
            err_y = cat(1, avg(3,:).(field));
            g = cat(1, avg(1,:).condition);
        end
    end
    
    % Find the unique conditions and initialize variables
    unordered_groups = unique(g, 'stable');
    groups = unordered_groups(x_label_order);
    x = zeros(size(g,1), 1);
    c = zeros(size(g,1), 3);
    mean_x = zeros(size(groups,1), 1);

    if norm_true || iscell(avg)
        mean_y = zeros(size(groups,1), 1);
        err_y = zeros(size(groups,1), 1);
    end

    % For each condition
    for i = 1:size(groups, 1)
        % Create logical array where condition matches
        ind = strcmp(g, groups{i});
        mean_ind = strcmp(unordered_groups, groups{i});

        if norm_true || iscell(avg)
            % Calculate the mean
            mean_y(mean_ind) = mean(y(strcmp(g, groups{i})));
            
            % Calculate error for error bars
            err_y(mean_ind) = calc_error(y(strcmp(g, groups{i})), 'SD');%SD%SEM
        end
        
        % Save i for only the matched conditions to plot all points with
        % the same condition at the same x
        x(ind) = i;
        mean_x(mean_ind) = i;
        
        % Assign a color to each point. Color is picked by using mod in
        % case a color needs to be repeated. The loop saves each column of
        % colors
        for j = 1:size(colors,2)
            c(ind,j) = colors(mod(c_order(i)-1, size(colors,1))+1,j);
%             c(ind,j) = colors(8,j);
        end
    end

%     if any(strcmp(groups, '40'))
%         data_40 = y(strcmp(g, '40'));
%         data_YW = y(strcmp(g, 'YW'));
%         
%         combined_data = cat(1, data_YW, data_40);
% 
%         mean_y(strcmp(unordered_groups, '40')) = mean(combined_data);
%         mean_y(strcmp(unordered_groups, 'YW')) = mean(combined_data);
% 
%         err_y(strcmp(unordered_groups, '40')) = calc_error(combined_data, 'SEM');
%         err_y(strcmp(unordered_groups, 'YW')) = calc_error(combined_data, 'SEM');
% 
%         groups(strcmp(groups, '40')) = [];
%     end
    [stat, p] = statistical_analysis_all(y, g);
    
    % Make a figure and keep the axis for plotting the raw data, the mean,
    % and the error bars on the same plot
    figure;
    hold on
    plot([0, size(groups,1)+1], [mean_y(strcmp(unordered_groups, control_norm)),...
        mean_y(strcmp(unordered_groups, control_norm))], '--', 'Color', 'k', 'LineWidth', 2);
    
    if change_marker
        for j = 1:size(groups, 1)
            if any(j == (1:2:size(groups, 1)))
                % Plot raw data with jitter to offset the points with some transparency
                scatter(x(x==j), y(x==j), 100, c(x==j,:), 'filled', 'jitter', 'on', 'jitteramount', 0.2,...
                        'MarkerFaceAlpha', .5,'MarkerEdgeAlpha', .5);
            else
                % Plot raw data with jitter to offset the points with some transparency
                scatter(x(x==j), y(x==j), 100, c(x==j,:), '^', 'filled', 'jitter', 'on', 'jitteramount', 0.2,...
                        'MarkerFaceAlpha', .5,'MarkerEdgeAlpha', .5);
            end
        end
    else
        % Plot raw data with jitter to offset the points with some transparency
        scatter(x, y, 100, c, 'filled', 'jitter', 'on', 'jitteramount', 0.2,...
                'MarkerFaceAlpha', .5,'MarkerEdgeAlpha', .5);
    end
    
    % Plot the mean with the error bars and set properties
    h = errorbar(mean_x, mean_y, err_y, '.', 'Color', 'k');
    set(h, 'linewidth', 2, 'markersize', 25);

    hold off

    set(gca,'fontname','arial');
%     set(gca,'YScale','log');
    
    if isempty(y_limits)
        y_limits = [0.95 .* min(y), 1.05 .* max(y)];
    end

    % Set properties of axis
    set(gca, 'xlim', [0.5, size(groups,1)+0.5],...
             'xtick', 1:size(groups,1), ...
             'xticklabels', groups,...
             'XTickLabelRotation', 45,...
             'ylim', y_limits,...%[0.02,4] [-0.1,3] [40, 200] [0, 10000] [-3500, 7000]
             'fontsize', 20);  
end

function [stat, p, p_ttest] = statistical_analysis(avg, field, indiv_ttest)
%STATISTICAL_ANALYSIS Perform ANOVA to compare the means between conditions
% 
%   Input
%       avg: the structure returned from calc_means
% 
%   Output
%       stat: structure containing outputs from anova1 and multcompare
%       p: table of p-values for pairwise comparisons
% 
%   Overview
%       This function performs statistical analysis on the data.
%       Specifically, it performs one way ANOVA using anova1 and multiple
%       comparisons using Tukey's HSD using multcompare. It returns the
%       outputs from anova1 and multcompare in the structure stat and
%       a table of p-values, p, for pairwise comparison between conditions.
    
    % To calculate individual t tests between groups of two samples,
    % organize samples so samples to be tested are concatanated together
    % (for example column 1 and 2 will be tested, 3 and 4, etc)
    if indiv_ttest
        p_ttest = cell(2, size(avg,2)/2);
        ind = 1;
        
        for i = 1:2:size(avg,2)
            p_ttest{1,ind} = sprintf('%s & %s', avg(2,i).condition, avg(2,i+1).condition);
            [~, p_ttest{2,ind}] = ttest2(avg(1,i).(field), avg(1,i+1).(field), 'Vartype', 'unequal');
            ind = ind + 1;
        end
    else
        p_ttest = [];
    end

    condition = cat(1, avg(1,:).condition);
    data = cat(1, avg(1,:).(field));
    
    % Initialize a structure for storing the results of the statistical
    % analysis
    stat = struct('p', [],...
                  'tbl', [],...
                  'stats', [],...
                  'p_indiv', [],...
                  'means', [],... 
                  'names', []);
    
    % Perform ANOVA on the intensity data grouped by condition
    [stat.p, stat.tbl, stat.stats, stat.p_indiv, stat.means,...
        stat.names] = stat_test(log(data), condition);
    
    % Initialize variables for making comparison tables
    p = cell(size(stat.names, 1), size(stat.names, 1));

    % Save time point i in table of mutiple comparisons
    p{1, 1} = 'p-values';
    
    % Make row names of conditions for comparison
    p(2:end, 1) = stat.names(1:(end-1));
    
    % Make column names of conditions for comparison
    p(1, 2:end) = stat.names(2:end);
            
    % For each comparison
    for j = 1:size(stat.p_indiv, 1)
        % save the p-value in the p-value table
        p{stat.p_indiv(j,1) + 1,...
            stat.p_indiv(j,2)} = stat.p_indiv(j,6);
    end
end

function [stat, p] = statistical_analysis_all(y, g)
%STATISTICAL_ANALYSIS Perform ANOVA to compare the means between conditions
% 
%   Input
%       avg: the structure returned from calc_means
% 
%   Output
%       stat: structure containing outputs from anova1 and multcompare
%       p: table of p-values for pairwise comparisons
% 
%   Overview
%       This function performs statistical analysis on the data.
%       Specifically, it performs one way ANOVA using anova1 and multiple
%       comparisons using Tukey's HSD using multcompare. It returns the
%       outputs from anova1 and multcompare in the structure stat and
%       a table of p-values, p, for pairwise comparison between conditions.
    
    condition = g;
    data = y;
    
    % Initialize a structure for storing the results of the statistical
    % analysis
    stat = struct('p', [],...
                  'tbl', [],...
                  'stats', [],...
                  'p_indiv', [],...
                  'means', [],... 
                  'names', []);
    
    % Perform ANOVA on the intensity data grouped by condition
    [stat.p, stat.tbl, stat.stats, stat.p_indiv, stat.means,...
        stat.names] = stat_test(log(data), condition);
    
    % Initialize variables for making comparison tables
    p = cell(size(stat.names, 1), size(stat.names, 1));

    % Save time point i in table of mutiple comparisons
    p{1, 1} = 'p-values';
    
    % Make row names of conditions for comparison
    p(2:end, 1) = stat.names(1:(end-1));
    
    % Make column names of conditions for comparison
    p(1, 2:end) = stat.names(2:end);
            
    % For each comparison
    for j = 1:size(stat.p_indiv, 1)
        % save the p-value in the p-value table
        p{stat.p_indiv(j,1) + 1,...
            stat.p_indiv(j,2)} = stat.p_indiv(j,6);
    end
end

function [p, tbl, stats, p_indiv, means, names] = stat_test(data,...
    group)
%STAT_TEST Perform ANOVA to compare the means between conditions
% 
%   Input
%       data: data that anova will be performed on
%       group: identifier for data to correctly group it
% 
%   Output
%       p: p-value from the anova
%       tbl: a table returned from anova
%       stats: statistics for mutiple comparison tests
%       p_indiv: pairwise p-values from mutiple comparisons
%       means: estimated means
%       names: names of groups
% 
%   Overview
%       This function performs statistical analysis on the data. 
%       Specifically, it performs one way ANOVA using anova1 and multiple
%       comparisons using Tukey's HSD using multcompare. It returns the
%       outputs from anova1 and multcompare in the structure stat.

    % Perform ANOVA on the data grouped by condition in group
    [p, tbl, stats] = anova1(data, group, 'off');
    
    % Perform pairwise comparisons of data between conditions
    % using Tukey's HSD
    [p_indiv, means, ~, names] = multcompare(stats, 'display', 'off');
end