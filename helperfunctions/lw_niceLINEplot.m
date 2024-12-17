function [niceBarplot_figure] = lw_niceLINEplot(chanlocs,channame,meanstruct,ngroups,limit_xaxis,timevec,label_title,label_xaxis,label_yaxis,varargin)

    %% What does it do?
    %--------------------------------------
    % Feed the function a vector containing means (1,2,3; 5,6,7; 2,4,7) and an accourding errorbar vector (0.1,1.2,3.3; 2.3,2.4,3.5; 1.0,2.5,2.1)
    % Define the title ('title'), yaxis label ('RT in ms'), xtick labels {'stand','walk','pert'}
    % Define yaxis limits, e.g. [0,800] to plot a range from 0 to 800
    %
    % optional:
    % cols to define colors in hex-decimal, e.g. {'#272730','#B8357E','#235DB8','#457b9d','#1d3557','#808080'}
    % label_legend to add legend, e.g. {'0degree','18degree','40degree'}
    % locs_legend to define location of legend, e.g. 'northwest'
    %
    % !!! Make sure, that meanvec and semvec have the same dimensions and that all labeling elements have the correct number of entries (3 in the example)
    %
    % OUTPUT ARGUMENTS
    % niceBarplot_figure: a beautiful plot :)
    %
    % EXAMPLE:
    % meanvec = [1,2,3;4,5,6;5,6,7]
    % semvec = [.1,2.3,.2; .1,2.1,1.3; 1,2,3]
    % abclabel_title = 'title'
    % abclabel_yaxis = 'ms'
    % abclabel_xticks = {'a','b','c'}
    % abclabel_legend = {'c','d','e'}
    % limit_yaxis = [0,10]
    %
    % f1 = figure;
    % lw_niceBARplot(meanvec,semvec,abclabel_title,abclabel_yaxis,abclabel_xticks,abclimit_yaxis,'label_legend',abclabel_legend,'locs_legend','northwest')
    %
    %
    % Julian Elias Reiser, 17.05.2020

    % Check if any input args
    if nargin < 6
        error('Not enough input arguments... :)');
        return;
    end

    % input parser
    p = inputParser;

    % default settings
    % default colors for a 2x3 grouped design
    % black, red, blue, black, red, blue
    default_cols = {'#272730','#B8357E','#235DB8','#457b9d','#1d3557','#808080'};
    default_style = {'-',':','-.'};
    default_legend = {};
    default_legendlocs = 'northeast';
    default_polarity = -1;

    % parse inputs and set defaults
    p.FunctionName  = mfilename;
    p.CaseSensitive = false;
    p.addRequired('chanlocs',@isstruct);
    p.addRequired('channame',@iscell);
    p.addRequired('meanstruct',@isdoublep);
    p.addRequired('ngroups', @isnumeric);
    p.addRequired('limit_xaxis', @isvector);
    p.addRequired('timevec', @isvector);
    p.addRequired('label_title', @isstr);
    p.addRequired('label_xaxis', @isstr);
    p.addRequired('label_yaxis', @isstr);
    p.addParamValue('label_legend', default_legend, @iscell);
    p.addParamValue('locs_legend', default_legendlocs, @isstr);
    p.addParamValue('cols', default_cols, @isstring);
    p.addParamValue('polarity', default_polarity, @isnumeric);

    parse(p, chanlocs, channame, meanstruct, ngroups, limit_xaxis, timevec, label_title, label_xaxis, label_yaxis, varargin{:});
    
    rgbcol = [];
    hexcols = p.Results.cols;
    for colconv = 1:length(hexcols)
        hexcoli = hexcols{colconv};
        newcol = [];
        newcol = sscanf(hexcoli(2:end),'%2x%2x%2x',[1 3])/255;
        rgbcol(colconv,:) = [newcol];
    end

    if size(p.Results.channame,2) > 1
        plotchan = zeros(1,size(p.Results.channame,2));
        for chanidx = 1:size(p.Results.channame,2)
            plotchan(chanidx) = DO_channum(p.Results.chanlocs,p.Results.channame{chanidx});
        end
    else [plotchan] = DO_channum(p.Results.chanlocs,p.Results.channame);
    end
            
    plotmat = squeeze(mean(p.Results.meanstruct(:,:,plotchan,:),[1,3]));
    limit_yaxis = [floor(min(plotmat(:))),ceil(max(plotmat(:)))];
    colnum = 1;
    
    %figure;
    
    for linenum = 1:p.Results.ngroups:size(plotmat,1)
        for groupnum = 1:p.Results.ngroups
            ploti = plot(p.Results.timevec,plotmat(linenum + (groupnum-1),:));
            ploti.Color = rgbcol(colnum,:)
            ploti.LineStyle = default_style{groupnum};
            ploti.LineWidth = 2;
            hold on
        end % groupnum
        colnum = colnum +1;
    end % linenum
    hold off
    title(['\fontsize{24}' p.Results.label_title])
    ylabel(['\fontsize{20}' p.Results.label_yaxis])
    xlabel(['\fontsize{20}' p.Results.label_xaxis])
    ylim(limit_yaxis)
    xlim(p.Results.limit_xaxis)
    xticks([p.Results.timevec(1):100:p.Results.timevec(end)])
    set(gcf,'color','w')
    set(gca,'TickDir','out')
    if p.Results.polarity < 0
        set(gca,'YDir','reverse')
    end
    a = gca;
    a.XAxisLocation = 'origin';
    a.YAxisLocation = 'origin';
    a.YGrid = 'on';
    a.Box = 'off';
    a.XAxis.FontSize = 20;
    a.XAxis.Color = '#828282';
    a.XLabel.Color = 'black'
    a.YAxis.FontSize = 20;
    a.YAxis.Color = '#828282';
    a.YLabel.Color = 'black'
    if not(isempty(p.Results.label_legend))
        legend(p.Results.label_legend,'Location',p.Results.locs_legend)
    end
