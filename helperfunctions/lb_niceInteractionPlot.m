function [niceBarplot_figure] = lw_niceBARplot(meanvec,semvec,label_title,label_yaxis,label_xticks,limit_yaxis,varargin)

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
    default_cols = {'#272730','#B8357E','#235DB8','#457b9d','#1d3557','#808080','#272730','#B8357E','#235DB8'};
    default_legend = {};
    default_legendlocs = 'best';
    default_spacing = 0.1;

    % parse inputs and set defaults
    p.FunctionName  = mfilename;
    p.CaseSensitive = false;
    p.addRequired('meanvec', @ismatrix);
    p.addRequired('semvec', @ismatrix);
    p.addRequired('label_title', @isstr);
    p.addRequired('label_yaxis', @isstr);
    p.addRequired('label_xticks', @iscell);
    p.addRequired('limit_yaxis', @isvector);
    p.addParamValue('label_legend', default_legend, @iscell);
    p.addParamValue('locs_legend', default_legendlocs, @isstr);
    p.addParamValue('cols', default_cols, @isstring);
    p.addParamValue('spacing',default_spacing,@isnumeric)

    parse(p, meanvec, semvec, label_title, label_yaxis, label_xticks, limit_yaxis, varargin{:});

    rgbcol = [];
    hexcols = p.Results.cols;
    for colconv = 1:length(hexcols)
        hexcoli = hexcols{colconv};
        newcol = [];
        newcol = sscanf(hexcoli(2:end),'%2x%2x%2x',[1 3])/255;
        rgbcol(colconv,:) = [newcol];
    end

    ba = plot(p.Results.meanvec,'LineWidth',3,'Marker','o','MarkerSize',10);
    % insert errorbars
    hold on
        % Find the number of groups and the number of bars in each group
        [ngroups,nbars] = size(p.Results.meanvec);
        % Get the x coordinate of the bars
        x = nan(nbars, ngroups);
        x_offset = (0:p.Results.spacing/(nbars-1):p.Results.spacing) - mean(0:p.Results.spacing/(nbars-1):p.Results.spacing)
        for i = 1:nbars
            ba(i).XData = ba(i).XData + x_offset(i)
            x(i,:) = ba(i).XData;
        end
        % Plot the errorbars
        e = errorbar(x',p.Results.meanvec,p.Results.semvec,'k','linestyle','none','LineWidth',3);
    hold off
    barcolcount = 1;
    for barcol = 1:nbars
        set(ba(barcol),'color',rgbcol(barcol,:),'LineStyle','-','MarkerEdgeColor','none','MarkerFaceColor',rgbcol(barcol,:))
        set(e(barcol),'color',rgbcol(barcol,:))
    end % barcol
    title(['\fontsize{24}' p.Results.label_title])
    ylabel(['\fontsize{20}' p.Results.label_yaxis])
    ylim(p.Results.limit_yaxis)
    yticks([p.Results.limit_yaxis(1):(p.Results.limit_yaxis(2)-p.Results.limit_yaxis(1))/10:p.Results.limit_yaxis(2)])
    set(gca,'XTick',x(1,:))
    set(gca,'XTickLabel',p.Results.label_xticks,'fontsize',15)
    xlim([0.5,size(x,2)+0.5])
    set(gcf,'color','w')
    set(gca,'TickDir','out')
    a = gca;
    a.YGrid = 'on';
    % delete mirrored x/y ticks
    set(a,'box','off','color','none')
    b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
    axes(a)
    if not(isempty(p.Results.label_legend))
        legend(p.Results.label_legend,'Location',p.Results.locs_legend)
    end

