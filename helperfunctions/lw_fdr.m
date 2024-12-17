function [fdrCorrTable] = lw_fdr(stattable,varargin)

    % Takes a matlab rmANOVA table and adds FDR-correction to the table
    % -----------------------------------------------------------------
    % Inputs:
    % stattable: must be rmANOVA table from matlab-function
    % corrMethod: defines which values are taken to compute the FDR correction
    %   - 0: unadjusted pValues
    %   - 1: GG-adjusted pValues
    %   - 2: HF-adjusted pValues
    

    % input parser
    p = inputParser;

    % default settings
    % default colors for a 2x3 grouped design
    % black, red, blue, black, red, blue
    default_corrMethod = 0;

    % parse inputs and set defaults
    p.FunctionName  = mfilename;
    p.CaseSensitive = false;
    p.addRequired('stattable', @istable);
    p.addParamValue('corrMethod', default_corrMethod, @isnumeric);

    parse(p, stattable, varargin{:});

    pvec = zeros(length(p.Results.stattable.pValue),1);
    logivec = zeros(length(p.Results.stattable.pValue),1);
    critvec = zeros(length(p.Results.stattable.pValue),1);
    
    if p.Results.corrMethod == 0
        corrName = 'uncorrPVal';
        sortmat = zeros(size(p.Results.stattable.pValue(3:2:end,:),1),4);
        sortmat(:,1) = 1:1:size(p.Results.stattable.pValue(3:2:end),1);
        sortmat(:,2) = p.Results.stattable.pValue(3:2:end);
        sortmat = sortrows(sortmat,2,'ascend');
        sortmat(:,3) = linspace(.05,.05/size(p.Results.stattable.pValue(3:2:end),1),size(p.Results.stattable.pValue(3:2:end),1));
        sortmat(find(sortmat(:,2) <= sortmat(:,3)),4) = 1;
        sortmat = sortrows(sortmat,1,'ascend');
    elseif p.Results.corrMethod == 1
        corrName = 'GGpVal';
        sortmat = zeros(size(p.Results.stattable.pValueGG(3:2:end,:),1),4);
        sortmat(:,1) = 1:1:size(p.Results.stattable.pValueGG(3:2:end),1);
        sortmat(:,2) = p.Results.stattable.pValueGG(3:2:end);
        sortmat = sortrows(sortmat,2,'ascend');
        sortmat(:,3) = linspace(.05,.05/size(p.Results.stattable.pValueGG(3:2:end),1),size(p.Results.stattable.pValueGG(3:2:end),1));
        sortmat(find(sortmat(:,2) <= sortmat(:,3)),4) = 1;
        sortmat = sortrows(sortmat,1,'ascend');
    elseif p.Results.corrMethod == 2
        corrName = 'HFpVal';
        sortmat = zeros(size(p.Results.stattable.pValueHF(3:2:end,:),1),4);
        sortmat(:,1) = 1:1:size(p.Results.stattable.pValueHF(3:2:end),1);
        sortmat(:,2) = p.Results.stattable.pValueHF(3:2:end);
        sortmat = sortrows(sortmat,2,'ascend');
        sortmat(:,3) = linspace(.05,.05/size(p.Results.stattable.pValueHF(3:2:end),1),size(p.Results.stattable.pValueHF(3:2:end),1));
        sortmat(find(sortmat(:,2) <= sortmat(:,3)),4) = 1;
        sortmat = sortrows(sortmat,1,'ascend');
    end % if p.Results.corrMethod
        
    pvalcounter = 1;
    for pvalidx = 3:2:length(pvec)
        pvec(pvalidx) = sortmat(pvalcounter,2);
        critvec(pvalidx) = sortmat(pvalcounter,3);
        logivec(pvalidx) = sortmat(pvalcounter,4);
        pvalcounter = pvalcounter + 1;
    end % pvalidx
    
    corrTable = array2table([pvec,critvec,logivec],'VariableNames',{corrName,'critVal','isSignificant'});
    fdrCorrTable = [stattable,corrTable];
    
    
    
   