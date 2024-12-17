function [etaTable] = lw_eta(stattable,varargin)

    % calculates adjusted eta squared
    % -----------------------------------------------------------------
    % Inputs:
    % stattable: must be rmANOVA table from matlab-function
    

    % input parser
    p = inputParser;

    % parse inputs and set defaults
    p.FunctionName  = mfilename;
    p.CaseSensitive = false;
    p.addRequired('stattable', @istable);

    parse(p, stattable, varargin{:});

    calctable = table2array(p.Results.stattable);
    etatable = zeros(size(calctable,1),1);
    
    for etaidx = 1:2:size(calctable,1)
        etatable(etaidx) = (calctable(etaidx,4)*calctable(etaidx,2)) / (calctable(etaidx,4)*calctable(etaidx,2) + calctable(etaidx+1,2)) - ...
                           (1 - (calctable(etaidx,4)*calctable(etaidx,2)) / (calctable(etaidx,4)*calctable(etaidx,2) + calctable(etaidx+1,2))) * (calctable(etaidx,2)/calctable(etaidx+1,2));
    end % etaidx
    
    corrTable = array2table([etatable],'VariableNames',{'eta'});
    etaTable = [p.Results.stattable,corrTable];