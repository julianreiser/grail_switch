function lw_qq(datamat,factorlevels1,factorlevels2,mauchlytable)

% input parser
p = inputParser;

% plot rmANOVA assumptions using 
% - levene variance homogeinity test,
% - lilliefors normal distribution test (with qq-plot)
% - report mauchly (if added)
%
% input: datamatrix for rmANOVA compuatation, number of factorlevels of
% factor 1, number of factorlevels of factor 2 (& mauchly results)

% parse inputs and set defaults
p.FunctionName  = mfilename;
p.CaseSensitive = false;
p.addRequired('datamat', @ismatrix);
p.addRequired('factorlevels1', @isnumeric);
p.addRequired('factorlevels2', @isnumeric);
p.addRequired('mauchlytable', @istable);

parse(p, datamat, factorlevels1, factorlevels2, mauchlytable);

[levenep, levenestats] = vartestn(p.Results.datamat,'TestType','LeveneAbsolute');

figure('visible','of');
for plotnum = 1:(p.Results.factorlevels1*p.Results.factorlevels2)
    lillieval = lillietest(p.Results.datamat(:,plotnum));
    
    subplot(p.Results.factorlevels1,p.Results.factorlevels2,plotnum)
    qqplot(p.Results.datamat(:,plotnum))
    title(['lilliefors = ' num2str(lillieval)]);
end
sgtitle(['levene: F(' num2str(levenestats.df(1)) ',' num2str(levenestats.df(2)) ')= ' num2str(levenestats.fstat) ', p=' num2str(levenep) ...
         newline 'mauchly: W=' num2str(p.Results.mauchlytable.W) ', Chi(' num2str(p.Results.mauchlytable.DF) ')= ' num2str(p.Results.mauchlytable.ChiStat) ', p=' num2str(p.Results.mauchlytable.pValue)])