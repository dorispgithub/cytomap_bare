function tdat = func_tscan(Dat_Pre)
    addpath([fileparts(mfilename('fullpath')) filesep 'TSCAN']);
    % Scale the data to be in the 0-100 range
    IND = range(Dat_Pre)<100;
    Dat_Pre(:, IND) = 100 * Dat_Pre(:, IND);
    % Shift the data so that it is positive
    IND = min(Dat_Pre)<100;
    Dat_Pre(:, IND) = Dat_Pre(:, IND) - min(Dat_Pre(:, IND));
    
    procdata = preprocess(Dat_Pre);
    lpsmclust = exprmclust(procdata');
    tdat(:,1) = lpsmclust.pcareduceres(:,1);
    tdat(:,2) = lpsmclust.pcareduceres(:,2);
    plotmclust(lpsmclust);
end