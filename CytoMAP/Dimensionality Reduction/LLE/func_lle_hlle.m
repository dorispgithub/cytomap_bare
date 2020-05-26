function tdat = func_lle_hlle(Dat_Pre)
    Dreduc_opts = cell(3,2);
    % Options for U-MAP
    Dreduc_opts(:,1) = {'function',...
        'number of neighbors', ...
        'max embedding dimensionality'};
    % Defaults
    Dreduc_opts(:,2) = {'LLE',...
        '3', ...
        '2'};
    dlg_title = 'Input Channel Arithmatics';
    num_lines = 1;
    vAnswer = inputdlg(Dreduc_opts(:,1),dlg_title,num_lines,Dreduc_opts(:,2));
    if isempty(vAnswer)
        return
    end
    
    Dreduc_opts(:,2) = vAnswer;
    
    Dat_Pre= 10*Dat_Pre + 0.02*rand(size(Dat_Pre));
    if strcmp(Dreduc_opts{1, 2}, 'LLE')
        [Y] = lle(Dat_Pre',str2double(Dreduc_opts{2, 2}), str2double(Dreduc_opts{3, 2}));
        tdat = Y';
    elseif strcmp(Dreduc_opts{1, 2},'HLLE')
        [Y, ~] = HLLE(Dat_Pre',str2double(Dreduc_opts{2, 2}), str2double(Dreduc_opts{3, 2}));
        tdat = Y';
    end