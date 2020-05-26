function tdat = func_pca(Dat_Pre)
  
    Dreduc_opts = cell(2,2);
    % Options for U-MAP
    Dreduc_opts(:,1) = {'Name', ...
        'Value'};
    % Defaults
    Dreduc_opts(:,2) = {'Rows', ...
        'complete'};
    
    dlg_title = 'Input Channel Arithmatics';
    num_lines = 1;
    vAnswer = inputdlg(Dreduc_opts(:,1),dlg_title,num_lines,Dreduc_opts(:,2));
    if isempty(vAnswer)
        return
    end
    Dreduc_opts(:,2) = vAnswer;

    [coeff,score,~] = pca(Dat_Pre,...
                    Dreduc_opts{1, 2}, ...
                    Dreduc_opts{2, 2});
                
    tdat = score(:, 1:2);
    figure;
    biplot(coeff(:,1:2),'scores',score(1:300,1:2));

