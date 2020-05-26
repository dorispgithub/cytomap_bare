function tdat = func_phate(Dat_Pre)
    
    ListOption = {'PHATE_EB','PHATE_DLA_TREE','PAHTE_TREE','PHATE_mESC'};
    [Selection, ~] = listdlg('ListString',ListOption, 'SelectionMode','single');
    
    addpath([fileparts(mfilename('fullpath')) filesep 'PHATE']);

    %%
    switch Selection
        case 1 
            % PHATE_EB
            Dat_Pre= 10*Dat_Pre + 0.02*rand(size(Dat_Pre));
            tdat = phate(Dat_Pre, 't', 20);
            figure;
            scatter(tdat(:,1), tdat(:,2), 3, 'filled');
            set(gca,'xticklabel',[]);
            set(gca,'yticklabel',[]);
            axis tight
            xlabel 'PHATE1'
            ylabel 'PHATE2'
            title 'PHATE'
            drawnow
        case 2
            % PHATE_DLA_TREE
            Dat_Pre= 10*Dat_Pre + 0.03*rand(size(Dat_Pre));
            tdat = phate(Dat_Pre, 't', 20, 'gamma', 0);
            figure;
            scatter(tdat(:,1), tdat(:,2), 10, 'filled');
            set(gca,'xticklabel',[]);
            set(gca,'yticklabel',[]);
            axis tight
            xlabel 'PHATE1'
            ylabel 'PHATE2'
            drawnow
        case 3
            % PAHTE_TREE
            Dat_Pre= 10*Dat_Pre + 0.03*rand(size(Dat_Pre));
            tdat = phate(Dat_Pre, 't', 32);
            figure;
            scatter(tdat(:,1), tdat(:,2), 10, 'filled');
        case 4
            % PHATE_mESC
            Dat_Pre= 20*Dat_Pre + 0.04*rand(size(Dat_Pre));
            tdat = phate(Dat_Pre);
            figure;
            scatter(tdat(:,1),tdat(:,2),20,'filled')
            axis tight
            xlabel('PHATE 1')
            ylabel('PHATE 2')
            set(gca,'xticklabel',[])
            set(gca,'yticklabel',[])
            title 'PHATE 2D'
            drawnow
    end