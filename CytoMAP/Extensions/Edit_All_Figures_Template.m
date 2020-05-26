function Edit_All_Figures_Template(app)
    % Pulls the names of all of the currently open figures
    FigN = get(groot, 'Children');
    if ~isempty(FigN)
        FigN = [FigN.Number];
        for i=1:numel(FigN)
            fig = figure(FigN(i));
            ax=gca;
%%%%%%%%%%%
%%%%%%%%%%%
            %% Apply these figure settings to all open figures
            
%             % Delete all of the buttons and controls in the figure
            h = findobj(gcf,'type','UIControl');
            delete(h)

            % Define your axes limits
            ax.XLim= [-450 450];
            ax.YLim= [-450 450];
            
            
            box off
            % Define the colors of your axes labels
            ax.XColor = 'k';
            ax.YColor = 'k';
            
            % Turn the color bar off/on
            colorbar off
            
            % Get rid of the title
            ax.Title = [];

        end
    end    
end
