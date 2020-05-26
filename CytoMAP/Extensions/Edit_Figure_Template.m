function Edit_Figure_Template(app)
    ax = gca;
    
    %% Apply the following figure settings to the currently selected figure
    
   % Delete all of the buttons and controls in the figure
    h = findobj(gcf,'type','UIControl');
    delete(h)

    % Define your axes limits
    ax.XLim= [0 2500];
    ax.YLim= [0 2500];

    % Define the colors of your axes labels
    ax.XColor = 'k';
    ax.YColor = 'k';

    % Turn the color bar off/on
    colorbar off

    % Get rid of the title
    ax.Title = [];



end