%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Multidimensional Analysis Pipeline for Histocytometry (CytoMAP)
%%%
%%% This program is meant to make advanced data analytic techniques
%%% accesible for Histo-Cytometry data. The goal is to take established
%%% analytical techniques, such as neural network based clustering
%%% algorithms, and package them in a user friendly way, such that researchers
%%% can use them to explore complex datasets.
%%%
%%% Any dataset structured as a .csv, with at least x-y positional data can
%%% be loaded into this program. The data should be formatted as Columns of parameters and
%%% rows of individual cells/objects.
%%%
%%% Authored by: 
%%% Dr. Caleb Stoltzfus - lead author
%%% University of Washington, Department of Immunology
%%% Jakub Filipek: 
%%% University of Washington, Paul G. Allen School of Computer Science and Engineering
%%% technical assistant from October 2018 to August 2019
%%% Yajun Wu: 
%%% University of Washington, Mechanical Engineering
%%% technical assistant from September 2019
%%% University of Washington, School of Medicine
%%% Immunology Department, Gerner Lab
%%% Version Beta - 07-17-2019
%%% Written in MATLAB 2018b/2019a
%%%
% % % MIT License
% % % 
% % % Copyright (c) 2019 Caleb Stoltzfus
% % % 
% % % Permission is hereby granted, free of charge, to any person obtaining a copy
% % % of this software and associated documentation files (the "Software"), to deal
% % % in the Software without restriction, including without limitation the rights
% % % to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% % % copies of the Software, and to permit persons to whom the Software is
% % % furnished to do so, subject to the following conditions:
% % % 
% % % The above copyright notice and this permission notice shall be included in all
% % % copies or substantial portions of the Software.
% % % 
% % % THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% % % IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% % % FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% % % AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% % % LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% % % OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% % % SOFTWARE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef CytoMAP < handle
    properties
        % These are my global variables that I will use in the program
        % This initializes the properties of the "app" object
        GUI;                     % Graphics handles
        GUIOPTS;                 % Options for the window interfaces
        Logo;                    % Make a handle for the logo object
        figOpts;                 % Make a handle for the options in the figure window
        PTable;
        PhList;
        data;
        DataN;
        bgclr;
        txtclr;
        rwindowRSN = 50;
        rwindowCCN = 50;
        RegAlg = 'Davies Bouldin (Default)';
        NReg = 2;
        pTH = 0.01;
        map = struct;
        net = struct;
        weights = [];
        CellInfo;
        points = struct;
        polygons = struct;
        PlotMenu;
        sortindex;
        NewGate;
        TMP; % a place to store temporary global variables
        DefaultGUI

    end

    
    methods
        % Build the CytoMAP graphical user interface
        function app = CytoMAP
            % APP Set the current Path
            CurrentPath = fileparts(mfilename('fullpath'));
            addpath(CurrentPath);
            % Import the various Class Definitions/Groups of functions
            import Plotting;
            import Analysis;
            import Statistics;
            import Helper;
            import tree;
            import IO;
            import IO_wsp;

            % This is the "constructor" for the class
            % It runs when an object of this class is created
            %% Create a UI figure window and initialize some variables
            % Set the default interface options
            app.GUIOPTS = struct;
            % Color Scheme
            app.GUIOPTS.Scheme = 'Light';
            % background main windo color            
            app.GUIOPTS.bgclr = [1,1,1];
            % background window color
            app.GUIOPTS.mainbgclr = [0.99,0.99,0.99];
            % button background color
            app.GUIOPTS.btnbgclr = [1,1,1];
            % text color
            app.GUIOPTS.txtclr = [0,0,0];
            % text font
            app.GUIOPTS.FontName = 'Helvetica';%'Arial Narrow'; % Use listfonts command to see what fonts are available
            % button height
            app.GUIOPTS.btnheight = 20;
            % Font weight
            app.GUIOPTS.FontWeight = 'normal';
            % Store available fonts
            app.GUIOPTS.Fonts = listfonts;
            % Store available sizes
            app.GUIOPTS.Alpha = {0.75,1,1.25,1.5};
            app.GUIOPTS.AlphaF = {8,11,16,20};
            app.GUIOPTS.AlphaL ={'Small','Medium','Large','Larger'};
            app.GUIOPTS.AlphaI = 2;
            % store a default continious resampled red-blue colormap
            cmap = redbluecmap;
            cmap(2,3) = cmap(2,3)-0.05;
            m = 120; % turn from 12 to 120 colorpoints
            upsmp = round(m/size(cmap,1));
            r = abs(interp(cmap(:,1),upsmp));
            g = abs(interp(cmap(:,2),upsmp));
            b = abs(interp(cmap(:,3),upsmp));
            r = r(1:(end-10));
            g = g(1:(end-10));
            b = b(1:(end-10));
            app.GUIOPTS.redbluecmap = [r, g, b];
% %             app.GUIOPTS.redbluecmap = redbluecmap;

            % Create the Main User interface component handles
            app.GUI = struct;
%             app.GUI.Main
            app.GUI.Dropdowns = struct; 
            app.GUI.Labels = struct;
            app.GUI.Buttons = struct;
            app.GUI.Options = app.GUIOPTS; % this should replace app.GUIOPTS

            % Create the main window
            Size = [255 700];      % width = 255;    height = 700;
            app.GUI.Main  = uifigure;
            app.GUI.Main.CloseRequestFcn = @(~, event) Helper.func_exit(app);
            app.GUI.Main.Position = [0 0 Size];
            app.GUI.Main.NumberTitle = 'on';
            app.GUI.Main.Name = 'CytoMAP';
            app.GUI.Main.Resize = 'off';
            app.GUI.Main.Color = app.GUIOPTS.mainbgclr;

% % %             app.bgclr = 'w';
% % %             app.txtclr = 'k';
            app.map.jet = jet;
            app.figOpts.bgclr = struct;
            app.figOpts.txtclr = struct;
            
            addpath(CurrentPath)
            % add path for third party functions
            addpath(fullfile([CurrentPath filesep '3rdPartyFunctions']));  
            app.figOpts.icons = Helper.Import_Icons();
            set(groot, 'DefaultTextInterpreter', 'none')
            set(groot, 'DefaultLegendInterpreter', 'none')
            %% Put in an CytoMAP label
            app.GUI.Labels.Cyto = uilabel(app.GUI.Main);
            app.GUI.Labels.Cyto.Text = 'CytoMAP';
            app.GUI.Labels.Cyto.Position = [0 (700-580), 255, 80];
            Helper.func_SetCLR(app, app.GUI.Labels.Cyto, 'label')
            app.GUI.Labels.Cyto.FontSize = 60;
            app.GUI.Labels.Cyto.FontName = 'Arial Narrow';
            app.GUI.Labels.Cyto.HorizontalAlignment = 'center';
            
            % put in subtext
            app.GUI.Labels.subtxt = uilabel(app.GUI.Main);
            app.GUI.Labels.subtxt.Text = sprintf(['Multidimensional Analysis Pipeline for Histocytometry\n\n'...
                                'Dr. Caleb Stoltzfus, Jakub Filipek, Yajun Wu\n' ...
                                'University of Washington, Gerner Lab\n' ...
                                'Version 1.3']);
            app.GUI.Labels.subtxt.FontWeight = 'normal';                
            app.GUI.Labels.subtxt.Position = [0 (700-655), 255, 70];
            Helper.func_SetCLR(app, app.GUI.Labels.subtxt, 'label')
            app.GUI.Labels.subtxt.FontSize = 10;
            app.GUI.Labels.subtxt.HorizontalAlignment = 'center';
            
            app.GUI.Buttons.LabWbs = uibutton(app.GUI.Main);
            app.GUI.Buttons.LabWbs.Text= 'Our Lab Website';
            app.GUI.Buttons.LabWbs.ButtonPushedFcn =@(button,event) web('http://depts.washington.edu/myglab/', '-browser');
            app.GUI.Buttons.LabWbs.Position = [10 (700-690), 117.5, 30];
            Helper.func_SetCLR(app, app.GUI.Buttons.LabWbs, 'button')
            app.GUI.Buttons.LabWbs.FontColor = 'b';
            
            app.GUI.Buttons.WikiWbs = uibutton(app.GUI.Main);
            app.GUI.Buttons.WikiWbs.Text= 'CytoMAP Wiki';
            app.GUI.Buttons.WikiWbs.ButtonPushedFcn =@(button,event) web('https://gitlab.com/gernerlab/cytomap/wikis/home/', '-browser');
            app.GUI.Buttons.WikiWbs.Position = [10+117.5 (700-690), 117.5, 30];
            Helper.func_SetCLR(app, app.GUI.Buttons.WikiWbs, 'button')
            app.GUI.Buttons.WikiWbs.FontColor = 'b';

            
            %% Create File Menu bar options
            % Create Menu bar
            app.GUI.Dropdowns.FileMenu.Main = uimenu(app.GUI.Main);
            app.GUI.Dropdowns.FileMenu.Main.Text = 'File';

            % Create a load data option
            app.GUI.Dropdowns.FileMenu.load = uimenu(app.GUI.Dropdowns.FileMenu.Main);
            app.GUI.Dropdowns.FileMenu.load.MenuSelectedFcn = @(btn,event) IO.func_load(app, 0);
            app.GUI.Dropdowns.FileMenu.load.Text = 'Load .csv or .mat';

            % Create a load multiple samples
            app.GUI.Dropdowns.FileMenu.loadmulti = uimenu(app.GUI.Dropdowns.FileMenu.Main);
            app.GUI.Dropdowns.FileMenu.loadmulti.MenuSelectedFcn = @(btn,event) IO.func_load(app, 1);
            app.GUI.Dropdowns.FileMenu.loadmulti.Text = 'Import Multiple Samples';

            % Create import .wsp info
            app.GUI.Dropdowns.FileMenu.loadwsp = uimenu(app.GUI.Dropdowns.FileMenu.Main);
            app.GUI.Dropdowns.FileMenu.loadwsp.MenuSelectedFcn = @(btn,event) IO_wsp.func_WSPload(app);
            app.GUI.Dropdowns.FileMenu.loadwsp.Text = 'Import data from wsp';

            % Create a save data option
            app.GUI.Dropdowns.FileMenu.savewrksp = uimenu(app.GUI.Dropdowns.FileMenu.Main);
            app.GUI.Dropdowns.FileMenu.savewrksp.MenuSelectedFcn = @(btn,event) IO.func_save(app);
            app.GUI.Dropdowns.FileMenu.savewrksp.Text = 'Save Workspace';
            
            % Create an export data option
            app.GUI.Dropdowns.FileMenu.ExportCSV = uimenu(app.GUI.Dropdowns.FileMenu.Main);
            app.GUI.Dropdowns.FileMenu.ExportCSV.MenuSelectedFcn = @(btn,event) IO.func_ExportCSV(app);
            app.GUI.Dropdowns.FileMenu.ExportCSV.Text = 'Export Full Data Tables as .csv';

            % Create an export data option
            app.GUI.Dropdowns.FileMenu.ExportPrsm = uimenu(app.GUI.Dropdowns.FileMenu.Main);
            app.GUI.Dropdowns.FileMenu.ExportPrsm.MenuSelectedFcn = @(btn,event) IO.func_ExportPrism(app);
            app.GUI.Dropdowns.FileMenu.ExportPrsm.Text = 'Export Data for Prism';

% % % % % % % % % DEPRECATED - need to clean up some of the associated functions
% % %             % Create import .fcs file 
% % %             btn = uimenu(app.GUI.Dropdowns.FileMenu.Main);
% % %             btn.MenuSelectedFcn = @(btn,event) IO.func_FCSLoad(app);
% % %             btn.Text = 'Import data from fcs file';

% % %             % Import .ims confocal image
% % %             btn = uimenu(app.GUI.Dropdowns.FileMenu.Main);
% % %             btn.MenuSelectedFcn = @(btn,event) IO.func_IMSload(app);
% % %             btn.Text = 'Import .ims';

            %% Create View Menu bar options
            app.GUI.Dropdowns.ViewMenu.Main = uimenu(app.GUI.Main);
            app.GUI.Dropdowns.ViewMenu.Main.Text = 'View';
            
            % Create a view data option
            app.GUI.Dropdowns.ViewMenu.tbl = uimenu(app.GUI.Dropdowns.ViewMenu.Main);
            app.GUI.Dropdowns.ViewMenu.tbl.MenuSelectedFcn = @(btn,event) Helper.func_dataWDW(app);
            app.GUI.Dropdowns.ViewMenu.tbl.Text = 'View data table';

            % Create a customize interface
            app.GUI.Dropdowns.ViewMenu.inter = uimenu(app.GUI.Dropdowns.ViewMenu.Main);
            app.GUI.Dropdowns.ViewMenu.inter.MenuSelectedFcn = @(btn,event) Helper.func_EditGUIInterface(app);
            app.GUI.Dropdowns.ViewMenu.inter.Text = 'Customize Interface';

            %% Create Edit Menu bar options
            EditMenu = uimenu(app.GUI.Main);
            EditMenu.Text = 'Edit';
            % Create a view data option
            btn = uimenu(EditMenu);
            btn.MenuSelectedFcn = @(btn,event) Plt_Helper.func_EditColormap(app);
            btn.Text = 'Edit Region Colors';
            
            % Create a view data option
            btn = uimenu(EditMenu);
            btn.MenuSelectedFcn = @(btn,event) Plt_Helper.func_EditChannels(app);
            btn.Text = 'Edit Channels';
            
            % Create a view data option
            btn = uimenu(EditMenu);
            btn.MenuSelectedFcn = @(btn,event) Plt_Helper.func_EditNeighbors(app);
            btn.Text = 'Edit Neighborhoods';
            
            %% Create extensions menu bar option
            ExtenMenu = uimenu(app.GUI.Main);
            ExtenMenu.Text = 'Extensions';
            Helper.func_RunExtensions(app, [], CurrentPath, ExtenMenu);

            %% Create Phenotype table 
            %%% deprecated -need to remove dependencies on PhList and DataN
            app.PhList = uifigure('Name', 'Data');
            app.PhList.CloseRequestFcn = @(~, event) Helper.func_dataWDW(app);
            app.PhList.Color = app.GUIOPTS.bgclr;
            app.PhList.Position(3) = 670;
            app.PhList.Visible = 'off';
            app.DataN = uidropdown(app.PhList);
            app.DataN.Items = {'No Data Loaded'}; %Initialize with no data loaded
            app.DataN.Position = [10 10 235 30];
            app.DataN.BackgroundColor = app.GUIOPTS.bgclr;
            
            %% Top Group of Main Functions
            top = 20; % Define the position of this group
            
            % New figure button
            app.GUI.Buttons.nfig = uibutton(app.GUI.Main,'push', 'ButtonPushedFcn', @(btn,event) Plotting.func_newfig(app));
            app.GUI.Buttons.nfig.Position = [10, (700-top-30), 235, 40];
            app.GUI.Buttons.nfig.Text = 'New Figure';
            Helper.func_SetCLR(app, app.GUI.Buttons.nfig, 'button')
            app.GUI.Buttons.nfig.FontSize = round(app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI}*1.25);

            %% Group of Cell specific functions
            top = top+40+25+5; % Define the position of this group
            
            app.GUI.Labels.btngrp4 = uilabel(app.GUI.Main); 
            app.GUI.Labels.btngrp4.Text = 'Cells';
            Helper.func_SetCLR(app, app.GUI.Labels.btngrp4, 'label')
            app.GUI.Labels.btngrp4.FontWeight = 'bold';
            app.GUI.Labels.btngrp4.FontSize = round(app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI}*1.5);
            app.GUI.Labels.btngrp4.Position = [0 (700-top) 255 30];
            app.GUI.Labels.btngrp4.HorizontalAlignment = 'center';

            % Calculate distance
            app.GUI.Buttons.calcdistance = uibutton(app.GUI.Main,'push', 'ButtonPushedFcn', @(btn,event) Analysis.dist_func(app));
            app.GUI.Buttons.calcdistance.Position = [10 (700-top-30), 235, 30];
            app.GUI.Buttons.calcdistance.Text = 'Calculate Distance';
            Helper.func_SetCLR(app, app.GUI.Buttons.calcdistance, 'button')
            
            % Clustering
            app.GUI.Buttons.clustercells = uibutton(app.GUI.Main,'push', 'ButtonPushedFcn', @(btn,event) Analysis.TNN_func(app, 'Cells'));
            app.GUI.Buttons.clustercells.Position = [10 (700-top-2*30), 235, 30];
            Helper.func_SetCLR(app, app.GUI.Buttons.clustercells, 'button')
            app.GUI.Buttons.clustercells.Text = 'Cluster Cells Into Phenotypes';
            
            % Plot cellularity
            app.GUI.Buttons.cellularity = uibutton(app.GUI.Main,'push', 'ButtonPushedFcn', @(btn,event) Statistics.cellularity_func(app));
            app.GUI.Buttons.cellularity.Position = [10 (700-top-3*30), 117.5, 30];
            app.GUI.Buttons.cellularity.Text = 'Cellularity';
            Helper.func_SetCLR(app, app.GUI.Buttons.cellularity, 'button')
            
            % Generate Box Plots with p-value for population comparisons
            app.GUI.Buttons.boxplt = uibutton(app.GUI.Main,'push', 'ButtonPushedFcn', @(btn,event) Statistics.boxplt_func(app));
            app.GUI.Buttons.boxplt.Position = [127.5 (700-top-3*30), 117.5, 30];
            app.GUI.Buttons.boxplt.Text = 'Compare MFI';
            Helper.func_SetCLR(app, app.GUI.Buttons.boxplt, 'button')


            %% Group of Neighborhood Specific Functions
            top = top+2*60+15; % Define the position of this group
            
            app.GUI.Labels.btngrp1 = uilabel(app.GUI.Main); 
            app.GUI.Labels.btngrp1.Text = 'Neighborhoods';
            Helper.func_SetCLR(app, app.GUI.Labels.btngrp1, 'label')
            app.GUI.Labels.btngrp1.FontWeight = 'bold';
            app.GUI.Labels.btngrp1.FontSize = round(app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI}*1.5);
            app.GUI.Labels.btngrp1.Position = [0 (700-top) 255 30];
            app.GUI.Labels.btngrp1.HorizontalAlignment = 'center';
            
            % Make a Generate Neighborhoods button
            app.GUI.Buttons.genneigh = uibutton(app.GUI.Main,'push', 'ButtonPushedFcn', @(btn,event) Analysis.gen_neigh_func(app));
            app.GUI.Buttons.genneigh.Position = [10 (700-top-30), 235, 30];
            app.GUI.Buttons.genneigh.Text = 'Define Neighborhoods';
            Helper.func_SetCLR(app, app.GUI.Buttons.genneigh, 'button')
            
            % Clustering
            app.GUI.Buttons.cluster = uibutton(app.GUI.Main,'push', 'ButtonPushedFcn', @(btn,event) Analysis.TNN_func(app, 'Neigh'));
            app.GUI.Buttons.cluster.Position = [10 (700-top-2*30), 235, 30];
            Helper.func_SetCLR(app, app.GUI.Buttons.cluster, 'button')
            app.GUI.Buttons.cluster.Text = 'Cluster Neighborhoods Into Regions';
            
            % Generate Pseudo-Space
            app.GUI.Buttons.pseudospace = uibutton(app.GUI.Main,'push', 'ButtonPushedFcn', @(btn,event) Statistics.func_Pseudospace(app));
            app.GUI.Buttons.pseudospace.Position = [10 (700-top-3*30), 117.5, 30];
            app.GUI.Buttons.pseudospace.Text = 'Pseudo-Space';
            Helper.func_SetCLR(app, app.GUI.Buttons.pseudospace, 'button')
            
            % cell-celll interaction maps
            app.GUI.Buttons.correlation = uibutton(app.GUI.Main,'push', 'ButtonPushedFcn', @(btn,event) Statistics.CellInteraction_func(app));
            app.GUI.Buttons.correlation.Position = [127.5 (700-top-3*30), 117.5, 30];
            app.GUI.Buttons.correlation.Text = 'Cell-Cell Correlation';
            Helper.func_SetCLR(app, app.GUI.Buttons.correlation, 'button')
            
            %% Group of region Specific Functions
            top = top+2*60+15; % Define the position of this group

            app.GUI.Labels.btngrp2 = uilabel(app.GUI.Main); 
            app.GUI.Labels.btngrp2.Text = 'Regions';
            Helper.func_SetCLR(app, app.GUI.Labels.btngrp2, 'label')
            app.GUI.Labels.btngrp2.FontWeight = 'bold';
            app.GUI.Labels.btngrp2.FontSize = round(app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI}*1.5);
            app.GUI.Labels.btngrp2.Position = [0 (700-top) 255 30];
            app.GUI.Labels.btngrp2.HorizontalAlignment = 'center';
            
            % Generate Heatmap
            app.GUI.Buttons.neighheatmaps = uibutton(app.GUI.Main,'push', 'ButtonPushedFcn', @(btn,event) Statistics.heatmap_func(app));
            app.GUI.Buttons.neighheatmaps.Position = [10 (700-top-1*30), 117.5, 30];
            app.GUI.Buttons.neighheatmaps.Text = 'Region Heatmaps';
            Helper.func_SetCLR(app, app.GUI.Buttons.neighheatmaps, 'button')
            
            % Generate Region Statistics
            app.GUI.Buttons.regstats = uibutton(app.GUI.Main,'push', 'ButtonPushedFcn', @(btn,event) Statistics.reg_stats_func(app));
            app.GUI.Buttons.regstats.Position = [127.5 (700-top-1*30), 117.5, 30];
            app.GUI.Buttons.regstats.Text = 'Region Statistics';
            Helper.func_SetCLR(app, app.GUI.Buttons.regstats, 'button')
            
            % Generate Interaction Map
            app.GUI.Buttons.reginteract = uibutton(app.GUI.Main,'push', 'ButtonPushedFcn', @(btn,event) Statistics.interact_func(app));
            app.GUI.Buttons.reginteract.Position = [10 (700-top-2*30), 117.5, 30];
            app.GUI.Buttons.reginteract.Text = 'Region Interactions';
            Helper.func_SetCLR(app, app.GUI.Buttons.reginteract, 'button')
            
            % Edit Regions
            app.GUI.Buttons.editreg = uibutton(app.GUI.Main,'push', 'ButtonPushedFcn', @(btn,event) Plt_Helper.func_EditColormap(app));
            app.GUI.Buttons.editreg.Position = [127.5 (700-top-2*30), 117.5, 30];
            app.GUI.Buttons.editreg.Text = 'Edit Regions';
            Helper.func_SetCLR(app, app.GUI.Buttons.editreg, 'button')
                        
            %% Group of Advanced Functions
            top = top+3*30+15; % Define the position of this group

            app.GUI.Labels.btngrp3 = uilabel(app.GUI.Main); 
            app.GUI.Labels.btngrp3.Text = 'Advanced';
            Helper.func_SetCLR(app, app.GUI.Labels.btngrp3, 'label')
            app.GUI.Labels.btngrp3.FontWeight = 'bold';
            app.GUI.Labels.btngrp3.FontSize = round(app.GUIOPTS.AlphaF{app.GUIOPTS.AlphaI}*1.5);
            app.GUI.Labels.btngrp3.Position = [0 (700-top) 255 30];
            app.GUI.Labels.btngrp3.HorizontalAlignment = 'center';
            
            % Generate Dimensionality Reductionplots like t-SNE, UMAP, etc
            app.GUI.Buttons.dimreduc = uibutton(app.GUI.Main,'push', 'ButtonPushedFcn', @(btn,event) Statistics.ReduceDims_func(app));
            app.GUI.Buttons.dimreduc.Position = [10 (700-top-1*30), 235, 30];
            app.GUI.Buttons.dimreduc.Text = 'Dimensionality Reduction';
            Helper.func_SetCLR(app, app.GUI.Buttons.dimreduc, 'button')
            
        end % End of Building GUI
    end % End of methods
end % End of class definition