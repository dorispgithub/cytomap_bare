classdef Constants

% This class defines some constant parameters that are used across the
% various CytoMAP functions
%
% Last edited by JF 2019/04/02

    properties (Constant = true)
        CURR_VER = "19.04"
        LEGACY_VER = "18_11" % do not change.

        %% Debug Variables
        ellipse_debug = true;

        %% Main Layout constants:
        plot_figure = "New Figure";

        % Neighborhood Analysis
        neighborhood_sect = "Neighborhood Analysis";
        first_scan = "Raster Scan Neighborhoods";
        second_scan = "Cell Centered Neighborhoods";
        distance = "Calculate";

        % Region Analysis
        region_sect = "Region Analysis";
        neigh_to_regions = "Classify Neighborhoods Into Regions";
        use_nn = "Use NN";
        save_nn = "Save NN";
        load_nn = "Load NN";

        % Statistics
        stat_sect = "Plots and Statistics";
        heatmap = "Plot Heatmaps";
        interact = "Interaction Map";
        reg_stats = "Region Statistics";
        t_sne = "t-SNE"; % What is it? (todo)
        cellularity = "Cellularity";
        box_plt = "Box Plot";

        %% Some other Layout
        % Here constants specific to given layouts will be kept

        %% Plot Constants
        min_fig_size = [725 400];
        fig_tag = 'fig';

        % Defining menu constants
        plt_menu_dx = 5;  % Distance between 2 cols in menu
        plt_menu_dy = 5;  % Distance between 2 rows in menu
        plt_menu_x0 = 5;  % Initial displacement in x from corner
        plt_menu_y0 = 5;  % Initial displacement in y from corner

        %% Data Tags
        gate_tag = 'Gate_';  % Tag for gates used in AllCells.
        channel_tag = 'Ch_';  % Tag for channels.
        other_tag = 'Other';  % Anything else tag.
        neigh_tag = 'Nei_';  % Tag for neighborhoods (created by gating on MFI's)

        other_names = {'Region', 'RegionRSN', 'RegionCCN', 'Model'};
        channel_starts = {'Ch_', 'Ch', 'Channel_', 'Channel'};
        ignore__names__internal = {'Neigh Area', 'Neigh Volume', 'Effective Neigh Area', 'Effective Neigh Volume', 'NCells', 'X', 'Y', 'Z'};
        ignore_names = [Constants.ignore__names__internal, Helper.valid_var(Constants.ignore__names__internal)];

        %% File Name related constants
        worksp = ".cym";

        %% Other Constants
        new_model = 'Create New Model';
        model_required_fields = {'name', 'NR', 'userdata', 'type', 'cmap'};
    end
end

