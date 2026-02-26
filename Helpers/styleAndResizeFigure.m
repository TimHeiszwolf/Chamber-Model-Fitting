function styleAndResizeFigure(fig_handle, width_cm, height_cm)
    % STYLEANDRESIZEFIGURE Formats a MATLAB figure for publication or reports.
    %
    %   Syntax:
    %       styleAndResizeFigure(fig_handle)
    %       styleAndResizeFigure(fig_handle, width_cm, height_cm)
    %
    %   Description:
    %       Sets the figure background to white, axis colors to black, and snaps 
    %       the figure to a fixed physical size in centimeters. This ensures 
    %       consistent scaling and resolution when exporting graphics.
    %
    %   Inputs:
    %       fig_handle - Handle to the figure to be formatted.
    %       width_cm   - (Optional) Width of the figure in centimeters
    %       height_cm  - (Optional) Height of the figure in centimeters

    if nargin < 2, width_cm = 25*0.66; end %30, 25*0.66
    if nargin < 3, height_cm = 22*0.66; end %20, 22*0.66

    % --- Configuration ---
    legend_font_size = 11;        
    axis_line_width = 1.5;      % Width for Axis lines and Grid lines (default of MATlab is 0.5)

    %% 1. Consistent Sizing and Scaling
    set(fig_handle, 'Units', 'centimeters');
    pos = get(fig_handle, 'Position');
    
    % Maintain current screen position but enforce fixed dimensions
    set(fig_handle, 'Position', [pos(1), pos(2), width_cm, height_cm]);
    
    % Ensure exported graphics match the on-screen dimensions
    set(fig_handle, 'PaperPositionMode', 'auto');

    %% 2. Color and Axis Formatting
    set(fig_handle, 'Color', 'w');
    
    allAxes = findall(fig_handle, 'Type', 'axes');
    for ax = allAxes' 
        set(ax, ...
            'Color',     'w', ...         
            'XColor',    'k', ...         
            'YColor',    'k', ...         
            'ZColor',    'k', ...         
            'GridColor', [0.5 0.5 0.5], ...
            'LineWidth', axis_line_width); 
        
        ax.Title.Color = 'k';
        ax.XLabel.Color = 'k';
        ax.YLabel.Color = 'k';
        ax.ZLabel.Color = 'k';
    end
    
    %% 3. Legend Formatting
    allLegends = findall(fig_handle, 'Type', 'legend');
    for lgd = allLegends'
        set(lgd, ...
            'Color',     'w', ...         % Solid white background
            'EdgeColor', 'k', ...         % Black border
            'Box',       'on', ...        % Ensure box is visible
            'TextColor', 'k', ...         % Black text.
            'FontSize',  legend_font_size);
    end
end