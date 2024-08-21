classdef Plotter < handle
    %Utility class for the plots.
    
    properties
        LineWidth = 2;
        MarkerSize = 10;
        HandleVisibility = 'on';
        FontSize = 14;
        AxesFontSize = 10;
        Interpreter = 'latex';
        XLabel = 'Time [s]';
        YLabel = '';
        Hold = 'on';
        Box = 'on';
        Grid = 'on';
        PBaspect = [2.4, 1, 1];
        LegendColumns = 1;
        LegendLocation = 'best';
        %Change this to change the length of the lines in the legend
        LegendItemTokenSize = [20, 12];
        LogScale = false;
        YScale   = 'linear'; % or 'log'
        XScale   = 'linear'; % or 'log'
        DrawNow  = false;
        Colors   = [];
    end
    
    methods
        function obj = Plotter()
            obj.Colors = colororder;
        end

        % Get a value of the colors
        function c = GetColor(obj, idx)
            lColors = length(obj.Colors);
            if idx > lColors
                idx = idx - (fix(idx/lColors))*lColors;
                if idx == 0
                    idx = lColors;
                end
            end
            c = obj.Colors(idx, :);
        end

        function SetColors(obj, colors)
            [~, s2] = size(colors);
            if s2 ~= 3
                error("The colors must be a N x 3 matrix where each row is a different color.");
            end
            obj.Colors = colors;
        end

        function y = formatdata(obj, x)
            % Function that formats a given data to be used by the plot function. 
            % In particular, the dimension of x is reduced to the minimum one and organized in columns
            arguments (Input)
                obj (1, 1) Plotter
                x double
            end
            
            % Squeeze the input and assign it to y
            y = squeeze(x);

            [s1, s2] = size(y);
            if s1 > s2
                y = y';
            end


        end
        
        %Custom fill function
        function fill(obj, XData, YData, options)
            arguments (Input)
               obj (1, 1) Plotter
               XData = []
               YData = []
               options.LineSpec = '-'
               options.LineWidth = obj.LineWidth
               options.FaceColor = ''
               options.EdgeColor = ''
               options.MarkerSize = obj.MarkerSize
               options.HandleVisibility (:, 1) = obj.HandleVisibility
               options.FontSize = obj.FontSize
               options.Interpreter = obj.Interpreter
               options.XLabel = obj.XLabel
               options.YLabel = obj.YLabel
               options.XTicks = []
               options.YTicks = []
               options.XLim = []
               options.YLim = []
               options.Hold = obj.Hold
               options.Box = obj.Box
               options.Grid = obj.Grid
               options.PBaspect = obj.PBaspect
               options.DisplayName = ''
               options.Legend = 'on'
               options.LegendLabels = {}
               options.LegendOrientation = 'horizontal'
               options.LegendLocation = obj.LegendLocation
               options.LegendColumns = obj.LegendColumns
               options.YScale = obj.YScale
               options.XScale = obj.XScale
               options.Drawnow = obj.DrawNow
               options.LegendItemTokenSize = obj.LegendItemTokenSize
               options.AxesFontSize = obj.AxesFontSize
               options.FaceAlpha  = 1;
               options.LineStyle  = "none";
            end

          
           if ~isempty(XData) && ~isempty(YData)
               pl = fill(XData, YData, options.LineSpec, 'LineWidth', options.LineWidth, 'MarkerSize', options.MarkerSize, 'FaceAlpha', options.FaceAlpha, "LineStyle", options.LineStyle);
               set(gca, 'XScale', options.XScale);
               set(gca, 'YScale', options.YScale);
           end
           % Set the diplayname
           DisplayNames = options.DisplayName;
           if ~iscell(DisplayNames)
               DisplayNames = {DisplayNames};
           end
           [s1, ~] = size(DisplayNames);
           if s1 == 1%Enforce a column cell array
                DisplayNames = DisplayNames';
           end
           if exist("pl", "var")
                set(pl, {'DisplayName'}, DisplayNames);
           end
           % Set the handle visibility
           HandleVisibilities = options.HandleVisibility;
           if ~iscell(HandleVisibilities)
               HandleVisibilities = {HandleVisibilities};
           end
           [s1, ~] = size(HandleVisibilities);
           if s1 == 1%Enforce a column cell array
                HandleVisibilities = HandleVisibilities';
           end
           if exist("pl", "var")
                set(pl, {'HandleVisibility'}, HandleVisibilities);
           end
           %
           if ~isempty(char(options.FaceColor))
                pl.FaceColor = options.FaceColor;
           else
                colororder(gca(), obj.Colors);% If no color is specified use the color order
           end

           if ~isempty(char(options.EdgeColor))
                pl.EdgeColor = options.EdgeColor;
           else
                colororder(gca(), obj.Colors);% If no color is specified use the color order
           end

           %Set the font of the axes
           ax = gca();
           ax.FontSize = options.AxesFontSize;
           ax.TickLabelInterpreter = options.Interpreter;

           %Set the x- and y-ticks
           if ~isempty(options.XTicks)
               xticks(options.XTicks);
           end
           if ~isempty(options.YTicks)
               yticks(options.YTicks);
           end

           %Set the x- and y-lim
           if ~isempty(options.XLim)
               xlim(options.XLim);
           end
           if ~isempty(options.YLim)
               ylim(options.YLim);
           end

           %Make the plot
           xlabel(options.XLabel, 'Interpreter', options.Interpreter, 'FontSize', options.FontSize);
           ylabel(options.YLabel, 'Interpreter', options.Interpreter, 'FontSize', options.FontSize);
           pbaspect(options.PBaspect);
           grid(options.Grid);
           hold(options.Hold);
           box(options.Box);

           
           if strcmp(options.Legend, 'on')
                lg = legend(options.LegendLabels, ...
                        'Interpreter', options.Interpreter, ...
                        'FontSize', options.FontSize, ...
                        'Orientation', options.LegendOrientation, ...
                        'Location', options.LegendLocation, ...
                        'NumColumns', options.LegendColumns);
                lg.ItemTokenSize = options.LegendItemTokenSize;
           end

           if options.Drawnow
               drawnow;
           end
        end


        %Custom plot function
        function plot(obj, XData, YData, options)
            arguments (Input)
               obj (1, 1) Plotter
               XData = []
               YData = []
               options.LineSpec = '-'
               options.LineWidth = obj.LineWidth
               options.Color = ''
               options.MarkerSize = obj.MarkerSize
               options.HandleVisibility (:, 1) = obj.HandleVisibility
               options.FontSize = obj.FontSize
               options.Interpreter = obj.Interpreter
               options.XLabel = obj.XLabel
               options.YLabel = obj.YLabel
               options.XTicks = []
               options.YTicks = []
               options.XLim = []
               options.YLim = []
               options.Hold = obj.Hold
               options.Box = obj.Box
               options.Grid = obj.Grid
               options.PBaspect = obj.PBaspect
               options.DisplayName = ''
               options.Legend = 'on'
               options.LegendLabels = {}
               options.LegendOrientation = 'horizontal'
               options.LegendLocation = obj.LegendLocation
               options.LegendColumns = obj.LegendColumns
               options.LogScale = obj.LogScale
               options.YScale = obj.YScale
               options.XScale = obj.XScale
               options.Drawnow = obj.DrawNow
               options.LegendItemTokenSize = obj.LegendItemTokenSize
               options.AxesFontSize = obj.AxesFontSize
            end

          
           if ~isempty(XData) && ~isempty(YData)
                if options.LogScale
                   pl = loglog(XData, YData, options.LineSpec, 'LineWidth', options.LineWidth, 'MarkerSize', options.MarkerSize);
               else
                   pl = plot(XData, YData, options.LineSpec, 'LineWidth', options.LineWidth, 'MarkerSize', options.MarkerSize);
                   set(gca, 'XScale', options.XScale);
                   set(gca, 'YScale', options.YScale);
               end
           end
           % Set the diplayname
           DisplayNames = options.DisplayName;
           if ~iscell(DisplayNames)
               DisplayNames = {DisplayNames};
           end
           [s1, ~] = size(DisplayNames);
           if s1 == 1%Enforce a column cell array
                DisplayNames = DisplayNames';
           end
           if exist("pl", "var")
                set(pl, {'DisplayName'}, DisplayNames);
           end
           % Set the handle visibility
           HandleVisibilities = options.HandleVisibility;
           if ~iscell(HandleVisibilities)
               HandleVisibilities = {HandleVisibilities};
           end
           [s1, ~] = size(HandleVisibilities);
           if s1 == 1%Enforce a column cell array
                HandleVisibilities = HandleVisibilities';
           end
           if exist("pl", "var")
                set(pl, {'HandleVisibility'}, HandleVisibilities);
           end
           %
           if ~isempty(char(options.Color))
                pl.Color = options.Color;
           else
                colororder(gca(), obj.Colors);% If no color is specified use the color order
           end

           %Set the font of the axes
           ax = gca();
           ax.FontSize = options.AxesFontSize;
           ax.TickLabelInterpreter = options.Interpreter;

           %Set the x- and y-ticks
           if ~isempty(options.XTicks)
               xticks(options.XTicks);
           end
           if ~isempty(options.YTicks)
               yticks(options.YTicks);
           end

           %Set the x- and y-lim
           if ~isempty(options.XLim)
               xlim(options.XLim);
           end
           if ~isempty(options.YLim)
               ylim(options.YLim);
           end

           %Make the plot
           xlabel(options.XLabel, 'Interpreter', options.Interpreter, 'FontSize', options.FontSize);
           ylabel(options.YLabel, 'Interpreter', options.Interpreter, 'FontSize', options.FontSize);
           pbaspect(options.PBaspect);
           grid(options.Grid);
           hold(options.Hold);
           box(options.Box);

           
           if strcmp(options.Legend, 'on')
                lg = legend(options.LegendLabels, ...
                        'Interpreter', options.Interpreter, ...
                        'FontSize', options.FontSize, ...
                        'Orientation', options.LegendOrientation, ...
                        'Location', options.LegendLocation, ...
                        'NumColumns', options.LegendColumns);
                lg.ItemTokenSize = options.LegendItemTokenSize;
           end

           if options.Drawnow
               drawnow;
           end
        end

        %Function to save the plots
        function save(~, figs, path, varargin)
            path = char(path);
            
            pathcheck = @(s) s(end) == '/';

            p = inputParser;
            addRequired(p, 'figs', @(x) iscell(x));
            addRequired(p, 'path', pathcheck);
            addParameter(p, 'AppendName', "");
            addParameter(p, 'SaveAsMatlabFig', false);%Save also the image as a Matlab figure
            addParameter(p, "PreventAxisScaling", true);
            parse(p, figs, path, varargin{:});
           
            
            for i = 1:length(figs)
                FigAxes = gca(p.Results.figs{i});
                %Prevent any scaling of the axis
                if p.Results.PreventAxisScaling == true
                    FigAxes.XTickMode = 'manual';
                    FigAxes.YTickMode = 'manual';
                    FigAxes.ZTickMode = 'manual';
                end
                %Get the name of the figure
                FigName  = p.Results.figs{i}.FileName;
                if p.Results.AppendName == ""
                    FileName = string(p.Results.path) + string(FigName);
                else
                    FileName = string(p.Results.path) + string(p.Results.AppendName) + "_" + string(FigName);
                end

                [fPath, ~, ~]    = fileparts(FileName);
                if ~exist(fPath, "dir")
                    mkdir(fPath);
                end
                
                if contains(FigName, "pdf")
                    exportgraphics(FigAxes, FileName, ...
                    'BackgroundColor', 'none', ...
                    'ContentType', 'vector');
                else
                    exportgraphics(FigAxes, FileName, 'Resolution', 400);
                end

                % Check if the user asked to save as MATLAB fig
                if p.Results.SaveAsMatlabFig == true
                    % Remove extension from FileName
                    [pa, fi, ~]    = fileparts(FileName);
                    FullFilename = fullfile(pa, fi);
                    savefig(p.Results.figs{i}, FullFilename + ".fig");
                end
                
            end
        end

    end
end
