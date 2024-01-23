classdef Plotter
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
    end
    
    methods
        function obj = Plotter()
        end
        
        %Custom plot function
        function plot(obj, varargin)
           % Parse the input
           p = inputParser;
           addOptional(p,  'XData', []);
           addOptional(p,  'YData', []);
           addParameter(p, 'LineSpec', '-');
           addParameter(p, 'LineWidth', obj.LineWidth);
           addParameter(p, 'Color', '');
           addParameter(p, 'MarkerSize', obj.MarkerSize);
           addParameter(p, 'HandleVisibility', obj.HandleVisibility);
           addParameter(p, 'FontSize', obj.FontSize);
           addParameter(p, 'Interpreter', obj.Interpreter);
           addParameter(p, 'XLabel', obj.XLabel);
           addParameter(p, 'YLabel', obj.YLabel);
           addParameter(p, 'XTicks', []);
           addParameter(p, 'YTicks', []);
           addParameter(p, 'XLim', []);
           addParameter(p, 'YLim', []);
           addParameter(p, 'Hold', obj.Hold);
           addParameter(p, 'Box', obj.Box);
           addParameter(p, 'Grid', obj.Grid);
           addParameter(p, 'PBaspect', obj.PBaspect);
           addParameter(p, 'DisplayName', '');
           addParameter(p, 'Legend', 'on');
           addParameter(p, 'LegendLabels', {});
           addParameter(p, 'LegendOrientation', 'horizontal');
           addParameter(p, 'LegendLocation', obj.LegendLocation);
           addParameter(p, 'LegendColumns', obj.LegendColumns);
           addParameter(p, 'LogScale', obj.LogScale);
           addParameter(p, 'YScale', obj.YScale);
           addParameter(p, 'XScale', obj.XScale);
           addParameter(p, 'Drawnow', obj.DrawNow);
           addParameter(p, 'LegendItemTokenSize', obj.LegendItemTokenSize);
           addParameter(p, 'AxesFontSize', obj.AxesFontSize);

           parse(p, varargin{:});
           
           X_data = p.Results.XData;
           Y_data = p.Results.YData;

           if ~isempty(X_data) && ~isempty(Y_data)
               if p.Results.LogScale
                   pl = loglog(X_data, Y_data, p.Results.LineSpec, 'LineWidth', p.Results.LineWidth, 'MarkerSize', p.Results.MarkerSize, 'HandleVisibility', p.Results.HandleVisibility, 'DisplayName', p.Results.DisplayName);
               else
                   pl = plot(X_data, Y_data, p.Results.LineSpec, 'LineWidth', p.Results.LineWidth, 'MarkerSize', p.Results.MarkerSize, 'HandleVisibility', p.Results.HandleVisibility, 'DisplayName', p.Results.DisplayName);
                   set(gca, 'XScale', p.Results.XScale);
                   set(gca, 'YScale', p.Results.YScale);
               end
           end
           %
           if ~isempty(char(p.Results.Color))
                pl.Color = p.Results.Color;
           end

           %Set the font of the axes
           ax = gca();
           ax.FontSize = p.Results.AxesFontSize;

           %Set the x- and y-ticks
           if ~isempty(p.Results.XTicks)
               xticks(p.Results.XTicks);
           end
           if ~isempty(p.Results.YTicks)
               yticks(p.Results.YTicks);
           end

           %Set the x- and y-lim
           if ~isempty(p.Results.XLim)
               xlim(p.Results.XLim);
           end
           if ~isempty(p.Results.YLim)
               ylim(p.Results.YLim);
           end

           %Make the plot
           xlabel(p.Results.XLabel, 'Interpreter', p.Results.Interpreter, 'FontSize', p.Results.FontSize);
           ylabel(p.Results.YLabel, 'Interpreter', p.Results.Interpreter, 'FontSize', p.Results.FontSize);
           pbaspect(p.Results.PBaspect);
           grid(p.Results.Grid);
           hold(p.Results.Hold);
           box(p.Results.Box);

           
           if strcmp(p.Results.Legend, 'on')
                lg = legend(p.Results.LegendLabels, ...
                        'Interpreter', p.Results.Interpreter, ...
                        'FontSize', p.Results.FontSize, ...
                        'Orientation', p.Results.LegendOrientation, ...
                        'Location', p.Results.LegendLocation, ...
                        'NumColumns', p.Results.LegendColumns);
                lg.ItemTokenSize = p.Results.LegendItemTokenSize;
           end

           if p.Results.Drawnow
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
            parse(p, figs, path, varargin{:});
           
            
            for i = 1:length(figs)
                FigAxes = gca(p.Results.figs{i});
                %Prevent any scaling of the axis
                FigAxes.XTickMode = 'manual';
                FigAxes.YTickMode = 'manual';
                FigAxes.ZTickMode = 'manual';
                %Get the name of the figure
                FigName  = p.Results.figs{i}.FileName;
                if p.Results.AppendName == ""
                    FileName = string(p.Results.path) + string(FigName);
                else
                    FileName = string(p.Results.path) + string(p.Results.AppendName) + "_" + string(FigName);
                end
                
                if contains(FigName, "pdf")
                    exportgraphics(FigAxes, FileName, ...
                    'BackgroundColor', 'none', ...
                    'ContentType', 'vector');
                else
                    exportgraphics(FigAxes, FileName, 'Resolution', 400);
                end
                
            end
        end

    end
end
