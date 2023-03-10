classdef uiinfolabels < TComponent
    properties (Constant)
        Type = "axes"
    end

    methods % CONSTRUCTOR
        function obj = uiinfolabels()
            set(obj.Handle, ...
                ... Ticks
                'XTick', [], ...
                'YTick', [], ...
                ... Rulers
                'XLim', [0.5 1.5], ...
                'YLim', [0.5 1.5], ...
                'YDir', 'reverse', ...
                ... Colour
                'Color', 'w')
            text(obj.Handle, ...
                1, 1, 'Slow Wave Information', ...
                ... Text
                'Color', 'k', ...
                ... Font
                'FontSize', 10, ...
                'FontName', 'FixedWidth', ...
                'FontWeight', 'bold', ...
                ... Text Box
                'Rotation', 90, ...
                ... Position
                'HorizontalAlignment', 'center', ...
                ... Callbacks
                'PickableParts', 'none')
        end
    end     % CONSTRUCTOR
    methods (Access = protected)
        function initialise(obj)
            obj.Width = 20;
        end
    end     % INITIALISE
end