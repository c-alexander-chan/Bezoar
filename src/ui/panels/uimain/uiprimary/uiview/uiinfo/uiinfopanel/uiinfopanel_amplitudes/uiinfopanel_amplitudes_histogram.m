classdef uiinfopanel_amplitudes_histogram < TComponent
    properties (Constant)
        Type = "axes"
    end
    properties (SetAccess = immutable)
        cl_cursor; rc_cursor
        hs_all; hs_selected
        tx_min; tx_max
    end
    properties (Access = private)
        amp_range           (1,:)   double  = double.empty(1,0)
        amp_range_last      (1,:)   double  = double.empty(1,0)
        selection_last              logical = []
    end

    methods
        function updateuiprvw(obj)
            obj.amp_range = [];
            obj.amp_range_last = [];
        end
        function updateuiinfo_amplitudes(obj)
            if ~isequal(obj.selection_last, obj.Data.uiinfo_selection.selection)
                obj.amp_range = [];
                obj.amp_range_last = [];
                obj.selection_last = [];
            end

            y = [-0.25 1.10];
            if any(obj.Data.uiinfo_amplitudes.amplitudes, "all")
                set(obj.Handle, ...
                    "XLim", [1 obj.Data.uiinfo_amplitudes.bins], ...
                    "YLim", y * max(obj.Data.uiinfo_amplitudes.counts))
            else
                set(obj.Handle, ...
                    "XLim", [1 2], ...
                    "YLim", y)
            end
            set(obj.tx_min, ...
                "String", obj.Data.uiinfo_amplitudes.str_min, ...
                "Position", [1 0])
            set(obj.tx_max, ...
                "String", obj.Data.uiinfo_amplitudes.str_max, ...
                "Position", [obj.Data.uiinfo_amplitudes.bins 0])
            set(obj.hs_all, ...
                "BinEdges", 1:obj.Data.uiinfo_amplitudes.bins, ...
                "BinCounts", obj.Data.uiinfo_amplitudes.counts);
            set(obj.hs_selected, ...
                "BinEdges", 1:obj.Data.uiinfo_amplitudes.bins, ...
                "BinCounts", obj.Data.uiinfo_amplitudes.counts_selected);

            if (length(obj.amp_range) == 2) || ...
                    any(obj.Data.uiinfo_selection.selection, "all")
                set(obj.hs_all, ...
                    "FaceAlpha", 0.3)
            else
                set(obj.hs_all, ...
                    "FaceAlpha", 0.7)
            end

            obj.updateCursor()
        end
    end
    methods (Access = protected)
        function traverseFcn(obj)
            switch length(obj.amp_range)
                case {0, 1}                   
                    obj.amp_range = obj.CurrentPosition(1);
            end
            obj.updateCursor()
        end
        function exitFcn(obj)
            switch length(obj.amp_range)
                case 1
                    obj.amp_range = [];
            end
            obj.updateCursor()
        end

        function mouseDrag(obj)
            p0 = obj.ClickedPosition(1);
            p1 = obj.CurrentPosition(1);
            
            e = obj.Data.uiinfo_amplitudes.edges;
            a = obj.Data.uiinfo_amplitudes.amplitudes;

            r = sort([p0 p1]);
            r = [floor(r(1)) ceil(r(2))];
            r = max(min(r, length(e)), 1);

            s = (a > e(r(1))) & (a < e(r(2)));

            obj.amp_range = r;
            obj.selection_last = s;

            obj.Data.uiinfo_selection.setSelection(s)
        end
        function mouseRelease(obj)
            p0 = obj.ClickedPosition(1);
            p1 = obj.CurrentPosition(1);
            
            e = obj.Data.uiinfo_amplitudes.edges;

            r = sort([p0 p1]);
            r = [floor(r(1)) ceil(r(2))];
            r = max(min(r, length(e)), 1);

            switch length(obj.amp_range_last)
                case 0
                case 1
                case 2
                    if (range(r) == 1) && isequal(r, obj.amp_range)
                        obj.amp_range = p1;
                        obj.Data.uiinfo_selection.setSelection([])
                    end
            end

            obj.amp_range_last = obj.amp_range;
        end
    end
    methods (Access = private)
        function updateAxes(obj)
            y = [-0.25 1.10];
            set(obj.Handle, ...
                "XLim", [1 obj.Data.uiinfo_amps_.bins], ...
                "YLim", y*max(obj.Data.uiinfo_amps_.counts))
        end
        function updateCursor(obj)
            e = obj.Data.uiinfo_amplitudes.edges;
            n = obj.Data.uiinfo_amplitudes.bins;

            s = obj.amp_range;
            switch length(s)
                case 0
                    p = 0;
                    set(obj.cl_cursor, ...
                        "Visible", "off")
                    set(obj.rc_cursor, ...
                        "Visible", "off")
                case 1
                    p = s(1);
                    set(obj.cl_cursor, ...
                        "Value", p, ...
                        "Label", num2str(e(max(min(round(s), n), 1)), '%.1f'), ...
                        "Visible", "on")
                    set(obj.rc_cursor, ...
                        "Visible", "off")
                case 2
                    p = s(2);
                    set(obj.cl_cursor, ...
                        "Value", p, ...
                        "Label", num2str(e(s(1)), '%.1f') + " to " + num2str(e(s(2)), '%.1f'), ...
                        "Visible", "on")
                    set(obj.rc_cursor, ...
                        "Position", [s(1) -256 range(s) 512], ...
                        "Visible", "on")
            end

            if p > (obj.Data.uiinfo_amplitudes.bins - 1)/2
                set(obj.cl_cursor, ...
                    "LabelHorizontalAlignment", "left")
            else
                set(obj.cl_cursor, ...
                    "LabelHorizontalAlignment", "right")
            end
        end
        function updateLabels(obj)
            set(obj.tx_min, ...
                "String", obj.Data.uiinfo_amplitudes.str_min, ...
                "Position", [1 0])
            set(obj.tx_max, ...
                "String", obj.Data.uiinfo_amplitudes.str_max, ...
                "Position", [obj.Data.uiinfo_amplitudes.bins 0])
        end
        function updateHistograms(obj)
            set(obj.hs_all, ...
                "BinEdges", 1:obj.Data.uiinfo_amplitudes.bins, ...
                "BinCounts", obj.Data.uiinfo_amplitudes.counts);
            set(obj.hs_selected, ...
                "BinEdges", 1:obj.Data.uiinfo_amplitudes.bins, ...
                "BinCounts", obj.Data.uiinfo_amplitudes.counts_selected);
        end
    end
    methods % CONSTRUCTOR
        function obj = uiinfopanel_amplitudes_histogram()
            set(obj.Handle, ...
                ... Ticks
                'TickLength', [0 0], ...
                ... Rulers
                'XLim', [0 1], ...
                'YLim', [0 1], ...
                'XLimMode', 'manual', ...
                'YLimMode', 'manual', ...
                'XColor', 'none', ...
                'YColor', 'none', ...
                ... Multiple Plots
                'NextPlot', 'add', ...
                ... Box Styling
                'Color', [1 1 1], ...
                ... Position
                'Position', [0 0 1 1], ...
                'PositionConstraint', 'innerposition', ...
                ... Interactivity
                'Visible', 'off', ...
                ... Callback Execution Control
                'PickableParts', 'all');
            yline(obj.Handle, 0, ...
                ... Color and Styling
                'Color', [0.8 0.8 0.8])
            obj.cl_cursor = xline(obj.Handle, 0, ...
                ... Labels
                'LabelOrientation', 'horizontal', ...
                ... Color and Styling
                'Color', 'w', ...
                ... Font
                'FontName', 'fixedwidth', ...
                'FontSize', 8, ...
                'FontWeight', 'bold');

            obj.hs_all = histogram(obj.Handle, [],  ...
                ... Data
                'Normalization', 'count', ...
                ... Color and Styling
                'FaceColor', 'w', ...
                'FaceAlpha', 0.3, ...
                'LineStyle', 'none');
            obj.hs_selected = histogram(obj.Handle, [],  ...
                ... Data
                'Normalization', 'count', ...
                ... Color and Styling
                'FaceColor', 'w', ...
                'FaceAlpha', 1, ...
                'LineStyle', 'none');
            obj.rc_cursor = rectangle(obj.Handle, ...
                ... Color and Styling
                'FaceColor', [256 256 256 032]/256, ...
                'EdgeColor', 'none', ..., ...
                'AlignVertexCenters', 'on', ...
                ... Position
                'Position', [0 0 1 1]);
            obj.tx_min = text(obj.Handle, 0, 0, '', ...
                ... Color and Styling
                'Color', 'w', ...
                ... Font
                'FontName', 'FixedWidth', ...
                'FontSize', 8, ...
                'FontWeight', 'bold', ...
                ... Position
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top');
            obj.tx_max = text(obj.Handle, 0, 0, '', ...
                ... Color and Styling
                'Color', 'w', ...
                ... Font
                'FontName', 'FixedWidth', ...
                'FontSize', 8, ...
                'FontWeight', 'bold', ...
                ... Position
                'HorizontalAlignment', 'right', ...
                'VerticalAlignment', 'top');

            set(allchild(obj.Handle), ...
                'PickableParts', 'none')
        end
    end     % CONSTRUCTOR
end