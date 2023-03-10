classdef uimaps_cbars_velocity < TComponent
    properties (Constant)
        Type = "axes"
    end
    properties (SetAccess = immutable)
        im
        tx_min; tx_max; tx_limit; tx_label; tx_highlight
        ln
        cm
    end

    methods
        function updateuimaps_highlight(obj)
            v = obj.Data.uimaps_highlight.speed;

            if isfinite(v)
                cmap = obj.Data.uimaps_settings.cmap_speed;

                n = size(cmap, 1);
                p = round(n/8);

                if v > obj.Data.uimaps_settings.max_speed
                    y = n + p/2;
                    set(obj.tx_highlight, ...
                        "Color", "r")
                    set(obj.ln, ...
                        "Color", "k")
                else
                    y = n * (v/obj.Data.uimaps_settings.max_speed);
                    set(obj.tx_highlight, ...
                        "Color", "k")
                    set(obj.ln, ...
                        "Color", "k")
                end

                set(obj.tx_highlight, ...
                    "String", string(round(v, 1)), ...
                    "Position", [1.2 y] + 0.5)
                set(obj.tx_label, ...
                    "Color", [224 224 224]/256)
                set(obj.ln, ...
                    "XData", [0 1] + 0.5, ...
                    "YData", [y y] + 0.5)
            else
                set(obj.tx_highlight, ...
                    "String", "")
                set(obj.tx_label, ...
                    "Color", [128 128 128]/256)
                set(obj.ln, ...
                    "XData", [], ...
                    "YData", [] + 0.5)
            end

        end
        function updateuimaps_settings(obj)
            cmap = obj.Data.uimaps_settings.cmap_speed;
            n = size(cmap, 1);
            p = round(n/8);

            set(obj.Handle, ...
                "YLim", [0 n + p] + 0.5, ...
                "XLim", [0 4] + 0.5, ...
                "Colormap", cmap, ...
                "Visible", "off")
            set(obj.im, ...
                "CData", ([1:n, ones(1, p) * n])')
            set(obj.tx_min, ...
                "String", string(0), ...
                "Position", [1.2 0] + 0.5)
            set(obj.tx_limit, ...
                "String", string(obj.Data.uimaps_settings.max_speed), ...
                "Position", [1.2 n] + 0.5)
            set(obj.tx_max, ...
                "String", "Inf", ...
                "Position", [1.2 n + p] + 0.5)
            set(obj.tx_label, ...
                "Position", [1.2 (n + p)/2] + 0.5)
        end
    end
    methods (Access = protected)
        function mouseAxesClickRight(obj)
            set(obj.cm, ...
                "Position", get(obj.Window, "CurrentPoint"), ...
                "Visible", "on")
        end
    end
    methods (Access = private)
        function menuFcn(obj, e)
            switch e.Source.Text
                case "5 mm/s"
                    obj.Data.uimaps_settings.setMaxSpeed(5)
                case "10 mm/s"
                    obj.Data.uimaps_settings.setMaxSpeed(10)
                case "20 mm/s"
                    obj.Data.uimaps_settings.setMaxSpeed(20)
                case "50 mm/s"
                    obj.Data.uimaps_settings.setMaxSpeed(50)
                otherwise 
                    return
            end
            set(e.Source.Parent.Children, ...
                "Checked", "off")
            set(e.Source, ...
                "Checked", "on")
        end
    end
    methods % CONSTRUCTOR
        function obj = uimaps_cbars_velocity()
            set(obj.Handle, ...
                ... Ticks
                'XTick', [], ...
                'YTick', [], ...
                ... Rulers
                'XLim', [0 4] + 0.5, ...
                'YLim', [0 1], ...
                'YDir', 'reverse', ...
                ... Colour
                'Color', 'w', ...
                ... Position
                'Units', 'normalized', ...
                'InnerPosition', [0.01 0.01 0.98 0.98], ...
                ... Interactivity
                'Visible', 'off', ...
                ... Callback Execution Control
                'PickableParts', 'all')

            obj.im = image(obj.Handle, ...
                ... Image Data and Quality
                "CData", [], ...
                "CDataMapping", "direct");
            obj.tx_min = text(obj.Handle, 0, 0, "", ...
                ... Text
                "Color", "k", ...
                ... Font
                "FontName", "fixedwidth", ...
                "FontWeight", "bold", ...
                "FontSize", 8, ...
                ... Text Box
                "Rotation", 0, ...
                "HorizontalAlignment", "left", ...
                "VerticalAlignment", "cap");
            obj.tx_limit = text(obj.Handle, 0, 0, "", ...
                ... Text
                "Color", "k", ...
                ... Font
                "FontName", "fixedwidth", ...
                "FontWeight", "bold", ...
                "FontSize", 8, ...
                ... Text Box
                "Rotation", 0, ...
                "HorizontalAlignment", "left", ...
                "VerticalAlignment", "middle");
            obj.tx_max = text(obj.Handle, 0, 0, "", ...
                ... Text
                "Color", "k", ...
                ... Font
                "FontName", "fixedwidth", ...
                "FontWeight", "bold", ...
                "FontSize", 8, ...
                ... Text Box
                "Rotation", 0, ...
                "HorizontalAlignment", "left", ...
                "VerticalAlignment", "baseline");
            obj.tx_label = text(obj.Handle, ...
                1.7, 0.5, "Velocity (mm/s)", ...
                ... Text
                "Color", [128 128 128]/256, ...
                ... Font
                "FontName", "fixedwidth", ...
                "FontWeight", "bold", ...
                "FontSize", 10, ...
                ... Text Box
                "Rotation", -90, ...
                "HorizontalAlignment", "center", ...
                "VerticalAlignment", "bottom");

            obj.ln = line(obj.Handle, ...
                ... Data
                'XData', double.empty(1,0), ...
                'YData', double.empty(1,0));
            obj.tx_highlight = text(obj.Handle, 0, 0, "", ...
                ... Text
                "Color", "k", ...
                ... Text Box
                "BackgroundColor", [obj.Handle.Parent.BackgroundColor 0.9], ...
                "Margin", 1, ...
                ... Font
                "FontName", "fixedwidth", ...
                "FontWeight", "bold", ...
                "FontSize", 8, ...
                ... Text Box
                "Rotation", 0, ...
                "HorizontalAlignment", "left", ...
                "VerticalAlignment", "middle");

            set(allchild(obj.Handle), "HitTest", "off")

            obj.cm = uicontextmenu(obj.Window);
            ti = uimenu(obj.cm, ...
                "Text", "Maximum velocity");
            uimenu(ti, ...
                "Text", "5 mm/s");
            uimenu(ti, ...
                "Text", "10 mm/s", ...
                "Checked", "on");
            uimenu(ti, ...
                "Text", "20 mm/s");
            uimenu(ti, ...
                "Text", "50 mm/s");

            obj.addlistener(findobj(obj.cm, "Children", gobjects(0)), ...
                "Action", @(~, e) obj.menuFcn(e))
        end
    end     % CONSTRUCTOR
end