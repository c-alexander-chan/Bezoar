classdef uipreviewinfo < TComponent
    properties (Constant)
        Type = "panel"
    end
    properties (SetAccess = immutable)
        uc
    end

    methods
        function updateuiprvw(obj)
            w = obj.Data.uiprvw.at;
            n = obj.Data.uiprvw.i;

            w_min = min(w, [], "all");
            w_max = max(w, [], "all");

            if n == 0
                strWave = "-";
            else
                strWave = string(n);
            end

            if all(isnan(w))
                str = "";
            else
                str = ...
                    strWave + ...
                    newline + ...
                    sprintf("%d", sum(~isnan(w(:)))) + " / " + num2str(numel(w)) + ...
                    newline + ...
                    sprintf("%.2f", w_max - w_min) + " seconds" + ...
                    newline + ...
                    "" + char(duration(0, 0, w_min), "mm:ss") + ...
                    "   " + sprintf("%.1f", w_min) + " s" + ...
                    newline + ...
                    "" + char(duration(0, 0, w_max), "mm:ss") + ...
                    "   " + sprintf("%.1f", w_max) + " s";
            end

            set(obj.uc, "String", str)
        end
    end
    methods % CONSTRUCTOR
        function obj = uipreviewinfo()
            obj.Height = 80;

            set(obj.Handle, ...
                "BackgroundColor", [048 048 048]/256, ...
                "Padding", 3)

            h = uix.HBox( ...
                "Parent", obj.Handle, ...
                "BackgroundColor", [048 048 048]/256, ...
                "Padding", 0, ...
                "Spacing", 7);

            uicontrol(h, ...
                ... Type of Control
                "Style", "Text", ...
                ... Text and Styling
                "String", propertyNames(), ...
                "BackgroundColor", [048 048 048]/256, ...
                "ForegroundColor", [240 240 240]/256, ...
                ... Font
                "FontName", "FixedWidth", ...
                "FontSize", 8, ...
                ... Position
                "HorizontalAlignment", "right")
            obj.uc = uicontrol(h, ...
                ... Type of Control
                "Style", "Text", ...
                ... Text and Styling
                "String", "", ...
                "BackgroundColor", [048 048 048]/256, ...
                "ForegroundColor", [240 240 240]/256, ...
                ... Font
                "FontName", "FixedWidth", ...
                "FontSize", 8, ...
                ... Position
                "HorizontalAlignment", "left");

            set(h, "Widths", [120 -1])
        end
    end     % CONSTRUCTOR
end

function str = propertyNames()
str = ...
    "   Wave number:" + newline + ...
    "Valid channels:" + newline + ...
    "    Time range:" + newline + ...
    "     Min. time:" + newline + ...
    "     Max. time:" + newline;
end