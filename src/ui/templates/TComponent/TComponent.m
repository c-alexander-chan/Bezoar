classdef TComponent < ...
        TComponentHeader & ...
        TCAxes & ...
        TCClick & ...
        TCClickAxes & ...
        TCKeyPress & ...
        TCModifiers & ...
        TCMouseOver & ...
        TCMouseScroll & ...
        TCResize & ...
        matlab.mixin.Heterogeneous

    methods
        function obj = TComponent()
            %disp(class(obj))
        end
    end
end
