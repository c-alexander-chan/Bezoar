classdef uimaps_velocities_axes_sphere < TComponent
    properties (Constant)
        Type = 'axes'
    end
    properties (SetAccess = immutable)
        index
    end
    properties (SetAccess = immutable)
        sf_bg
        sf
        qv
    end
    properties (Access = private)
        clicked
    end

    methods
        function updatecnfg(obj)
            switch obj.Data.cnfg.mode
                case 'flat'
                    set(obj.sf_bg, ...
                        'Visible', 'off')
                    set(obj.sf, ...
                        'Visible', 'off')
                    set(obj.qv, ...
                        'Visible', 'off')
                    set(obj.Handle, ...
                        'Visible', 'off', ...
                        'PickableParts', 'none')
                case 'sphere'
                    [x, y, z] = sph2cart( ...
                        obj.Data.cnfg.xg, ...
                        obj.Data.cnfg.yg, ...
                        1);
                    set(obj.sf_bg, ...
                        'XData', x * 0.98, ...
                        'YData', y * 0.98, ...
                        'ZData', z * 0.98, ...
                        'EdgeColor', [128 128 128]/256, ...
                        'Visible', 'on')
                    set(obj.sf, ...
                        'XData', x, ...
                        'YData', y, ...
                        'ZData', z, ...
                        'EdgeColor', 'none', ...
                        'Visible', 'on')
                    set(obj.qv, ...
                        'Visible', 'on')
                    set(obj.Handle, ...
                        'Visible', 'off', ...
                        'PickableParts', 'all')
            end
        end
        function updateuimaps(obj)
            if strcmp('sphere', obj.Data.cnfg.mode) && ...
                    (obj.index <= obj.Data.uimaps.n)

                ig = obj.Data.uimaps.selection(obj.index);
                at = obj.Data.wave.waves(:, :, ig);

                [v, b] = obj.Data.a_velo.calculateVelocities(at);

                set(obj.sf, ...
                    'EdgeColor', 'none')
                set(obj.sf, ...
                    'CData', padarray(v, [1 1], NaN, 'post'))

                xv = obj.Data.a_velo.xv;
                yv = obj.Data.a_velo.yv;

                nv = numel(v);
                [xq, yq, zq, uq, vq, wq] = deal(NaN(1, nv));
                for i = 1:nv
                    azi = xv(i);
                    ele = yv(i);

                    N = rotz(pi - azi) * ...
                        roty(pi/2 - ele) * ...
                        [ ...
                            1   0   0; ...
                            0  -1   0; ...
                            0   0   1];
                    vec = N * ...
                        [...
                        cos(b(i) + pi); ...
                        sin(b(i) + pi); ...
                        0];

                    [xq(i), yq(i), zq(i)] = sph2cart(azi, ele, 1);
                    [uq(i), vq(i), wq(i)] = deal(vec(1), vec(2), vec(3));
                end

                set(obj.qv, ...
                    'UData', uq, ...
                    'VData', vq, ...
                    'WData', wq, ...
                    'XData', xq, ...
                    'YData', yq, ...
                    'ZData', zq)
                set(obj.Handle, ...
                    'CLim', [0 10], ...
                    'Colormap', cool, ...
                    'Visible', 'off')
            end
        end
        function updateuimaps_view(obj)
            set(obj.Handle, ...
                'View', obj.Data.uimaps_view.view)
        end
    end
    methods (Access = protected)
        function initialise(obj)
            
        end
        function mouseClickLeft(obj)
            obj.clicked = get(0, 'PointerLocation');
        end
        function mouseDrag(obj)
            p = get(0, 'PointerLocation');
            v = obj.Data.uimaps_view.view;

            obj.Data.uimaps_view.setView(v + obj.clicked - p)

            obj.clicked = p;
        end
    end
    methods % CONSTRUCTOR
        function obj = uimaps_velocities_axes_sphere()
            set(obj.Handle, ...
                ... Rulers
                'XLim', [-1 1], ...
                'YLim', [-1 1], ...
                'ZLim', [-1 1], ...
                'XColor', [192 192 192]/256, ...
                'YColor', [192 192 192]/256, ...
                'YDir', 'reverse', ...
                ... Box Styling
                'Color', 'w', ...
                'Box', 'on', ...
                ... Position
                'Units', 'normalized', ...
                'InnerPosition', [0.01 0.01 0.98 0.98], ...
                'DataAspectRatio', [1 1 1], ...
                ... View
                'View', [60 30], ...
                'CameraViewAngle', 10, ...
                ... Callback Execution Control
                'PickableParts', 'all')

            obj.sf_bg = surf([], ...
                'Parent', obj.Handle, ...
                ... Faces
                'FaceColor', [224 224 224]/256, ...
                'FaceAlpha', 0.8, ...
                ... Edges
                'EdgeColor','none', ...
                ... Callback Execution Control
                'PickableParts', 'none');
            obj.sf = surf([], 'Parent', obj.Handle, ...
                ... Faces
                'FaceColor', 'flat', ...
                ... Edges
                'EdgeColor', 'none', ...
                ... Callback Execution Control
                'PickableParts', 'none');
            obj.qv = quiver3(obj.Handle, ...
                [], [], [], [], ...
                ... Arrows
                'Color', [032 032 032]/256, ...
                'AutoScaleFactor', 0.6, ...
                ... Callback Execution Control
                'PickableParts', 'none');

            obj.index = obj.Parent.index;
        end
    end     % CONSTRUCTOR
end

function [Rx] = rotx(a)
Rx = [ ...
    1       0       0       ;
    0       cos(a)  sin(a)  ;
    0       -sin(a) cos(a)  ;];
end
function [Rz] = rotz(a)
Rz = [ ...
    cos(a)  sin(a)  0       ;
    -sin(a) cos(a)  0       ;
    0       0       1       ;];
end
function [Ry] = roty(a)
Ry = [ ...
    cos(a)  0       -sin(a) ;
    0       1       0       ;
    sin(a)  0       cos(a)  ;];
end