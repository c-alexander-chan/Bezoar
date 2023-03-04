classdef uimenus < TComponent
    properties (Constant)
        Type = 'placeholder'
    end
    properties (SetAccess = immutable)
        mn_file_load
        mn_file_save
        mn_file_saveas
        mn_file_import

        mn_edit_detect
        mn_edit_clear

        mn_view
        mn_view_trace
        mn_view_movie
        mn_view_wave
        mn_view_maps

        mn_window_navigator
        mn_window_properties
    end
    
    methods
        function updatefile(obj)
            set(obj.mn_file_save, 'Enable', obj.Data.file.isloaded)
            set(obj.mn_file_saveas, 'Enable', obj.Data.file.isloaded)
            set(obj.mn_file_load, 'Enable', obj.Data.file.isloaded)
        end

        function updateuiview(obj)
            obj.update()
        end
        function update(obj)
            set(allchild(obj.mn_view), 'Checked', 'off')
            switch obj.Data.uiview.viewname
                case 'sign'
                    set(obj.mn_view_trace, 'Checked', 'on')
                case 'flmv'
                    set(obj.mn_view_movie, 'Checked', 'on')
                case 'info'
                    set(obj.mn_view_wave, 'Checked', 'on')
                case 'maps'
                    set(obj.mn_view_maps, 'Checked', 'on')
            end
            
            set(obj.mn_window_navigator, 'Checked', obj.Data.uiview.navigator)
            set(obj.mn_window_properties, 'Checked', obj.Data.uiview.property)
        end
    end
    methods % CONSTRUCTOR
        function obj = uimenus()
            d = obj.Data;
            w = obj.Window;

            m_file = uimenu(w, 'Text', '&File');
            obj.mn_file_save = uimenu(m_file, ...
                'Text', 'Save', ...
                'Accelerator', 'S', ...
                'Enable', 'off', ...
                'MenuSelectedFcn', @(~,~) d.save.save());
            obj.mn_file_saveas = uimenu(m_file, ...
                'Text', 'Save as...', ...
                'Enable', 'off', ...
                'MenuSelectedFcn', @(~,~) d.save.saveas());
            obj.mn_file_load = uimenu(m_file, ...
                'Text', 'Load', ...
                'Accelerator', 'L', ...
                'Enable', 'off', ...
                'MenuSelectedFcn', @(~,~) d.save.load());
            obj.mn_file_import = uimenu(m_file, ...
                'Text', 'Import file', ...
                'Separator', 'on', ...
                'Accelerator', 'I', ...
                'MenuSelectedFcn', @(~,~) d.file.import());

            m_edit = uimenu(w, 'Text', '&Edit');
            obj.mn_edit_detect = uimenu(m_edit, ...
                'Text', 'Detect waves', ...
                'Accelerator', 'D', ...
                'MenuSelectedFcn', @(~,~) d.evnt.detect());
            obj.mn_edit_clear = uimenu(m_edit, ...
                'Separator', 'on', ...
                'Text', 'Clear ALL markers...', ...
                'MenuSelectedFcn', @(~,~) d.evnt.clear());

            obj.mn_view = uimenu(w, 'Text', '&View');
            obj.mn_view_trace = uimenu(obj.mn_view, ...
                'Text', 'Signal', ...
                'Accelerator', '1', ...
                'MenuSelectedFcn', @(~,~) d.uiview.setView(1));
            obj.mn_view_movie = uimenu(obj.mn_view, ...
                'Text', 'Animation', ...
                'Accelerator', '2', ...
                'MenuSelectedFcn', @(~,~) d.uiview.setView(2));
            obj.mn_view_wave = uimenu(obj.mn_view, ...
                'Text', 'Slow wave', ...
                'Accelerator', '3', ...
                'MenuSelectedFcn', @(~,~) d.uiview.setView(3));
            obj.mn_view_maps = uimenu(obj.mn_view, ...
                'Text', 'Maps', ...
                'Accelerator', '4', ...
                'MenuSelectedFcn', @(~,~) d.uiview.setView(4));

            m_window = uimenu(w, 'Text', '&Window');
            obj.mn_window_navigator = uimenu(m_window, ...
                'Text', 'Navigator', ...
                'Accelerator', 'N', ...
                'MenuSelectedFcn', @(~,~) d.uiview.toggleVisibilityNavigator());
            obj.mn_window_properties = uimenu(m_window, ...
                'Text', 'Properties', ...
                'Accelerator', 'P', ...
                'MenuSelectedFcn', @(~,~) d.uiview.toggleVisibilityProperty());

            drawnow
            jf = get(handle(w), 'JavaFrame'); %#ok<JAVFM> 
            jmb = jf.fHG2Client.getMenuBar;
            jm = jmb.getComponent(0);
            jm.doClick;
            drawnow
            javax.swing.MenuSelectionManager.defaultManager().clearSelectedPath();
            drawnow
            jm = jm.getMenuComponent(1);
            ja = javax.swing.KeyStroke.getKeyStroke('shift ctrl S');
            jm.setAccelerator(ja);
        end
    end     % CONSTRUCTOR
end