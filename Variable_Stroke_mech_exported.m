classdef Variable_Stroke_mech_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        Image_2                       matlab.ui.control.Image
        Image                         matlab.ui.control.Image
        TYPE2Button                   matlab.ui.control.Button
        TYPE1Button                   matlab.ui.control.Button
        VARIABLESTROKEMECHANISMLabel  matlab.ui.control.Label
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: TYPE1Button
        function TYPE1ButtonPushed(app, event)
            Variable_Stroke_mechanism_t1_exported
        end

        % Button pushed function: TYPE2Button
        function TYPE2ButtonPushed(app, event)
            Variable_stroke_mechanism_t2_exported
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 837 542];
            app.UIFigure.Name = 'MATLAB App';

            % Create VARIABLESTROKEMECHANISMLabel
            app.VARIABLESTROKEMECHANISMLabel = uilabel(app.UIFigure);
            app.VARIABLESTROKEMECHANISMLabel.Position = [344 489 190 34];
            app.VARIABLESTROKEMECHANISMLabel.Text = 'VARIABLE STROKE MECHANISM';

            % Create TYPE1Button
            app.TYPE1Button = uibutton(app.UIFigure, 'push');
            app.TYPE1Button.ButtonPushedFcn = createCallbackFcn(app, @TYPE1ButtonPushed, true);
            app.TYPE1Button.Position = [142 30 135 32];
            app.TYPE1Button.Text = 'TYPE 1';

            % Create TYPE2Button
            app.TYPE2Button = uibutton(app.UIFigure, 'push');
            app.TYPE2Button.ButtonPushedFcn = createCallbackFcn(app, @TYPE2ButtonPushed, true);
            app.TYPE2Button.Position = [577 30 135 32];
            app.TYPE2Button.Text = 'TYPE 2';

            % Create Image
            app.Image = uiimage(app.UIFigure);
            app.Image.Position = [27 77 389 425];
            app.Image.ImageSource = 'TYPE-1.jpeg';

            % Create Image_2
            app.Image_2 = uiimage(app.UIFigure);
            app.Image_2.Position = [437 69 387 454];
            app.Image_2.ImageSource = 'TYPE-2.jpeg';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Variable_Stroke_mech_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end