%% Basic Image viewer rewritten by Mira M. Liu November 2021 (from July 2019)
%  current functions: include 'permute' as input if the stack is of
% the form [slice, x,y]. if it is of the form [x,y,slice] 'permute' is not
% needed. 
% can include name as second input to put a name of the stack of images
% sets to middle slice, able to scroll or click through. select image to be able to use arrow keys to move through. 
% able to change the minimum and maximum pixel value
% figure scales with window size
% can zoom in and out 

% any questions direct to mira liu at liusarkarm@uchicago.edu
% https://miramliu.com/

classdef imagestack < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure           matlab.ui.Figure
        %MaxEditFieldLabel  matlab.ui.control.Label
        MaxColor           matlab.ui.control.NumericEditField
        %MinEditFieldLabel  matlab.ui.control.Label
        MinColor           matlab.ui.control.NumericEditField
        SliceSliderLabel   matlab.ui.control.Label
        Slice              matlab.ui.control.Slider
        UIAxes             matlab.ui.control.UIAxes
        Slider            matlab.ui.control.Slider %added to replace the textbok and make it a slider
        InputImage % Input stack of images
        Ss
        dcm
        MaxValue
        MinValue
        color
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, varargin)
            if isempty(varargin)
                error('No input images')
            elseif length(size(varargin{1})) ==3 %if it's three dimensional
                app.InputImage = varargin{1};
                name = 'Figure';
                app.color = 'jet';
                if size(varargin,2) ==2
                    if strcmp(string(varargin{2}),'permute')
                        app.InputImage = permute(varargin{1},[2 3 1]); %if it is given as [slice,x,y]
                    elseif strcmp(string(varargin{2}),'gray')
                        app.color = 'gray';
                    else
                        name = varargin{2};
                    end
                end
            elseif length(size(varargin{1}))==2 %if it's 2 dimensional, just gotta add one fake dimension to be compatible with other code...
                nx = size(varargin{1},2); %the dimensions of the image
                B = reshape(varargin{1},[1,nx,nx]); %assuming it's a square image...
                app.InputImage=permute(B,[2 3 1]);
            else
                error('Wrong Image Input: data must be 3d image [x,y,slice] or input [slice,x,y], "permute"')
            end
            

            % make datacursor usable
            %datacursormode(app.UIFigure, 'on')
            %app.dcm = datacursormode(app.UIFigure);
            %app.color = 'gray';
            %app.color = 'jet';

            app.MaxValue = max(app.InputImage,[],'all'); %set range across ALL slices
            app.MinValue = 0.0; %set to 0

            colormapname = app.color; %set colormap
            data = app.InputImage(:,:,round(app.Ss/2));
            app.UIAxes.reset
            imshow(data,[],'parent',app.UIAxes)
            title(name,'parent',app.UIAxes)
            colormap(app.UIAxes,colormapname),colorbar(app.UIAxes);
        end

        % Value changed function: MaxColor
        function SliderColorValueChanged(app,event)
            caxis(app.UIAxes, [double(0) double(app.Slider.Value)])
        end

        % Value changed function: Slice
        function SliceValueChanged(app, event)
            colormapname = app.color; %set colormap
            data = app.InputImage(:,:,round(double(app.Slice.Value)));
            imshow(data,[double(app.MinValue) double(app.Slider.Value)],'parent',app.UIAxes)
            colormap(app.UIAxes,colormapname),colorbar(app.UIAxes);
            app.SliceSliderLabel.Text = string('Slice ' + string(round(double(app.Slice.Value))));
            %app.MaxValue = max(app.InputImage(:,:,round(app.Slice.Value)),[],'all'); %update range
        end

        function SliceValueChangedKey(app,event)
            colormapname = app.color; %set colormap
            %colormapname = 'gray'; %set colormap
            if strcmp(event.Key, 'leftarrow')
                if round(app.Slice.Value)~= 1
                    app.Slice.Value = round(app.Slice.Value) - 1;
                    app.SliceSliderLabel.Text = ['Slice' ' ' num2str(round(app.Slice.Value))];
                    data = app.InputImage(:,:,round(double(app.Slice.Value)));
                    imshow(data,[double(app.MinValue) double(app.Slider.Value)],'parent',app.UIAxes)
                    colormap(app.UIAxes,colormapname),colorbar(app.UIAxes);
                    %app.MaxValue = max(app.InputImage(:,:,round(app.Slice.Value)),[],'all'); %update range
                end
            end
            if strcmp(event.Key, 'rightarrow')
                if round(app.Slice.Value)~= app.Ss
                    app.Slice.Value = round(app.Slice.Value) + 1;
                    app.SliceSliderLabel.Text = ['Slice' ' ' num2str(round(app.Slice.Value))];
                    data = app.InputImage(:,:,round(double(app.Slice.Value)));
                    imshow(data,[double(app.MinValue) double(app.Slider.Value)],'parent',app.UIAxes)
                    colormap(app.UIAxes,colormapname),colorbar(app.UIAxes);
                    %app.MaxValue = max(app.InputImage(:,:,round(app.Slice.Value)),[],'all'); %update range
                end
            end
            if strcmp(event.Key,'v') % vertical flip (over horizontal axis)
                app.InputImage = flip(app.InputImage,1);
                data = app.InputImage(:,:,round(double(app.Slice.Value)));
                imshow(data,[double(app.MinValue) double(app.Slider.Value)],'parent',app.UIAxes)
                colormap(app.UIAxes,colormapname),colorbar(app.UIAxes);
                %app.MaxValue = max(app.InputImage(:,:,round(app.Slice.Value)),[],'all'); %update range
            end
            if strcmp(event.Key, 'h')
                app.InputImage = flip(app.InputImage,2); % horizontal flip (over vertical axis)
                data = app.InputImage(:,:,round(double(app.Slice.Value)));
                imshow(data,[double(app.MinValue) double(app.Slider.Value)],'parent',app.UIAxes)
                colormap(app.UIAxes,colormapname),colorbar(app.UIAxes);
                %app.MaxValue = max(app.InputImage(:,:,round(app.Slice.Value)),[],'all'); %update range
            end
            
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app,varargin)
            app.Ss = size(varargin{1},3); % number of slices (assuming x,y,slice)
            if size(varargin,2) ==2 %checking order of matrix
                if strcmp(string(varargin{2}),'permute')
                	app.Ss = size(varargin{1},1); %if it is given as [slice,x,y]
                end
            end
            if length(size(varargin{1})) ==3 %if it's three dimensional
                app.InputImage = varargin{1};
                name = 'Figure';
                if size(varargin,2) ==2
                    if strcmp(string(varargin{2}),'permute')
                        app.InputImage = permute(varargin{1},[2 3 1]); %if it is given as [slice,x,y]
                    else
                        %name = varargin{2}(82:end-4);
                        %name=strrep(name,'_','-');
                        name = varargin{2};
                    end
                end
            elseif length(size(varargin{1}))==2 %if it's 2 dimensional, just gotta add one fake dimension to be compatible with other code...
                nx = size(varargin{1},2); %the dimensions of the image
                B = reshape(varargin{1},[1,nx,nx]); %assuming it's a square image...
                app.InputImage=permute(B,[2 3 1]);
            end

            app.MaxValue = real(max(app.InputImage,[],'all'));
            app.MinValue = max(real(min(app.InputImage,[],'all')),0); %if min is < 0, set min to 0?

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 652 527];
            app.UIFigure.Name = 'imagestack';
            app.UIFigure.KeyPressFcn = createCallbackFcn(app, @SliceValueChangedKey, true);
            %app.UIFigure.KeyPressFcn = createCallbackFcn(app, @Flip, true);

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            app.UIAxes.Position = [68 49 577 471];
            % Show the figure after all components are created
            %app.UIFigure.Visible = 'on';

            % Create Slider
            app.Slider = uislider(app.UIFigure);
            app.Slider.Orientation = 'vertical';
            app.Slider.Position = [18 57 3 446];
            app.Slider.Limits = [app.MinValue app.MaxValue]; %the slider limits set here, arbitrarily to 1/4 of the maximum value just because of noise spikes dominating
            app.Slider.Value= app.MaxValue/2;%250; %just standard max for perfusion MRI
            app.Slider.ValueChangedFcn = createCallbackFcn(app, @SliderColorValueChanged, true);

            % Create Slice
            if app.Ss>1
                app.Slice = uislider(app.UIFigure);
                app.Slice.Limits = [1 app.Ss];
                app.Slice.MinorTicks = [];
                app.Slice.ValueChangedFcn = createCallbackFcn(app, @SliceValueChanged, true);
                if app.Ss> 5
                    app.Slice.MajorTicks = [1 round(app.Ss/4) round(app.Ss/2) round(app.Ss*3/4) app.Ss];
                    app.Slice.Position = [210 37 279 3];
                    app.Slice.Value = round(app.Ss/2);
                else
                    app.Slice.MajorTicks = [1 2 3 4 5];
                    app.Slice.Position = [210 37 279 3];
                    app.Slice.Value = round(app.Ss/2);
                end
            end
                

            % Create SliceSliderLabel
            
            if app.Ss>1
                app.SliceSliderLabel = uilabel(app.UIFigure);
                app.SliceSliderLabel.HorizontalAlignment = 'right';
                app.SliceSliderLabel.Position = [69 18 95 22];
                app.SliceSliderLabel.Text = ['Slice' ' ' num2str(round(app.Slice.Value))];
            end
            
            
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = imagestack(varargin)

            % Create UIFigure and components
            createComponents(app,varargin{:})

            % Execute the startup function
            runStartupFcn(app, @(app)startupFcn(app, varargin{:}))
            app.UIFigure.Visible = 'on';

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
