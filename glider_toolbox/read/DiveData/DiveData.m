% NOTE - to update version text displayed in main window, search for VERSION_NUMBER

% ********************************************************************************
% DiveData M-file for DiveData.fig
%      DiveData, by itself, creates a new DiveData or raises the existing
%      singleton*.
%
%      H = DiveData returns the handle to a new DiveData or the handle to
%      the existing singleton*.
%
%      DiveData('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DiveData.M with the given input arguments.
%
%      DiveData('Property','Value',...) creates a new DiveData or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DiveData_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DiveData_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% ********************************************************************************
function varargout = DiveData(varargin)

    % Edit the above text to modify the response to help DiveData

    % Last Modified by GUIDE v2.5 14-Jun-2014 17:58:41

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @DiveData_OpeningFcn, ...
                       'gui_OutputFcn',  @DiveData_OutputFcn, ...
                       'gui_LayoutFcn',  [] , ...
                       'gui_Callback',   []);
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end

%     try
        if nargout
            [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
        else
            gui_mainfcn(gui_State, varargin{:});
        end
%     catch exception
%         for i = 1 : length(exception.stack)
%             disp(exception.message);
%             disp(exception.stack(i).file);
%             disp(exception.stack(i).name);
%             disp(exception.stack(i).line);
%         end
%     end
    % End initialization code - DO NOT EDIT
end

% ********************************************************************************
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DiveData (see VARARGIN)
%
% J.E. Stranzl Jr.
% 29-Jul-2012
% Updated version number to reflect incorporation of speed-of-sound changes
% 
% ********************************************************************************
function DiveData_OpeningFcn(hObject, eventdata, handles, varargin)

    % Choose default command line output for DiveData
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);

    initialize_gui(hObject, handles, false);

    set(handles.figure1, 'MenuBar', 'none');

    % UIWAIT makes DiveData wait for user response (see UIRESUME)
    % uiwait(handles.figure1);
end

% ********************************************************************************
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ********************************************************************************
function varargout = DiveData_OutputFcn(hObject, eventdata, handles)
    % Get default command line output from handles structure
    varargout{1} = handles.output;
end

% ********************************************************************************
% Initialize the GUI
% ********************************************************************************
function initialize_gui(fig_handle, handles, isreset)
    % If the metricdata field is present and the pushbuttonRefresh flag is false, it means
    % we are we are just re-initializing a GUI by calling it from the cmd line
    % while it is up. So, bail out as we dont want to pushbuttonRefresh the data.
    if isfield(handles, 'metricdata') && ~isreset
        return;
    end
    
    global NUM_PLOTS;
    global REGRESSION_ANALYSIS;
    global DIVE_ANALYSIS;
    global analysisMode;
    global gSavePlotsInProgress;
    global bRotate_3D;
    global gDisablePopups;
    global gSpeedOfSound;
    global gTrackLineColor;

    REGRESSION_ANALYSIS = 0;
    DIVE_ANALYSIS = 1;
    analysisMode = DIVE_ANALYSIS;
    NUM_PLOTS = 17;
    gSavePlotsInProgress = 0;
    bRotate_3D = 0;
    gDisablePopups = 0;
    gSpeedOfSound = 0;

    % Google Earth track line color.
    if ~isreset
        % Only reset this when app first start up.
        % Specifically do not reset it when Refresh button is pressed,
        % otherwise Refresh followed by Show Track will plot new
        % dive positions with different colors...
        gTrackLineColor = 1; % red - see createKML.m
    end
    
    % Disable buttons and other UI features that KUTI does not yet(?)
    % support.
	%set(handles.pushbuttonFlightViz,  'Enable', 'off');
	%set(handles.pushbuttonFaster,  'Enable', 'off');
	%set(handles.pushbuttonSlower,  'Enable', 'off');
	set(handles.checkboxUniqueFileName,  'Enable', 'off');
	set(handles.checkboxSaveDiveData,  'Enable', 'off');

    set(handles.pushbuttonRegressionVBD,  'Enable', 'off');
	set(handles.pushbuttonGenerateKML,  'Enable', 'off');
    set(handles.pushbuttonList_C_VBD_and_D_TGT,  'Enable', 'off');
    
    
    set(handles.figure1, 'MenuBar', 'figure');  % Display standard menu bar menus.

    %VERSION_NUMBER = 'Version 1.0-2';  % jf  Apr  3, 2014 - SG127, SG558, SG602-SG606    
    %VERSION_NUMBER = 'Version 1.0-4';  % jf  Jun  1, 2014 - SG506 (wetlabs .cnf file prefix change)    
    %VERSION_NUMBER = 'Version 1.0-5';  % jf  Jun  5, 2014 - WET Labs bb3 and fl3 support    
    %VERSION_NUMBER = 'Version 1.0-6';  % jf  Jun  6, 2014 - SG576 updated to UW/KUTI software, then deployed with iRobot .cnf files    
    %VERSION_NUMBER = 'Version 1.0-7';  % RTK Jun 13, 2014 - Added regression plots
    %VERSION_NUMBER = 'Version 1.0-8';   % jf  Jun 20, 2014 - Support plotting when .nc file contains data from multiple WETLabs pucks 
    VERSION_NUMBER = 'Version 1.1';   % wrr  Jun 24, 2014 - Since this is the version # displayed to the user, let's increment 
    set(handles.textVersion, 'String', VERSION_NUMBER);

    %jf addFoldersToPath();
    baseDir = fileparts(mfilename('fullpath'));
    
    if (exist(fullfile(baseDir,'DiveData.mat')))
        load (fullfile(baseDir,'DiveData.mat'));

        if isfield(settings, 'diveNumber') == 0
            settings.diveNumber = 1;
        end

        if isfield(settings, 'pathname')
            if (exist(settings.pathname, 'dir') == 7)
                [tmpId] = determineSeagliderId(settings.pathname);
                if (tmpId < 0)
                    set(handles.pushbuttonGeneratePlots, 'Enable', 'off');
                    set(handles.pushbuttonGenerateKML,  'Enable', 'off');
                    set(handles.pushbuttonRegressionVBD,  'Enable', 'off');
                    set(handles.pushbuttonList_C_VBD_and_D_TGT,  'Enable', 'off');
                    vizMsgBox('Please select a data folder. The folder you specified does not contain any SeaGlider .nc data files.', 'Required Data Files Missing Error', 'modal');
                else
                    settings.seagliderID = tmpId;
 
                    %   [calibData] = readCalibConstants(fullfile(settings.pathname, 'sg_calib_constants.m'));
                    %   settings.seagliderID = calibData.id_str;
                    %   settings.volMax = calibData.volmax;

                    set(handles.editPathname,    'string', settings.pathname);
                    set(handles.editSeagliderID, 'string', settings.seagliderID);
                    %    set(handles.editVolMax,      'string', num2str(settings.volMax));

                    % guidata(hObject,handles);
                    populateDiveNumberListBox(handles, settings.pathname, settings.seagliderID, settings.diveNumber);
                    populateSeagliderSerialNumberEditBox(handles, settings.pathname);

                    save(fullfile(baseDir,'DiveData.mat'), 'settings');

                    configureSettings(handles);
                end
            else
                set(handles.pushbuttonGeneratePlots, 'Enable', 'off');
                set(handles.pushbuttonGenerateKML,  'Enable', 'off');
                set(handles.pushbuttonRegressionVBD,  'Enable', 'off');
                set(handles.pushbuttonList_C_VBD_and_D_TGT,  'Enable', 'off');
               %vizMsgBox('The folder that you specified does not appear to have any or all of the required data files.', 'Required Data Files Missing Error', 'modal');
                vizMsgBox('The item that you specified does not appear to be a folder.', 'Required Data Files Missing Error', 'modal');
            end
        else
            set(handles.pushbuttonGeneratePlots, 'Enable', 'off');
            set(handles.pushbuttonGenerateKML,  'Enable', 'off');
            set(handles.pushbuttonRegressionVBD,  'Enable', 'off');
            set(handles.pushbuttonList_C_VBD_and_D_TGT,  'Enable', 'off');
            vizMsgBox('Please select a data folder. The folder you specified does not contain any or all of the required data files.', 'Required Data Files Missing Error', 'modal');
        end
    else
        configureSettings(handles);

        msgstr = sprintf('%s\n\n', 'Click the Browse button and browse to a folder that contains SeaGlider dive data files in NetCDF format.');
        msgstr = [msgstr sprintf('%s\n\n', 'A default settings file will automatically be generated for you if valid data is in the folder that you specify.')];
        vizMsgBox(msgstr, 'Creating Default Settings File', 'modal');
    end
    
    if ~isreset
        % msg gets old after awhile - only show it at startup, not reset
        msgstr = sprintf('%s\n', '1. Click the Browse button to browse and select the folder containing the Seaglider data to be analyzed.');
        msgstr = sprintf('%s\n%s\n', msgstr, '2. Using the dive number listbox, select the Dive Number(s) to be analyzed.');
        msgstr = sprintf('%s\n%s\n', msgstr, '3. Use the Plot Selection panel to control which plots are displayed.');
        msgstr = sprintf('%s\n%s', msgstr, '4. Click Generate Plots to generate and display the plots.');
        updateResults(handles, msgstr);
    end

    initGUI(handles);
    
    guidata(handles.figure1, handles);
end

function outputToMsgBox(str)
    msgstr = sprintf('%s\n', str);
    updateResults(handles, msgstr);
end

% ********************************************************************************
% Callback when the generate plots button is pressed (only roll, pitch
% regression supported)
% ********************************************************************************
function pushbuttonGeneratePlots_Callback(~, eventdata, handles)

    global NUM_PLOTS;
    global DIVE_ANALYSIS;
    global analysisMode;
 
    set(handles.pushbuttonSavePlots,  'Enable', 'off');
    set(handles.pushbuttonClosePlots, 'Enable', 'off');
    set(handles.pushbuttonTilePlots,  'Enable', 'off');

    analysisMode = DIVE_ANALYSIS;

    setAnnotatePlots(get(handles.checkboxAnnotatePlots, 'value'));

    settings.seagliderID = str2num(get(handles.editSeagliderID, 'string'));
    settings.vbdBias     = str2num(get(handles.editVBDBias, 'string'));
    settings.pathname    = get(handles.editPathname, 'string');
    settings.diveNumber  = getDiveNumber(handles);
    settings.volMax      = str2num(get(handles.editVolMax, 'string'));

    diveNumbers = get(handles.listboxDiveNumbers, 'String');
    selectedDiveNumbers = get(handles.listboxDiveNumbers, 'Value');
    diveNumberList = str2num(char(diveNumbers(selectedDiveNumbers)));

%     settings.diveNumber  = diveNumberList(1);

    % Determine which plots were selected
    settings.plotSelectionState(1:NUM_PLOTS) = 0;
    for i=1:NUM_PLOTS
        settings.plotSelectionState(i) = getPlotSelectionState(handles, i);
    end
    
    baseDir = fileparts(mfilename('fullpath'));

    save(fullfile(baseDir,'DiveData.mat'), 'settings');
    
    set(handles.textResults, 'String', '');

    %% Generate plots for selected dives...
    for i=1:length(diveNumberList)
        if (get(handles.checkboxAutoClosePlots, 'Value') == 1)
            autoClosePlots();
        end

        strResults = sprintf('\nProcessing data for dive %d\n  ', diveNumberList(i));
        updateResults(handles, strResults);

        [filenames] = buildFilenames(settings.pathname, settings.seagliderID, diveNumberList(i));

        [engHeader, fp] = read_header(filenames.p_eng);
        if (fp > 0)
            fclose(fp);
        end

        if (isfield(engHeader, 'columns'))
            reportColumnHeader(handles, engHeader);
        end

        err = [];
        if (get(handles.checkboxSaveDiveData, 'Value') == 1)
            [logp, eng, logInfo, err] = diveplot_func(settings.pathname, settings.pathname, diveNumberList(i), settings.vbdBias, settings, filenames, handles);

        else
            [logp, eng, logInfo, err] = diveplot_func(settings.pathname, [], diveNumberList(i), settings.vbdBias, settings, filenames, handles);
        end
        
        if (err == 1)
            return;
        end
        
        str = sprintf('Reference pitch gain = %g', logInfo.PITCH_GAIN);
        updateResults(handles, str);
   
               
        %% Generate regression plots of selected dives TODO: Update waitbar
        reference_c_vbd = str2num(get(handles.editreference_C_VBD, 'String'));
        regress_dives = diveNumberList;
        dd_regress_graphs();
         

        %% Auto tile plots        
        % IF WE ARE GENERATING PLOTS FOR MULTIPLE DIVES, THEN THERE IS NO NEED TO DISPLAY THE PLOTS.  THEY CAN BE
        % REVIEWED IN AN IMAGE VIEWER WHEN THEY GET SAVED.
        % FUTURE FEATURE: PUT A "CONTINUE" OPTION IN THAT LETS THE USER REVIEW THE PLOTS BEFORE GENERATING PLOTS FOR
        % THE NEXT DIVE IN THE LIST.
        if (get(handles.checkboxAutoTilePlots, 'value') == 1) && (i == length(diveNumberList))
            autoTilePlots();
            hWaitbar = waitbar(1);
            close(hWaitbar);
            %return;
        end
        
        %if (isempty(logp) || isempty(eng))
            %! jf - non-issue - using .nc files instead
            %!strResults = sprintf('The engineering and/or log file for Dive %d are empty; therefore Dive %d cannot be visualized.', diveNumberList(i), diveNumberList(i));
            %!disp(strResults);
            %!updateResults(handles, strResults);
        %else
        if get(handles.checkboxAutoSavePlots, 'value') == 1
            strResults = sprintf('Saving results for dive %d', diveNumberList(i));
            updateResults(handles, strResults);
            [filenames] = buildFilenames(settings.pathname, settings.seagliderID, diveNumberList(i));

            set(handles.pushbuttonTilePlots, 'Value', 0);
            guidata(handles.pushbuttonTilePlots, handles);
            drawnow();
            savePlots([], fullfile(settings.pathname, ['p',filenames.root,'_DivePlots']));
            set(handles.pushbuttonTilePlots, 'Value', 1);

            guidata(handles.pushbuttonTilePlots, handles);

            drawnow();
        end
        %end
        
    end
    

    set(handles.pushbuttonSavePlots,  'Enable', 'on');
    set(handles.pushbuttonClosePlots, 'Enable', 'on');
    set(handles.pushbuttonTilePlots,  'Enable', 'on');

    hWaitbar = waitbar(1);
    close(hWaitbar);
end


% ********************************************************************************
% Perform regression of VBD data
% ********************************************************************************
function doRegressionVBD(handles)
% handles - Common GUI handles
%
    global REGRESSION_ANALYSIS;
    global analysisMode;

    baseDir = fileparts(mfilename('fullpath'));

    if exist(fullfile(baseDir,'DiveData.mat'))
        if (get(handles.checkboxAutoClosePlots, 'Value') == 1)
            autoClosePlots();
        end
        
        strResults = sprintf('\n\n%s', 'Regression VBD');
        updateResults(handles, strResults);

        analysisMode = REGRESSION_ANALYSIS;

        load (fullfile(baseDir,'DiveData.mat'));

        settings.pathname       = get(handles.editPathname, 'string');
        diveNumbers             = get(handles.listboxDiveNumbers, 'String');
        selectedDiveNumbers     = get(handles.listboxDiveNumbers, 'Value');
        diveNumberList          = str2num(char(diveNumbers(selectedDiveNumbers)));
        reference_c_vbd         = str2num(get(handles.editreference_C_VBD, 'String'));
        settings.seagliderID    = str2num(get(handles.editSeagliderID, 'string'));
        settings.volMax         = str2num(get(handles.editVolMax, 'string'));

        save(fullfile(baseDir,'DiveData.mat'), 'settings');

        strResults = sprintf('Reference C VBD = %10.4f', reference_c_vbd);
        updateResults(handles, strResults);

        [vbdBias, w_rms_final] = dd_regress_vbd(handles, reference_c_vbd, [diveNumberList], settings);
        
        settings.vbdBias = vbdBias;
        save(fullfile(baseDir,'DiveData.mat'), 'settings');

        %strResults = sprintf('vbd Bias = %10.4f\nrms(final) = %10.4f', vbdBias, w_rms_final);
        %updateResults(handles, strResults);
        if get(handles.checkboxAutoSavePlots, 'value') == 1
            %strResults = sprintf('Saving regression VBD plots...');
            %updateResults(handles, strResults);

            [filenames] = buildFilenames(settings.pathname, settings.seagliderID, diveNumberList(1));

            subFolder = sprintf('p%3.3d%4.4d-p%3.3d%4.4d_RegressionPlots',settings.seagliderID, diveNumberList(1), settings.seagliderID, diveNumberList(end));
            savePlots([], fullfile(settings.pathname, subFolder));
        end

        if (get(handles.checkboxAutoTilePlots, 'value') == 1)
            autoTilePlots();
        end
   else
        configureSettings(handles);

        msgstr = sprintf('%s\n\n', 'A default settings files does not exist for the dive data analysis program.');
        msgstr = [msgstr sprintf(' %s', 'One will be created for you.')];
        msgstr = [msgstr sprintf('%s\n\n', 'You must browse to a folder that contains valid and complete dive data.')];
        msgstr = [msgstr sprintf('%s\n\n', 'A default settings file will automatically be generated for you if valid data is in the folder that you specify.')];
        vizMsgBox(msgstr, 'CreateMode', 'modal');
    end
end


% ********************************************************************************
% Callback when the regression VBD button is pressed
% ********************************************************************************
function pushbuttonRegressionVBD_Callback(hObject, eventdata, handles)
    doRegressionVBD(handles);
end

% ********************************************************************************
% GUI Callback.
% ********************************************************************************
function editDiveNumber_CreateFcn(hObject, eventdata, handles)
    usewhitebg = 1;
    if usewhitebg
        set(hObject,'BackgroundColor','white');
    else
        set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
    end
end

% ********************************************************************************
% GUI Callback.
% ********************************************************************************
function editDiveNumber_Callback(hObject, eventdata, handles)
    density = str2double(get(hObject, 'String'));
    if isnan(density)
        set(hObject, 'String', 0);
        errordlg('Input must be a number','Error');
    end

    % Save the new editDiveNumber value
    handles.metricdata.density = density;
    guidata(hObject,handles)
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function editVBDBias_CreateFcn(hObject, eventdata, handles)
    usewhitebg = 1;
    if usewhitebg
        set(hObject,'BackgroundColor','white');
    else
        set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
    end
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function editVBDBias_Callback(hObject, eventdata, handles)
    volume = str2double(get(hObject, 'String'));
    if isnan(volume)
        set(hObject, 'String', 0);
        errordlg('Input must be a number','Error');
    end

    % Save the new editVBDBias value
    handles.metricdata.volume = volume;
    %set(handles.editVBDBias, 'Value', volume);
	set(handles.editVBDBias, 'string', num2str(volume));

    guidata(hObject,handles)
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function pushbuttonRefresh_Callback(hObject, eventdata, handles)
    initialize_gui(gcbf, handles, true);
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function unitgroup_SelectionChangeFcn(hObject, eventdata, handles)
    if (hObject == handles.english)
        set(handles.text4, 'String', 'lb/cu.in');
        set(handles.textVBDBiasUnits, 'String', 'cu.in');
        set(handles.text6, 'String', 'lb');
    else
        set(handles.text4, 'String', 'kg/cu.m');
        set(handles.textVBDBiasUnits, 'String', 'cu.m');
        set(handles.text6, 'String', 'kg');
    end
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function editPathname_Callback(hObject, eventdata, handles)

    [pathname] = get(handles.editPathname, 'String');
    baseDir = fileparts(mfilename('fullpath'));

    if (exist(pathname, 'dir') == 7)
        if (exist(fullfile(baseDir,'DiveData.mat')) > 0)
            load (fullfile(baseDir,'DiveData.mat'));
        end
        settings.pathname = pathname;
        %[calibData] = readCalibConstants(fullfile(settings.pathname, 'sg_calib_constants.m'));
        %settings.seagliderID = calibData.id_str;
        %settings.volMax = calibData.volmax;
        [tmpId] = determineSeagliderId(settings.pathname);
        if (tmpId < 0)
            set(handles.pushbuttonGeneratePlots, 'Enable', 'off');
            vizMsgBox('Please select a data folder. The folder you specified does not contain any SeaGlider .nc data files.', 'Required Data Files Missing Error', 'modal');
        else
            settings.seagliderID = tmpId;


            set(handles.editPathname,    'string', settings.pathname);
            set(handles.editSeagliderID, 'string', settings.seagliderID);
            %set(handles.editVolMax, 'string', num2str(settings.volMax));

            % guidata(hObject,handles);
            populateDiveNumberListBox(handles, settings.pathname, settings.seagliderID, settings.diveNumber);
            populateSeagliderSerialNumberEditBox(handles, settings.pathname);

            save(fullfile(baseDir,'DiveData.mat'), 'settings');
        end
    else
        set(handles.pushbuttonGeneratePlots, 'Enable', 'off');
        set(handles.pushbuttonGenerateKML,  'Enable', 'off');
        set(handles.pushbuttonRegressionVBD,  'Enable', 'off');
        set(handles.pushbuttonList_C_VBD_and_D_TGT,  'Enable', 'off');
        vizMsgBox('The item you selected does not appear to be a folder.', 'CreateMode', 'modal');
    end
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function editPathname_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function pushbuttonBrowseForDataFiles_Callback(hObject, eventdata, handles)
    global   seagliderID;
    global   gTrackLineColor;
    
    [pathname] = uigetdir( get(handles.editPathname,'string'), 'Select folder containing Seaglider .nc data files.' );
    baseDir = fileparts(mfilename('fullpath'));

    if (pathname ~= 0)
        if (exist(pathname, 'dir') == 7)
            if (exist(fullfile(baseDir,'DiveData.mat')) > 0)
                load (fullfile(baseDir,'DiveData.mat'));
            end

            [tmpId] = determineSeagliderId(pathname);
            if (tmpId < 0)
                set(handles.pushbuttonGeneratePlots, 'Enable', 'off');
                set(handles.pushbuttonGenerateKML,  'Enable', 'off');
                set(handles.pushbuttonRegressionVBD,  'Enable', 'off');
                set(handles.pushbuttonList_C_VBD_and_D_TGT,  'Enable', 'off');
                vizMsgBox('Please select a data folder. The folder you specified does not contain any SeaGlider .nc data files.', 'Required Data Files Missing Error', 'modal');
            else
                settings.seaGliderID = tmpId;
                
                % Update track-line color.
                gTrackLineColor = gTrackLineColor + 1;
                % 8 colors supported - if > 8 vehicles, colors will be reused
                % see createKML.m
                if (gTrackLineColor > 8)
                     gTrackLineColor = 1;
                end
 
            
            %[calibData] = readCalibConstants(fullfile(pathname, 'sg_calib_constants.m')); % THIS IS A KLUDGE TO GET THE SEAGLIDER ID
            %if (isempty(calibData) == 0)
                settings.pathname = pathname;
                %settings.seagliderID = calibData.id_str;
                %settings.volMax = calibData.volmax;

                save(fullfile(baseDir,'DiveData.mat'), 'settings');

                set(handles.editSeagliderID, 'string', settings.seaGliderID);
                set(handles.editPathname,    'string', pathname);
                %set(handles.editVolMax,      'string', num2str(settings.volMax));

                populateDiveNumberListBox(handles, settings.pathname, settings.seagliderID, []);
                populateSeagliderSerialNumberEditBox(handles, settings.pathname);
            %else
               % vizMsgBox('The Seaglider Calibration Constants file is either missing or corrupt. The data in this mission cannot be viewed.', 'CreateMode', 'modal');
            end
        else
            set(handles.pushbuttonGeneratePlots, 'Enable', 'off');
            set(handles.pushbuttonGenerateKML,  'Enable', 'off');
            set(handles.pushbuttonRegressionVBD,  'Enable', 'off');
            set(handles.pushbuttonList_C_VBD_and_D_TGT,  'Enable', 'off');
            vizMsgBox('The folder that you specified does not appear to have any or all of the required data files.', 'CreateMode', 'modal');
        end
    end
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function editSeagliderID_Callback(hObject, eventdata, handles)
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function editSeagliderID_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function editSalinityMin_Callback(hObject, eventdata, handles)
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function editSalinityMin_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function editSalinityMax_Callback(hObject, eventdata, handles)
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function editSalinityMax_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function editPressureMin_Callback(hObject, eventdata, handles)
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function editPressureMin_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function editPressureMax_Callback(hObject, eventdata, ~)
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function editPressureMax_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function editDensityMin_Callback(hObject, eventdata, handles)
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function editDensityMin_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function editDensityMax_Callback(hObject, eventdata, handles)
end

% ********************************************************************************
% ********************************************************************************
function editDensityMax_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function editTemperatureMin_Callback(hObject, eventdata, handles)
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function editTemperatureMin_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function editTemperatureMax_Callback(hObject, eventdata, handles)
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function editTemperatureMax_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function pushbuttonSavePlots_Callback(hObject, eventdata, handles)
    global REGRESSION_ANALYSIS;
    global DIVE_ANALYSIS;
    global analysisMode;

    baseDir = fileparts(mfilename('fullpath'));

    load(fullfile(baseDir,'DiveData.mat'));

    if (analysisMode == REGRESSION_ANALYSIS)
        subFolder = sprintf('p%3.3d%4.4d-p%3.3d%4.4d_RegressionPlots',settings.seagliderID, settings.diveNumber(1), settings.seagliderID, settings.diveNumber(end));
        savePlots([], fullfile(settings.pathname, subFolder));
    else
        for i=1:length(settings.diveNumber)
            [filenames] = buildFilenames(settings.pathname, settings.seagliderID, settings.diveNumber(i));
            dirName = fullfile(settings.pathname, ['p',filenames.root,'_DivePlots']);
            savePlots([], dirName); %fullfile(settings.pathname, ['p',filenames.root,'_DivePlots']));
        end
    end

end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function axesThumbPlot_2_ButtonDownFcn(hObject, eventdata, handles)
    % hObject    handle to axesThumbPlot_2 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function listboxDiveNumbers_Callback(hObject, eventdata, handles)
    % hObject    handle to listboxDiveNumbers (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = get(hObject,'String') returns listboxDiveNumbers contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from listboxDiveNumbers
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function listboxDiveNumbers_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to listboxDiveNumbers (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: listbox controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

% ********************************************************************************
% Populates the dive number list box.
% ********************************************************************************
function populateDiveNumberListBox(handles, pathname, seagliderID, defaultDiveNumber)

    diveNumbers = [];

    hWaitbar = waitbar(0, 'Fetching list of available dives...');

    tempSeagliderID = str2num(get(handles.editSeagliderID, 'String'));
    %kuti - jf - only need .nc file
    %fileFolder = sprintf('p%3.3d*.log', tempSeagliderID);
    fileFolder = sprintf('p%3.3d*.nc', tempSeagliderID);
    filelist = dir(fullfile(pathname, fileFolder));

    if (isempty(filelist) == 0)

        waitbar(0, hWaitbar, 'Validating existence of dive data files...');

        set(handles.listboxDiveNumbers, 'String', []);

        % EXTRACT DIVE NUMBER IF BOTH THE .LOG AND .ENG FILES ARE PRESENT FOR THE DIVE.
        % kuti - jf - only need .nc file
        for i=1:length(filelist)
            waitbar(i/length(filelist), hWaitbar);
            [token, remain] = strtok(filelist(i).name, '.');
            %if (exist(fullfile(pathname,[token,'.eng'])) > 0)
                diveNumbers{end+1} = num2str(str2num(token(end-3:end)));
            %end
            
            waitbar(i/length(filelist));
        end

        diveNumbers = strtrim(cellstr(num2str(sort(str2num(char(diveNumbers))))));

        set(handles.listboxDiveNumbers, 'String', diveNumbers);

        set(handles.listboxDiveNumbers, 'Value', 1);

        set(handles.pushbuttonGeneratePlots, 'Enable', 'on');
        set(handles.pushbuttonGenerateKML,  'Enable', 'on');
        set(handles.pushbuttonRegressionVBD,  'Enable', 'on');
        set(handles.pushbuttonList_C_VBD_and_D_TGT,  'Enable', 'on');

        set(handles.listboxDiveNumbers,'Max', length(diveNumbers));
    else
        set(handles.listboxDiveNumbers, 'String', []);
        set(handles.pushbuttonGeneratePlots, 'Enable', 'off');
        set(handles.pushbuttonGenerateKML,  'Enable', 'off');
        set(handles.pushbuttonRegressionVBD,  'Enable', 'off');
        set(handles.pushbuttonList_C_VBD_and_D_TGT,  'Enable', 'off');
        vizMsgBox('The folder that you specified does not appear to have any or all of the required data files.', 'CreateMode', 'modal');
    end
    
    close(hWaitbar);
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function populateSeagliderSerialNumberEditBox(handles, pathname)
    % kuti - jf - only need .nc file
    %filelist = dir(fullfile(pathname, 'p*.log'));
    filelist = dir(fullfile(pathname, 'p*.nc'));

    if (isempty(filelist) == 0)
        [token, remain] = strtok(filelist(1).name, '.');
        SeagliderID = token(2:4);
        set(handles.editSeagliderID, 'String', SeagliderID);

        set(handles.pushbuttonGeneratePlots, 'Enable', 'on');
        set(handles.pushbuttonGenerateKML,  'Enable', 'on');
        set(handles.pushbuttonRegressionVBD,  'Enable', 'on');
        set(handles.pushbuttonList_C_VBD_and_D_TGT,  'Enable', 'on');
    else
        set(handles.editSeagliderID, 'String', []);
        set(handles.pushbuttonGeneratePlots, 'Enable', 'off');
        set(handles.pushbuttonGenerateKML,  'Enable', 'off');
        set(handles.pushbuttonRegressionVBD,  'Enable', 'off');
        set(handles.pushbuttonList_C_VBD_and_D_TGT,  'Enable', 'off');
    end
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function pushbuttonReadHeader_Callback(hObject, eventdata, handles)
    [filenames] = buildFilenames(getPathname(handles), getSeagliderID(handles), getDiveNumber(handles));

    [engHeader, fp] = read_header(filenames.p_eng);
    fclose(fp);

    if (isfield(engHeader, 'columns'))
        reportColumnHeader(handles, engHeader);
    end
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function [diveNumber] = getDiveNumber(handles)
    diveNumbers = [get(handles.listboxDiveNumbers, 'String')];
    diveNumberListboxIndex = get(handles.listboxDiveNumbers, 'Value');
    diveNumber = str2num(char(diveNumbers(diveNumberListboxIndex)));
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function [seagliderID] = getSeagliderID(handles)
    seagliderID = str2num(get(handles.editSeagliderID, 'string'));
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function [pathname] = getPathname(handles)
    pathname = get(handles.editPathname, 'string');
end

% ********************************************************************************
% GUI callback - Plot 1 - Operational - Vehicle Attitude (diveplot)
% ********************************************************************************
function checkboxSelectPlot_1_Callback(hObject, eventdata, handles)
    savePlotSelectionState(handles, get(hObject, 'Tag'));
end

% ********************************************************************************
% GUI callback - Plot 2 - Operational - Path Through the Water.
% ********************************************************************************
function checkboxSelectPlot_2_Callback(hObject, eventdata, handles)
    savePlotSelectionState(handles, get(hObject, 'Tag'));
end

% ********************************************************************************
% GUI callback - Plot 3 - Operational - Vertical Velocities.
% ********************************************************************************
function checkboxSelectPlot_3_Callback(hObject, eventdata, handles)
    savePlotSelectionState(handles, get(hObject, 'Tag'));
end

% ********************************************************************************
% GUI callback - Plot 4 - Operational - Pitch Control
% ********************************************************************************
function checkboxSelectPlot_4_Callback(hObject, eventdata, handles)
    savePlotSelectionState(handles, get(hObject, 'Tag'));
end

% ********************************************************************************
% GUI callback - Plot 5 - Operational - Roll Control
% ********************************************************************************
function checkboxSelectPlot_5_Callback(hObject, eventdata, handles)
    savePlotSelectionState(handles, get(hObject, 'Tag')); 
end

% ********************************************************************************
% GUI callback - Plot 6 - Operational - Turn Rate
% ********************************************************************************
function checkboxSelectPlot_6_Callback(hObject, eventdata, handles)
    savePlotSelectionState(handles, get(hObject, 'Tag')); 
end

% ********************************************************************************
% GUI callback - Plot 7 - Operational - Buoyancy
% ********************************************************************************
function checkboxSelectPlot_7_Callback(hObject, eventdata, handles)
    savePlotSelectionState(handles, get(hObject, 'Tag')); 
end

% ********************************************************************************
% GUI callback - Plot 8 - Operational - VBD Configuration
% ********************************************************************************
function checkboxSelectPlot_8_Callback(hObject, eventdata, handles)
    %savePlotSelectionState(handles, get(hObject, 'Tag')); 
end

% ********************************************************************************
% GUI callback - Plot 9 - Sensor Data - Temperature and Salinity
% ********************************************************************************
function checkboxSelectPlot_9_Callback(hObject, eventdata, handles)
    savePlotSelectionState(handles, get(hObject, 'Tag')); 
end

% ********************************************************************************
% GUI callback - Plot 10 - Sensor Data - Temperature vs. Salinity
% ********************************************************************************
function checkboxSelectPlot_10_Callback(hObject, eventdata, handles)
    savePlotSelectionState(handles, get(hObject, 'Tag')); 
end

% ********************************************************************************
% GUI callback - Plot 11  - Sensor Data - Dissolved Oxygen
% ********************************************************************************
function checkboxSelectPlot_11_Callback(hObject, eventdata, handles)
    savePlotSelectionState(handles, get(hObject, 'Tag'));
end

% ********************************************************************************
% GUI callback - Plot 12 - Sensor Data - Fluorescence
% ********************************************************************************
function checkboxSelectPlot_12_Callback(hObject, eventdata, handles)
    savePlotSelectionState(handles, get(hObject, 'Tag'));
end

% ********************************************************************************
% GUI callback - Plot 13 - Sensor Data - Backscatter
% ********************************************************************************
function checkboxSelectPlot_13_Callback(hObject, eventdata, handles)
    savePlotSelectionState(handles, get(hObject, 'Tag'));
end

% ********************************************************************************
% GUI callback - Plot 14 - Sensor Data - PAR
% ********************************************************************************
function checkboxSelectPlot_14_Callback(hObject, eventdata, handles)
    savePlotSelectionState(handles, get(hObject, 'Tag'));
end

% ********************************************************************************
% GUI callback - Plot 15 - Regression - Pitch
% ********************************************************************************
function checkboxSelectPlot_15_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxSelectPlot_15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxSelectPlot_15
end

% ********************************************************************************
% GUI callback - Plot 16 - Regression - Roll
% ********************************************************************************
function checkboxSelectPlot_16_Callback(hObject, eventdata, handles)
end

% ********************************************************************************
% GUI callback - Plot 17 - Regression - VBD
% ********************************************************************************
function checkboxSelectPlot_17_Callback(hObject, eventdata, handles)
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function savePlotSelectionState(handles, uiControlTag)
    baseDir = fileparts(mfilename('fullpath'));

    load(fullfile(baseDir,'DiveData.mat'));
        strIndex = uiControlTag( (strfind(uiControlTag, '_')+1):end);
        %strEval = sprintf( 'settings.Plot(%s).SelectionState = get(handles.checkboxSelectPlot_%s, ''Value'');', strIndex, strIndex);
        strEval = sprintf( 'settings.plotSelectionState(%s) = get(handles.checkboxSelectPlot_%s, ''Value'');', strIndex, strIndex);
        eval(strEval);
    save(fullfile(baseDir,'DiveData.mat'), 'settings');
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function checkboxAutoSavePlots_Callback(hObject, ~, handles)
    saveAutoSavePlotsSetting(hObject, handles);
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function saveAutoSavePlotsSetting(hObject, handles)
    baseDir = fileparts(mfilename('fullpath'));

    load(fullfile(baseDir,'DiveData.mat'));
        strEval = sprintf( 'settings.autoSavePlots = get(handles.%s, ''Value'');', get(hObject, 'Tag'));
        eval(strEval);
    save(fullfile(baseDir,'DiveData.mat'), 'settings');
end


% ********************************************************************************
% GUI callback
% ********************************************************************************
function pushbuttonClosePlots_Callback(hObject, eventdata, handles)
    global gSavePlotsInProgress;

    if (gSavePlotsInProgress == 0)
        autoClosePlots();
    end
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function checkboxAutoClosePlots_Callback(hObject, eventdata, handles)
    saveAutoClosePlotsSetting(hObject, handles);
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function saveAutoClosePlotsSetting(hObject, handles)
    baseDir = fileparts(mfilename('fullpath'));

    load(fullfile(baseDir,'DiveData.mat'));
        strEval = sprintf( 'settings.autoClosePlots = get(handles.%s, ''Value'');', get(hObject, 'Tag'));
        eval(strEval);
    save(fullfile(baseDir,'DiveData.mat'), 'settings');
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function pushbuttonTilePlots_Callback(hObject, eventdata, handles)
    autoTilePlots();
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function checkboxAutoTilePlots_Callback(hObject, eventdata, handles)
    saveAutoTilePlotsSetting(hObject, handles);
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function saveAutoTilePlotsSetting(hObject, handles)
    baseDir = fileparts(mfilename('fullpath'));

    load(fullfile(baseDir,'DiveData.mat'));
        strEval = sprintf( 'settings.autoTilePlots = get(handles.%s, ''Value'');', get(hObject, 'Tag'));
        eval(strEval);
    save(fullfile(baseDir,'DiveData.mat'), 'settings');
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function autoTilePlots()
    [plotFigureHandles, plotFigureNames] = getPlotFigureHandles();

    if (isempty(plotFigureHandles) == 0)
        tilePlots(plotFigureHandles, [], 5);
    end
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function createDefaultSettingsFile()

    global NUM_PLOTS;

    % IF THE SETTINGS FILE DOES NOT HAVE THE PLOT SELECTION FIELD DEFINED, DEFINE ONE AND 
    for i=1:NUM_PLOTS
        strEval = sprintf( 'set(handles.checkboxSelectPlot_%d, ''Value'', 1)', i);
        eval(strEval);

        strEval = sprintf( 'settings.Plot(%d).SelectionState = get(handles.checkboxSelectPlot_%d, ''Value'');', i, i);
        eval(strEval);
    end

    settings.analysisMode   = 1;

    settings.autoSavePlots  = 1;
    settings.autoClosePlots = 1;
    settings.autoTilePlots  = 1;

    baseDir = fileparts(mfilename('fullpath'));
    save(fullfile(baseDir,'DiveData.mat'), 'settings');
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function configureSettings(handles)

    global NUM_PLOTS;

    baseDir = fileparts(mfilename('fullpath'));

    try
        load(fullfile(baseDir,'DiveData.mat'));
    catch
        fprintf('Generating new configuration files (DiveData.mat).\n');
    end

    if exist('settings') == 0
        settings = [];
        settings.seagliderID = 601;
        % initialize plot selection state - all ON.
        for i=1:NUM_PLOTS
            settings.plotSelectionState(i) = 1;
        end
        % except vbd config plot is off (always).
        settings.plotSelectionState(8) = 0;
    end
    
    % PLOT SELECTIONS
    if isfield(settings, 'plotSelectionState') == 1
        for i=1:NUM_PLOTS
            % Set plot selection checkbox value based on saved settings value
            strEval = sprintf( 'isfield(handles, ''checkboxSelectPlot_%d'')', i);
            if (eval(strEval) == 1)
                strEval = sprintf( 'set(handles.checkboxSelectPlot_%d, ''Value'', %d)', i, settings.plotSelectionState(i));
                eval(strEval);
            end
        end
        % except vbd config plot is off (always).
        settings.plotSelectionState(8) = 0;
    end
    
    % how to disable one of the plots
    %set(handles.checkboxSelectPlot_8, 'Enable', 'off');
    
    [settings] = disableProvisionalPlots(handles, settings);

    % AUTO SAVE PLOTS
    if (isfield(settings, 'autoSavePlots'))
        set(handles.checkboxAutoSavePlots,  'Value', settings.autoSavePlots);
    else
    	settings.autoSavePlots = 1;
        set(handles.checkboxAutoSavePlots,  'Value', settings.autoSavePlots);
    end

    % AUTO CLOSE PLOTS
    if (isfield(settings, 'autoClosePlots'))
        set(handles.checkboxAutoClosePlots,  'Value', settings.autoClosePlots);
    else
    	settings.autoClosePlots = 1;
        set(handles.checkboxAutoClosePlots,  'Value', settings.autoClosePlots);
    end

    % AUTO TILE PLOTS
    if (isfield(settings, 'autoTilePlots'))
        set(handles.checkboxAutoTilePlots,  'Value', settings.autoTilePlots);
    else
    	settings.autoTilePlots = 1;
        set(handles.checkboxAutoTilePlots,  'Value', settings.autoTilePlots);
    end

    % VBD BIAS
    if (isfield(settings, 'vbdBias'))
        set(handles.editVBDBias, 'string', num2str(settings.vbdBias));
    else
        settings.vbdBias = 0.5;
        set(handles.editVBDBias, 'string', num2str(settings.vbdBias));
    end

    save(fullfile(baseDir,'DiveData.mat'), 'settings');
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function pushbuttonDiscontinuities_Callback(hObject, eventdata, handles)
    base_file = sprintf('p%d%04d', getSeagliderID(handles), getDiveNumber(handles));

    eng_file = fullfile(get(handles.editPathname, 'string'), strcat(base_file, '.eng'));
    log_file = fullfile(get(handles.editPathname, 'string'), strcat(base_file, '.log'));

    if length(dir(log_file)) == 0
       disp('log file has zero length');
       return
    end

    logp = read_log(log_file);

    eng = read_eng(eng_file);

    [filenames] = buildFilenames(get(handles.editPathname, 'string'), getSeagliderID(handles), getDiveNumber(handles));
    
    if (exist(filenames.ppc_a_eng) > 0) || (exist(filenames.ppc_b_eng) > 0)
        eng.GPCTD = read_GPCTD_data(filenames.ppc_a_eng, filenames.ppc_b_eng);
        ctdData = eng.GPCTD;
%     elseif (exist(filenames.ppc_a_eng) > 0)
%         eng.GPCTD = read_GPCTD_data(filenames.ppc_a_eng, filenames.ppc_b_eng);
%         ctdData = eng.GPCTD;
%     elseif (exist(filenames.ppc_b_eng) > 0)
%         eng.GPCTD = read_GPCTD_data(filenames.ppc_a_eng, filenames.ppc_b_eng);
%         ctdData = eng.GPCTD;
    else
        ctdData = [];
    end

    GPCTD = [];
    if (isempty(ctdData) == 0)
%         [indexes] = subsampleCTD_Data(eng, ctdData);
        
        if (length(ctdData.downCast.raw) > 0)
            GPCTD = [GPCTD ctdData.downCast.raw'];
        end
        if (length(ctdData.upCast.raw) > 0)
            GPCTD = [GPCTD, ctdData.upCast.raw'];
        end
        GPCTD = [GPCTD]';
    end

    strLegendIndex = 0;
    strLegend = {};
    figure; clf;
    hold on;
    plot (1:length(ctdData.downCast.raw(1:end,1))-1, ctdData.downCast.raw(2:end,1)-ctdData.downCast.raw(1:end-1,1), 'red');
    strLegendIndex = strLegendIndex + 1;
    strLegend{strLegendIndex} = sprintf('Down cast: Pressure');

    plot (1:length(ctdData.downCast.raw(1:end,2))-1, ctdData.downCast.raw(2:end,2)-ctdData.downCast.raw(1:end-1,2), 'green');
    strLegendIndex = strLegendIndex + 1;
    strLegend{strLegendIndex} = sprintf('Down cast: Temperature');

    plot (1:length(ctdData.downCast.raw(1:end,3))-1, ctdData.downCast.raw(2:end,3)-ctdData.downCast.raw(1:end-1,3), 'blue');
    strLegendIndex = strLegendIndex + 1;
    strLegend{strLegendIndex} = sprintf('Down cast: Conductivity');

    xlabel( 'Sample Index' );
    ylabel( 'First Derivative');

    if(isempty(strLegend) == 0)
        legend( strLegend );
    end

    hold off;
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function pushbuttonSelectAllPlots_Callback(hObject, eventdata, handles)
    selectAllPlots(handles);
    disableProvisionalPlots(handles);
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function selectAllPlots(handles)

    global  NUM_PLOTS;

    for i=1:NUM_PLOTS
        strEval = sprintf( 'isfield(handles, ''checkboxSelectPlot_%d'')', i);
        if (eval(strEval) == 1)
            strEval = sprintf( 'set(handles.checkboxSelectPlot_%d, ''Value'', 1)', i);
            eval(strEval);
        end
    end
    
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function pushbuttonDeselectAllPlots_Callback(hObject, eventdata, handles)
    deselectAllPlots(handles);
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function deselectAllPlots(handles)

    global  NUM_PLOTS;

    for i=1:NUM_PLOTS
        strEval = sprintf( 'isfield(handles, ''checkboxSelectPlot_%d'')', i);
        if (eval(strEval) == 1)
            strEval = sprintf( 'set(handles.checkboxSelectPlot_%d, ''Value'', 0)', i);
            eval(strEval);
       end
    end
     
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function pushbuttonSoundSpeed_Callback(hObject, eventdata, handles)
    soundSpeed = calcSoundSpeed(temp_C, depth, salinity);
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function figure1_CloseRequestFcn(hObject, eventdata, handles)
    autoClosePlots();

    delete(hObject);
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function editreference_C_VBD_Callback(hObject, eventdata, handles)
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function editreference_C_VBD_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
%function checkboxRegression_VBD_ABC_Callback(hObject, eventdata, handles)
%    if (get(handles.checkboxRegression_VBD_ABC, 'Value') == 1)
%        if (get(handles.checkboxRegression_VBD_Therm, 'Value') == 1)
%            set(handles.checkboxRegression_VBD_Therm, 'Value', 0);
%        end
%    end
%end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function checkboxSaveDiveData_Callback(hObject, eventdata, handles)
end

% ********************************************************************************
% GUI callback
% ********************************************************************************
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%     hAxes       = findall(hFig,'type','axes');
%     hText       = findall(hAxes,'type','text');
%     hUIControls = findall(gcf,'type','uicontrol');
% 
%     set([hAxes; hText; hUIControls],'units','normalized','fontunits','normalized');

%     set(0,'Units','Normalized')
%     set(fig,'Units','Normalized')
%     set(get(fig,'Children'),'Units','Normalized')
%     old_pos=get(fig,'Position')

end

% ****************************************************************************************************
% GUI callback
% ****************************************************************************************************
function editVolMax_Callback(hObject, eventdata, handles)
end

% ****************************************************************************************************
% GUI callback
% ****************************************************************************************************
function editVolMax_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

% ****************************************************************************************************
% GUI callback
% ****************************************************************************************************
function [calibParam] = XgetCalibParam(calibData, strParam)
    calibParam = 0.0;

    idx = strfind(calibData, strParam);
    idx = find(not(cellfun('isempty', idx)));
    if isempty(idx) == 0
        eval(calibData{idx});
        calibParam = eval('volmax');
    end
end

% ******************************************************************************
% GUI callback
% ******************************************************************************
function pushbuttonList_C_VBD_and_D_TGT_Callback(hObject, eventdata, handles)
    diveNumbers = get(handles.listboxDiveNumbers, 'String');
    selectedDiveNumbers = get(handles.listboxDiveNumbers, 'Value');
    diveNumberList = str2num(char(diveNumbers(selectedDiveNumbers)));

    baseDir = fileparts(mfilename('fullpath'));
    
    if (exist(fullfile(baseDir,'DiveData.mat')))
        load(fullfile(baseDir,'DiveData.mat'));
        tempSeagliderID = str2num(get(handles.editSeagliderID, 'String'));

        vbdPitchList = [];
        vbdPitchList = sprintf('\n\nDIVE   C_VBD   D_TGT');
        vbdPitchList = sprintf('%s\n%s', vbdPitchList, '=========================');

        h = waitbar(0, 'Please wait while ALL dive data log files are searched...');
 
        for i=1:length(diveNumberList)
            waitbar(i/length(diveNumberList));
           
            % read the .nc file
            id_dive = sprintf('%03d%04d',tempSeagliderID,diveNumberList(i));
            [loginfo,eng,results] = dd_get_dive_data(id_dive,1,settings.pathname);
            %sg_cal = results.sg_calib_constants;

            % extract desired fields and add to output
            vbdPitchList = sprintf('%s\n%5d      %6d     %6d', vbdPitchList, diveNumberList(i), loginfo.C_VBD, loginfo.D_TGT);
        end
        
        updateResults(handles, char(vbdPitchList));
                
        close(h);
    end
end

% ******************************************************************************
% GUI callback
% ******************************************************************************
function pushbuttonClear_Callback(hObject, eventdata, handles)
    % Clear contents of big message window at bottom of DiveData main
    % window.
	set(handles.textResults, 'String', '');
end

% ******************************************************************************
% GUI callback
% ******************************************************************************
function reportColumnHeader(handles, engHeader)
    strMsg = strcat('Data Columns', 10, '============');
    
    strCol = char([cellstr(engHeader.columns)']);
    [m, n] = size(strCol);

    strMsg = sprintf('%s\n', strMsg);
    for i=1:m
        strMsg = sprintf('%s %s,', strMsg, strCol(i,:));
    end
    strMsg(end) = [];
    strMsg(ismember(strMsg,' ')) = [];

    updateResults(handles, strMsg);
end

% ******************************************************************************
%
% ******************************************************************************
function pushbuttonGenerateKML_Callback(hObject, eventdata, handles)

    global gTrackLineColor;
    
    diveNumbers = get(handles.listboxDiveNumbers, 'String');
    selectedDiveNumbers = get(handles.listboxDiveNumbers, 'Value');
    diveNumberList = str2num(char(diveNumbers(selectedDiveNumbers)));

    baseDir = fileparts(mfilename('fullpath'));
    pathname = get(handles.editPathname, 'string');

    if (exist(fullfile(baseDir,'GenKML.mat')))
        load(fullfile(baseDir,'GenKML.mat'));
    end
    tempSeagliderID = str2num(get(handles.editSeagliderID, 'String'));

    % Value between 1 and 8, see createKML.m
    iconColor = gTrackLineColor;

    diveInfoList = [];

    h = waitbar(0, 'Extracting GPS data from all selected dives...');
    for i=1:length(diveNumberList)
        
        waitbar(i/length(diveNumberList));
        
        % look in each selected dive's .nc file for the gps info.
        try
            id_dive = sprintf('%03d%04d',tempSeagliderID,diveNumberList(i));
            
            [loginfo,eng,results] = dd_get_dive_data(id_dive,1,pathname);
            %sg_cal = results.sg_calib_constants;
           
            % Save info about each dive for use in creating the .kml file.
            diveInfoList(i).diveNum = diveNumberList(i);
           
            % where and when vehicle initially surfaced before this dive
            % GPS1 - position where glider initially surfaced at end of
            % previous dive.  Surface drift occurs from here to GPS2 below.
            diveInfoList(i).prior_surface.time = loginfo.GPS1_t;
            diveInfoList(i).prior_surface.pos.valid = loginfo.GPS1_valid;
            if (loginfo.GPS1_valid)
                diveInfoList(i).prior_surface.pos.lat = loginfo.GPS1_lat;
                diveInfoList(i).prior_surface.pos.lon = loginfo.GPS1_lon;
            end
            
            % where and when vehicle left the surface on this dive
            % GPS2 - position where glider left surface on this dive
            diveInfoList(i).dove.time = loginfo.GPS2_t;
            diveInfoList(i).dove.pos.valid = loginfo.GPS2_valid;
            if (loginfo.GPS2_valid)
                diveInfoList(i).dove.pos.lat = loginfo.GPS2_lat;
                diveInfoList(i).dove.pos.lon = loginfo.GPS2_lon;
            end
           
            % where and when vehicle surfaced after the dive
            % GPS - position where glider surfaced after this dive
            diveInfoList(i).surfaced.time = loginfo.GPS_t;
            diveInfoList(i).surfaced.pos.valid = loginfo.GPS_valid;
            if (loginfo.GPS_valid)
                diveInfoList(i).surfaced.pos.lat = loginfo.GPS_lat;
                diveInfoList(i).surfaced.pos.lon = loginfo.GPS_lon;
            end
            
            % where vehicle was headed on the dive
            diveInfoList(i).target.name = loginfo.TGT_NAME;
            diveInfoList(i).target.location = loginfo.TGT_LATLONG;
            diveInfoList(i).target.radius = loginfo.TGT_RADIUS; % meters
            
            % errors that occurred during the dive
            diveInfoList(i).errors = loginfo.ERRORS;
            
        catch
            vizMsgBox('Unable to load Seaglider data file', 'ERROR: Data read', 'modal');
            closeOrphanedWaitbars();
            return;
        end

    end
    
    % Use same filename each time (do not differentiate based on start/end
    % dive number).  Google Earth will automatically update when DiveData user
    % requests to show the track after data from additional dives has arrived.
    %filename = fullfile(pathname, sprintf('p%03d%04d-%04d.kml', tempSeagliderID, diveNumberList(1), diveNumberList(end)));
    filename = fullfile(pathname, sprintf('sg%03d_track.kml', tempSeagliderID));
    name = sprintf('SG%03d', tempSeagliderID);

    createKML(handles, filename, name, tempSeagliderID, diveInfoList, iconColor);

    close(h);

    % jfaust - instead of prompting user each time, just check if
    % GoogleEarth is installed.  If so, invoke it without prompting...
    %szResponse = questdlg(sprintf('KML file generation is complete.\n\nThe file is located in\n%s.\n\nDo you want to view the file in Google Earth?', filename), 'Open KML', 'Yes', 'No', 'Yes');
    szResponse = 'Yes';
    
    % Put a msg out into the msg region saying where the .kml output file
    % was created.
    strResults = sprintf('SG%03d track for dives %d:%d written to %s\n', tempSeagliderID, diveNumberList(1), diveNumberList(end), filename);
    updateResults(handles, strResults);
    
    if strcmp(szResponse, 'Yes') == 1
        path = 'C:\Program Files\Google\Google Earth\client\googleearth.exe';
        if exist(path, 'file')
            system(sprintf('"%s" "%s" &', path, filename));
        else
            path = 'C:\Program Files (x86)\Google\Google Earth\client\googleearth.exe';
            if exist(path, 'file')
                system(sprintf('"%s" "%s" &', path, filename));
            else
                vizMsgBox(sprintf('Unable to find GoogleEarth.\n\nLeaving file at\n%s\n', filename), 'ERROR', 'modal');
            end
        end

    end

end

% ******************************************************************************
%
% ******************************************************************************
function figure1_KeyPressFcn(hObject, eventdata, handles)
    global keyModifier;
    global keyPressed;

    if strcmp(eventdata.Key, 'control') == 1 
        keyModifier = eventdata.Key;
    else
        keyPressed = eventdata.Key;
    end
end

% **********************************************************************************************************************************
%
% **********************************************************************************************************************************
function figure1_KeyReleaseFcn(hObject, eventdata, handles)
    global keyModifier;
    global keyPressed;
    global bRotate_3D;
    global gSpeedOfSound;

    if (strcmp(keyModifier, 'control') == 1) && (strcmp(keyPressed, 'k') == 1)
        pushbuttonGenerateKML_Callback(hObject, eventdata, handles);
    elseif (strcmp(keyModifier, 'control') == 1) && (strcmp(keyPressed, 'r') == 1)
        initGUI(handles);
    elseif (strcmp(keyModifier, 'control') == 1) && (strcmp(keyPressed, 's') == 1)
        gSpeedOfSound = ~gSpeedOfSound;

        set(handles.checkboxSelectPlot_8, 'Enable', toggleOnOff(get(handles.checkboxSelectPlot_8, 'Enable')));
        set(handles.checkboxSelectPlot_8, 'Value', 0);
    elseif (strcmp(keyModifier, 'control') == 1) && (strcmp(keyPressed, '3') == 1)
        bRotate_3D = ~bRotate_3D;

        strLabel = get(handles.checkboxSelectPlot_7, 'String');
        if isempty(strfind(strLabel,'</sub>'))
            strLabel = strrep(strLabel, '</html>', '<sub>3D</sub></html>');
        else
            strLabel = strrep(strLabel, '<sub>3D</sub>', '</html>');
        end
        set(handles.checkboxSelectPlot_7, 'String', strLabel);
    end
    
    keyModifier = [];
    keyPressed  = [];

    drawnow();
end

% **********************************************************************************************************************************
%
% **********************************************************************************************************************************
function [strToggledState] = toggleOnOff(strState)
    if strcmpi(strState, 'On') == 1
        strToggledState = 'Off';
    else
        strToggledState = 'On';
    end        
end

% **********************************************************************************************************************************
%
% **********************************************************************************************************************************
function checkboxUniqueFileName_Callback(hObject, eventdata, handles)
    %IMPORTFILE(FILETOREAD1)
    %  Imports data from the specified file
    %  FILETOREAD1:  file to read

    %  Auto-generated by MATLAB on 14-Sep-2012 17:08:46

    % Import the file
%     fileToRead1 = 'C:\Seahawks.png';
%     newData1 = importdata(fileToRead1);
% 
%     % Create new variables in the base workspace from those fields.
%     vars = fieldnames(newData1);
%     for i = 1:length(vars)
%         assignin('base', vars{i}, newData1.(vars{i}));
%     end

    if exist('C:\Seahawks.png', 'file')
        hFig = figure('toolbar', 'none', 'menubar', 'none');
        set(hFig, 'Name', 'Go Seahawks!');
        tmp = imread('C:\Seahawks.png');
        imshow(tmp);
        uiwait(hFig);
    end
end

% **********************************************************************************************************************************
%
% **********************************************************************************************************************************
function dead_pushbuttonFlightViz_Callback(hObject, eventdata, handles)
	vizMsgBox(sprintf('Sorry, not yet supported.'), 'ERROR', 'modal');
    % global gGliderSimState;
    % global gFlightVizHandle;
    % global gPlaybackSpeed;
    % 
    % if get(hObject, 'Value') == 1
    %     gGliderSimState = 1;
    %     gPlaybackSpeed = 0.125;
    % 
    %     diveNumbers = get(handles.listboxDiveNumbers, 'String');
    %     selectedDiveNumbers = get(handles.listboxDiveNumbers, 'Value');
    %     diveNumberList = str2num(char(diveNumbers(selectedDiveNumbers)));
    % 
    %     strDiveDataFolder = get(handles.editPathname, 'string');
    % 
    %     for i=1:length(diveNumberList)
    %         % jf - getting calib data from .nc files now.
    %         %[calibData] = readCalibConstants(fullfile(strDiveDataFolder , 'sg_calib_constants.m'));
    %         %[filenames] = buildFilenames(strDiveDataFolder, str2num(calibData.id_str), diveNumberList(i));
    % 
    %         [eng] = read_eng(filenames.p_eng);
    % 
    %         gFlightVizHandle = gliderSim(eng);
    %     end
    % else
    %     gGliderSimState = 0;
    %     if ~isempty(gFlightVizHandle)
    %         close(gFlightVizHandle);
    %     end
    % end
end

% **********************************************************************************************************************************
%
% **********************************************************************************************************************************
function dead_pushbuttonSlower_Callback(hObject, eventdata, handles)
    global gPlaybackSpeed;

    gPlaybackSpeed = gPlaybackSpeed + 0.025;
    
    if (gPlaybackSpeed > 0.5)
        gPlaybackSpeed = 0.5;
    end
end

% **********************************************************************************************************************************
%
% **********************************************************************************************************************************
function dead_pushbuttonFaster_Callback(hObject, eventdata, handles)
    global gPlaybackSpeed;

    gPlaybackSpeed = gPlaybackSpeed - 0.025;
    
    if (gPlaybackSpeed < 0)
        gPlaybackSpeed = 0;
    end
end

function checkboxAnnotatePlots_Callback(hObject, eventdata, handles)
	vizMsgBox(sprintf('Sorry, annotate plots not yet supported.'), 'ERROR', 'modal');
end

% **********************************************************************************************************************************
%
% **********************************************************************************************************************************
function callbackCalculateVelocity(~, eventdata, handles)
    %[results] = velocityCalculator();
    
	vizMsgBox(sprintf('Sorry, not yet supported.'), 'ERROR', 'modal');
    [results] = -1;
end

% **********************************************************************************************************************************
%
% **********************************************************************************************************************************
function initGUI(handles)

if 0
    hTemp = findall(gcf, 'tag', 'figMenuPilotingTools');

    if isempty(hTemp)
        hTemp = uimenu('Label', 'Piloting Tools', 'tag', 'figMenuPilotingTools');
        uimenu(hTemp, 'Label', 'Calculate Velocity', 'tag', 'figMenuPilotingToolsCalcVel', 'Callback', @callbackCalculateVelocity);
        %uimenu(hTemp, 'Label', 'Generate KML', 'tag', 'figMenuPilotingToolsGenKML', 'Callback', {@pushbuttonGenerateKML_Callback, handles});
%     else
%         uimenu(hTemp, 'Callback', @callbackCalculateVelocity);
%         uimenu(hTemp, 'Callback', @pushbuttonGenerateKML_Callback);
    end
end
end

% **********************************************************************************************************************************
%  Determine Seaglider ID 
% **********************************************************************************************************************************
function id = determineSeagliderId(dirPath)
	% Check to see if directory 'dirPath' contains any netCDF files.
	% If so, then the filename contains the seaglider ID.
	star_dot_nc = fullfile(dirPath,'p*.nc');
	nc_files_list = dir(star_dot_nc);
	if (length(nc_files_list) == 0)
        [id] = -1;
    else
        % Extract SeaGlider ID from the first .nc file name.
        % Format is pIIIDDDD.nc:
        %   III is the 3 digit seaglider id
        %   DDDD is the four digit dive number
        tmp_file = nc_files_list(1);
        tmp_file_name = tmp_file.name;
        tmp_idStr = tmp_file_name(2:4);
        [id] = str2num(tmp_idStr);
    end
end

function textResults_Callback(hObject, eventdata, handles)
% hObject    handle to textResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textResults as text
%        str2double(get(hObject,'String')) returns contents of textResults as a double
end


% --- Executes on button press in checkboxRegression_VBD_Therm.
%function checkboxRegression_VBD_Therm_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxRegression_VBD_Therm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of
% checkboxRegression_VBD_Therm
%    if (get(handles.checkboxRegression_VBD_Therm, 'Value') == 1)
%        if (get(handles.checkboxRegression_VBD_ABC, 'Value') == 1)
%            set(handles.checkboxRegression_VBD_ABC, 'Value', 0);
%        end
%    end

%end


function pushbuttonRegressDives_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRegressDives (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currentlySelectedDives = str2num(char(get(handles.listboxDiveNumbers, 'String')));

pathname = get(handles.editPathname, 'string');

settings.seagliderID    = str2num(get(handles.editSeagliderID, 'string')); % billr
id_dive = sprintf('%03d%04d',settings.seagliderID,currentlySelectedDives(1)); % billr: this is all to get the mission string - was it worth it? Maybe not.
[loginfo,eng,results] = dd_get_dive_data(id_dive,1,pathname); % billr 
settings.mission = results.GLOBALS.project; % billr

regress_dives = dd_regress_vbd_dives(handles, pathname, 1, settings); %billr: added settings param

% Update the listbox selection to match the regress_dives list
diveNumbers = [get(handles.listboxDiveNumbers, 'String')];

N = length(regress_dives);
diveNumberSelectIdx = zeros(1,N);

dives = length(diveNumbers);
diveNumbers = str2num(char(diveNumbers));

for n=1:N
    for m=1:dives
        if regress_dives(n) == diveNumbers(m)
            diveNumberSelectIdx(n) = m;
        end
    end
end


set(handles.listboxDiveNumbers, 'Value', diveNumberSelectIdx);
end
