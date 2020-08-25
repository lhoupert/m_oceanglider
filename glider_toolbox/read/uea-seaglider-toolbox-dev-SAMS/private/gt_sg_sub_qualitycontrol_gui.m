function data = gt_sg_sub_qualitycontrol_gui(data)
tmp = 0
tmp1 = 0
tmp2 = 0
close all

dives = intersect([data.log.dive],[data.eng.dive]);
if isfield(data.gt_sg_settings,'QC_settings')
    QC = data.gt_sg_settings.QC_settings;
end

variableList = fieldnames(data.hydrography);
variableList = setxor({'time','depth','used_samples','flagged_samples'},...
    variableList(...
    cellfun(@(x) ...
    all(size(data.hydrography(5).(x)) == size(data.hydrography(5).time)),...
    variableList)));

yaxisList = {'depth','time',variableList{:}};

for fstep=1:numel(variableList)
    QC.saved.(variableList{fstep}) = 1;
    blankFlags(variableList{fstep});
end

sampleNum = arrayfun(@(x) [1:numel(x.time)],data.hydrography,'Uni',0);
diveNum = arrayfun(@(x) ones(size(x.elaps_t))*x.dive,data.eng,'Uni',0);

%% Create panel
close all;
h_gui = figure('Visible','on');
% Set options
set(h_gui,...
    'MenuBar','none',...
    'Units','normalized',...
    'Position',[0.05 0.05 0.9 0.9]...
    );

%% Save buttons
h_save = uicontrol(h_gui,'Style','pushbutton','String','Save',...
    'Units','normalized','Position',[0.86 0.945 0.13 0.05],...
    'CallBack',@button_save);

uicontrol(h_gui,'Style','pushbutton','String','Save all',...
    'Units','normalized','Position',[0.86 0.89 0.13 0.05],...
    'CallBack',@button_save_all);

%% Display flagged samples?
h_showFlag = uicontrol(h_gui,'Style','checkbox','String','Show flagged samples?',...
    'Units','normalized','Position',[0.72 0.94 0.13 0.03],...
    'CallBack',@refreshPlots);

%% Variable Selection
uicontrol(h_gui,'Style','text','String','Variable',...
    'Units','normalized','Position',[0.01 0.96 0.10 0.03]);
h_variable = uicontrol(h_gui,'Style','popup','String',variableList,'Value',1,...
    'Min',1,'Max',1,'Units','normalized','Position',[0.01 0.92 0.10 0.03],...
    'CallBack',@changeVariable);

%% Profile y-axis dimension
uicontrol(h_gui,'Style','text','String','Y Axis',...
    'Units','normalized','Position',[0.22 0.96 0.10 0.03]);
h_yaxisvar = uicontrol(h_gui,'Style','popup','String',yaxisList,'Value',1,...
    'Units','normalized','Position',[0.22 0.92 0.10 0.03],...
    'CallBack',@refreshPlots);

%% Dive Selection
uicontrol(h_gui,'Style','text','String','Dive numbers',...
    'Units','normalized','Position',[0.01 0.86 0.10 0.03]);
h_list_dives = uicontrol(h_gui,'Style','listbox','String',dives,...
    'Units','normalized','Position',[0.01 0.09 0.10 0.76],...
    'Max',10000,'Min',1,'Value',dives,...
    'CallBack',@refreshPlots);
uicontrol(h_gui,'Style','pushbutton','String','Select all',...
    'Units','normalized','Position',[0.01 0.05 0.10 0.03],...
    'CallBack',@(x,y) cellfun(@(z)feval(z),...
    {@(x,y) set(h_list_dives,'Value',1:numel(dives)),...
    @refreshPlots}));
uicontrol(h_gui,'Style','pushbutton','String','Select none',...
    'Units','normalized','Position',[0.01 0.01 0.10 0.03],...
    'CallBack',@(x,y) cellfun(@(z)feval(z),...
    {@(x,y) set(h_list_dives,'Value',[]),...
    @refreshPlots}));

%% Cut off thresholds
uicontrol(h_gui,'Style','text','String','Minimum',...
    'Units','normalized','Position',[0.42 0.96 0.10 0.03]);
h_min = uicontrol(h_gui,'Style','edit','String',[],...
    'Units','normalized','Position',[0.42 0.92 0.10 0.03],...
    'CallBack',@refreshPlots);

uicontrol(h_gui,'Style','text','String','Maximum',...
    'Units','normalized','Position',[0.53 0.96 0.10 0.03]);
h_max = uicontrol(h_gui,'Style','edit','String',[],...
    'Units','normalized','Position',[0.53 0.92 0.10 0.03],...
    'CallBack',@refreshPlots);

%% Filters
h_tabs = uitabgroup(h_gui,'Position',[0.54 0.01 0.45 0.45]);

h_tab(1).tab = uitab(h_tabs,'Title','+',...
    'ButtonDownFcn',@createNewFilter);
h_tab(2).tab = uitab(h_tabs,'Title','Help');

filterTypes = {'Num. stddev beyond the mean','| x - mean | / stddev','| x - mean |','| x - median |'};
filterFuns = {'filtStddev','filtZscore','filtMeandiff','filtMediandiff'};

fwindowTypes = {'Sample Num.',yaxisList{:}};
fdiveTypes = {'All dives','Selected dives'};

createNewFilter(0,0);
activeFilters = [];

    function createNewFilter(~,~)
        % Move the "+" and "Help" tabs to the end
        h_tab(end+1) = h_tab(end);
        h_tab(end-1) = h_tab(end-2);
        tabNum = numel(h_tab)-2;
        
        h_tab(tabNum).tab = uitab(h_tabs,'Title',['Filt.',num2str(tabNum)],...
            'ButtonDownFcn',@changeTab);
        
        % Run Filter?
        h_tab(tabNum).check = uicontrol(h_tab(tabNum).tab,'Style','pushbutton','String','Apply all filters',...
            'Units','normalized','Position',[0.01 0.81 0.18 0.18],...
            'CallBack',@refreshPlots);
        
        % Info on which filters are active
        h_tab(tabNum).active = uicontrol(h_tab(tabNum).tab,'Style','text','String',['Active filters:'],...
            'Units','normalized','Position',[0.5 0.94 0.48 0.05],'HorizontalAlignment','left',...
            'CallBack',@refreshPlots);
        
        % Text bits
        uicontrol(h_tab(tabNum).tab,'Style','text','String','Filter type: ',...
            'Units','normalized','Position',[0.21 0.83 0.28 0.05],'HorizontalAlignment','right',...
            'CallBack',@refreshPlots);
        uicontrol(h_tab(tabNum).tab,'Style','text','String','Threshold: ',...
            'Units','normalized','Position',[0.01 0.73 0.48 0.05],'HorizontalAlignment','right',...
            'CallBack',@refreshPlots);
        uicontrol(h_tab(tabNum).tab,'Style','text','String','Window dimension: ',...
            'Units','normalized','Position',[0.01 0.63 0.48 0.05],'HorizontalAlignment','right',...
            'CallBack',@refreshPlots);
        uicontrol(h_tab(tabNum).tab,'Style','text','String','Extend window to: ',...
            'Units','normalized','Position',[0.01 0.53 0.48 0.05],'HorizontalAlignment','right',...
            'CallBack',@refreshPlots);
        uicontrol(h_tab(tabNum).tab,'Style','text','String','(Optional) Window range: ',...
            'Units','normalized','Position',[0.01 0.43 0.48 0.05],'HorizontalAlignment','right',...
            'CallBack',@refreshPlots);
        uicontrol(h_tab(tabNum).tab,'Style','text','String','(Optional) Window resolution: ',...
            'Units','normalized','Position',[0.01 0.33 0.48 0.05],'HorizontalAlignment','right',...
            'CallBack',@refreshPlots);
        
        % Filter settings
        h_tab(tabNum).filterType = uicontrol(h_tab(tabNum).tab,'Style','popup','String',filterTypes,...
            'Value',1,'Units','normalized','Position',[0.5 0.8 0.48 0.09],...
            'CallBack',@refreshPlots);
        h_tab(tabNum).threshType = uicontrol(h_tab(tabNum).tab,'Style','edit','String',[],...
            'Value',1,'Units','normalized','Position',[0.5 0.7 0.48 0.09],...
            'CallBack',@refreshPlots);
        h_tab(tabNum).fwindowType = uicontrol(h_tab(tabNum).tab,'Style','popup','String',fwindowTypes,...
            'Value',1,'Units','normalized','Position',[0.5 0.6 0.48 0.09],...
            'CallBack',@refreshPlots);
        h_tab(tabNum).fdiveType = uicontrol(h_tab(tabNum).tab,'Style','popup','String',fdiveTypes,...
            'Value',1,'Units','normalized','Position',[0.5 0.5 0.48 0.09],...
            'CallBack',@refreshPlots);
        h_tab(tabNum).frangeType = uicontrol(h_tab(tabNum).tab,'Style','edit','String',[],...
            'Value',1,'Units','normalized','Position',[0.5 0.4 0.48 0.09],...
            'CallBack',@refreshPlots);
        h_tab(tabNum).fresType = uicontrol(h_tab(tabNum).tab,'Style','edit','String',[],...
            'Value',1,'Units','normalized','Position',[0.5 0.3 0.48 0.09],...
            'CallBack',@refreshPlots);
        
        % Reorder tabs
        set(h_tabs,'Children',[h_tab.tab]','SelectedTab',h_tab(tabNum).tab);
    end

%% Profiles
h_prof = axes('YAxisLocation','left','Box','on',...
    'Units','normalized','Position',[0.16 0.06 0.35 0.82],...
    'XGrid','on','YGrid','on',...
    'YDir','reverse',...
    'Parent',h_gui);
h_prof_xmin = uicontrol(h_gui,'Style','edit','String',[],...
    'Units','normalized','Position',[0.15 0.03 0.03 0.03],...
    'CallBack',@changeAxes);
h_prof_xmax = uicontrol(h_gui,'Style','edit','String',[],...
    'Units','normalized','Position',[0.49 0.03 0.03 0.03],...
    'CallBack',@changeAxes);
h_prof_ymin = uicontrol(h_gui,'Style','edit','String',[],...
    'Units','normalized','Position',[0.13 0.865 0.03 0.03],...
    'CallBack',@changeAxes);
h_prof_ymax = uicontrol(h_gui,'Style','edit','String',[],...
    'Units','normalized','Position',[0.13 0.055 0.03 0.03],...
    'CallBack',@changeAxes);

%% Histogram
h_hist = axes('YAxisLocation','right','Box','on',...
    'Units','normalized','Position',[0.55 0.50 0.40 0.38],...
    'XGrid','on','YGrid','on',...
    'Parent',h_gui);
h_hist_xmin = uicontrol(h_gui,'Style','edit','String',[],...
    'Units','normalized','Position',[0.54 0.47 0.03 0.03],...
    'CallBack',@changeAxes);
h_hist_xmax = uicontrol(h_gui,'Style','edit','String',[],...
    'Units','normalized','Position',[0.93 0.47 0.03 0.03],...
    'CallBack',@changeAxes);

%% Initialise
selectedDives = dives;
variable = variableList{get(h_variable,'Value')};
changeVariable(0,0);

%% SUBFUNCTIONS

    function changeVariable(~,~)
        % Initialises GUI for a set variable
        if ~QC.saved.(variable)
            choice = questdlg({[variable,' flags not saved.'],'Would you like to proceed'}, ...
                'Flags not saved', ...
                'Yes','Cancel','Cancel');
            % Handle response
            if ~strcmp(choice,'Yes')
                return;
            end
        end
        
        variable = variableList{get(h_variable,'Value')};
        
        set(h_save,'String',['Save ',variable]);
        
        restoreFlagsSettings;
        
        populateAxesLim;  %TODO (set xmin/ax etc.)
        
        refreshPlots;
        
    end

    function refreshPlots(~,~)
        % Replots everything anytime there is a change.
        selectedDives = dives(h_list_dives.Value);
        
        % Mark bad values
        applyCutoffs
        applyFilters
        
        % Do plots
        doProfile;
        doHistogram;
        populateAxesLim;
        
        % Cosmetics
        trimTickLabels;
        
        QC.saved.(variable) = 0;
    end

    function applyCutoffs
        
        blankFlags(variable);
        
        tmp = arrayfun(@(x,y) ~isfinite(x.(variable))+y.(variable),...
            data.hydrography(selectedDives),QC.flag(selectedDives),'Uni',0);
        [QC.flag(selectedDives).(variable)] = tmp{:};
        
        if ~isempty(h_min.String)
            tmp = arrayfun(@(x,y) (x.(variable) < str2double(h_min.String))+y.(variable),...
                data.hydrography(selectedDives),QC.flag(selectedDives),'Uni',0);
            [QC.flag(selectedDives).(variable)] = tmp{:};
        end
        if ~isempty(h_max.String)
            tmp = arrayfun(@(x,y) (x.(variable) > str2double(h_max.String))+y.(variable),...
                data.hydrography(selectedDives),QC.flag(selectedDives),'Uni',0);
            [QC.flag(selectedDives).(variable)] = tmp{:};
        end
        
    end

    function applyFilters
        
        % Now start applying filters
        numFilters = numel(h_tab) - 2;
        activeFilters = find(arrayfun(@(x) x.check.Value,h_tab(1:numFilters)));
        
        if isempty(activeFilters)
            return;
        end
        disp('Starting filters... please wait.');
        for tstep = activeFilters            
            % Extract parameters
            range = str2double(h_tab(tstep).frangeType.String);
            
            threshold = str2double(h_tab(tstep).threshType.String);
            
            fdives = subsref({[1:numel(data.hydrography)],selectedDives},...
                struct('type','{}','subs',{{h_tab(tstep).fdiveType.Value}}));
            
            dim = h_tab(tstep).fwindowType.String{h_tab(tstep).fwindowType.Value};
            if strcmp(dim,'Sample Num.')
                fdimArray = indexarray([sampleNum{fdives}],~[QC.flag(fdives).(variable)]);
            else
                fdimArray = indexarray([data.hydrography(fdives).(dim)],~[QC.flag(fdives).(variable)]);
            end
            
            if isempty(range) || isnan(range)
                range = Inf;
            end
            
            dataArray = indexarray([data.hydrography(fdives).(variable)],~[QC.flag(fdives).(variable)]);
            
            % Define filter to run
            switch filterFuns{h_tab(tstep).filterType.Value}
                case 'filtStddev'
                    filtStddev;
                case 'filtZscore'
                    filtZscore;
                case 'filtMeandiff'
                    filtMeandiff;
                case 'filtMediandiff'
                    filtMediandiff;
            end
        end
        disp('Finished filtering');
        
        % Define filter functions
        function redistributeFlags(flagarray)
            diveInd = indexarray([diveNum{fdives}],~[QC.flag(fdives).(variable)]);
            sampInd = indexarray([sampleNum{fdives}],~[QC.flag(fdives).(variable)]);
            
            indices = find(flagarray);
            
            for istep = indices
                QC.flag(diveInd(istep)).(variable)(sampInd(istep)) = 1;
            end
        end
        function filtStddev %TODO write me
            tmp = arrayfun(@(x,y) abs(x - mean(dataArray(...
                y+range > fdimArray & y-range < fdimArray)))...
                > ...
                threshold*sqrt(var(dataArray(...
                y+range > fdimArray & y-range < fdimArray))),dataArray,fdimArray);
            redistributeFlags(tmp)
        end
        function filtZscore %TODO write me
            
        end
        function filtMeandiff %TODO write me
            
        end
        function filtMediandiff %TODO write me
            
        end
    end

    function button_save(~,~)
        saveFlags(variable);
        saveSettings(variable);
        gt_sg_sub_echo({['Saving flag data for ',variable,'.']});
        QC.saved.(variable) = 1;
    end

    function button_save_all(~,~)
        for fstep = 1:numel(variableList)
            saveFlags(variableList{fstep});
            saveSettings(variableList{fstep});
            QC.saved.(variableList{fstep}) = 1;
        end
        gt_sg_sub_echo({['Saving flag data for all variables.']});
    end

%% SUBSUBFUNCTIONS

    function blankFlags(var)
        tmp = arrayfun(@(x) zeros(size(x.(var))),data.hydrography,'Uni',0);
        [QC.flag(1:numel(data.hydrography),1).(var)] = tmp{:};
    end

    function saveFlags(var)
        [data.flag(1:numel(data.hydrography),1).(var)] = QC.flag(1:numel(data.hydrography)).(var);
    end

    function saveSettings(var)
        data.gt_sg_settings.QC.settings.(var).min = str2double(h_min.String);
        data.gt_sg_settings.QC.settings.(var).max = str2double(h_max.String);
        QC.prof_XLim = h_prof.XLim;
        QC.prof_YLim = h_prof.YLim;
        QC.hist_XLim = h_hist.XLim;
    end

    function restoreFlagsSettings
        set(h_min,'String',[])
        set(h_max,'String',[])
        
        % Restore saved values if exist
        try
            [QC.flag(1:numel(data.hydrography),1).(variable)] = data.flag(1:numel(data.hydrography)).(variable);
        end
        try
            set(h_min,'String',data.gt_sg_settings.QC.settings.(variable).min)
        end
        try
            set(h_max,'String',data.gt_sg_settings.QC.settings.(variable).max)
        end
    end

    function populateAxesLim
        set(h_prof_xmin,'String',h_prof.XLim(1));
        set(h_prof_xmax,'String',h_prof.XLim(2));
        
        set(h_prof_ymin,'String',h_prof.YLim(1));
        set(h_prof_ymax,'String',h_prof.YLim(2));
        
        
        set(h_hist_xmin,'String',h_hist.XLim(1));
        set(h_hist_xmax,'String',h_hist.XLim(2));
    end

    function changeAxes(~,~)
        set(h_prof,'XLim',[str2double(h_prof_xmin.String) str2double(h_prof_xmax.String)],...
            'YLim',[str2double(h_prof_ymin.String) str2double(h_prof_ymax.String)]);
        
        set(h_hist,'XLim',[str2double(h_hist_xmin.String) str2double(h_hist_xmax.String)]);
        
        trimTickLabels;
    end

    function changeTab(~,~)
        %restore settings
        %update list of active filters
    end

    function trimTickLabels
        set(h_prof,'XTickLabel',arrayfun(@(x) num2str(x),h_prof.XTick,'Uni',0)');
        h_prof.XTickLabel{1} = [];
        h_prof.XTickLabel{end} = [];
        
        set(h_prof,'YTickLabel',arrayfun(@(x) num2str(x),h_prof.YTick,'Uni',0)');
        h_prof.YTickLabel{1} = [];
        h_prof.YTickLabel{end} = [];
        
        set(h_hist,'XTickLabel',arrayfun(@(x) num2str(x),h_hist.XTick,'Uni',0)');
        h_hist.XTickLabel{1} = [];
        h_hist.XTickLabel{end} = [];
        
        set(h_hist,'YTickLabel',arrayfun(@(x) num2str(x),h_hist.YTick,'Uni',0)');
        h_hist.YTickLabel{1} = [];
    end

    function doHistogram
        cla(h_hist);
        if h_showFlag.Value
            dataArray = [data.hydrography(selectedDives).(variable)];
        else
            dataArray = indexarray([data.hydrography(selectedDives).(variable)],~[QC.flag(selectedDives).(variable)]);
        end
        set(h_hist,'NextPlot','add');
        hist(h_hist,dataArray,1000);
        
        tmp = nanmean(dataArray);
        tmp2 = nanstd(dataArray);
        tmp3 = h_hist.XLim;
        plot(h_hist,[tmp tmp],h_hist.YLim,':k','LineWidth',1);
        plot(h_hist,[tmp-tmp2 tmp-tmp2],h_hist.YLim,':r','LineWidth',1);
        plot(h_hist,[tmp+tmp2 tmp+tmp2],h_hist.YLim,':r','LineWidth',1);
        
        set(h_hist,'XGrid','on','YGrid','on','NextPlot','replace');
    end

    function doProfile
        if h_showFlag.Value
            plot(h_prof,...
                indexarray([data.hydrography(selectedDives).(variable)],~[QC.flag(selectedDives).(variable)]),...
                indexarray([data.hydrography(selectedDives).(yaxisList{h_yaxisvar.Value})],~[QC.flag(selectedDives).(variable)]),...
                'k','LineStyle','none','Marker','.');
            set(h_prof,'NextPlot','add');
            
            plot(h_prof,...
                indexarray([data.hydrography(selectedDives).(variable)],[QC.flag(selectedDives).(variable)]),...
                indexarray([data.hydrography(selectedDives).(yaxisList{h_yaxisvar.Value})],[QC.flag(selectedDives).(variable)]),...
                'r','LineStyle','none','Marker','.','Color',[1 0.5 0.5]);
            set(h_prof,'NextPlot','replace');
        else
            plot(h_prof,...
                indexarray([data.hydrography(selectedDives).(variable)],~[QC.flag(selectedDives).(variable)]),...
                indexarray([data.hydrography(selectedDives).(yaxisList{h_yaxisvar.Value})],~[QC.flag(selectedDives).(variable)]),...
                'k','LineStyle','none','Marker','.');
        end
        if any(strcmp(yaxisList{h_yaxisvar.Value},{'depth','pressure','sigma0','rho'}));
            set(h_prof,'YDir','reverse',...
                'XGrid','on','YGrid','on');
        else
            set(h_prof,'YDir','normal',...
                'XGrid','on','YGrid','on');
        end
    end

    function out = indexarray(array,flags)
        % I like this one :)
        out = ...
            subsref(array,...
            struct('type','()','subs', {{ find(flags) }} ));
    end

%% WAIT FOR GUI CLOSE
%uiwait(h_gui)
keyboard
end
