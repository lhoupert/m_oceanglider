function data = gt_sg_sub_flightmodelregression_gui(data)

HelpInfo = ...
    {'If you want the quick and easy option, just click the button on the right. ->',...
    ' ',...
    'Otherwise, use the first tab to select useful data and regress in the second tab. The two cost functions work better on different things. The one labelled hydro works best on hd_ parameters. The one labelled buoyancy copes best with volmax, abs_compress and therm_expan.',...
    ' ',...
    'Info on visualisation tab...'};


%% Create panel
close all;
h_gui = figure('Visible','on');
% Set options
set(h_gui,...
    'MenuBar','none',...
    'Units','normalized',...
    'Position',[0.05 0.05 0.9 0.9]...
    );

%% Extract basic variables and setup settings storage
dives = intersect([data.log.dive],[data.eng.dive]);
cast_dir = arrayfun(@(x) ((x.elaps_t*0)+1).*x.dive,data.eng,'Uni',0); % Dive no
for dstep=1:numel(cast_dir)
    cast_dir{dstep}(1:data.hydrography(dstep).max_pressure_index) = -cast_dir{dstep}(1:data.hydrography(dstep).max_pressure_index); % Sign for direction
end

if isfield(data.gt_sg_settings,'flight_regression_settings')
    FRS = data.gt_sg_settings.flight_regression_settings;
else
    % Data selection
    FRS.surface_cutoff = 15;
    FRS.apogee_cutoff = 15;
    FRS.roll_max = 4;
    FRS.pitch_max = 50;
    FRS.stall_max = 10;
    FRS.motor_trim = 1;
    FRS.dives = dives;
    FRS.use_regress = [];
    FRS.vertaccel_max = 0.01;
    FRS.vertvel_min = 0.02;
end
% Hull length
HullLength = {'Standard',1.8;...
    'Ogive',2.0;...
    'Deep (N/A)',2.2};
if isfield(data.gt_sg_settings,'length')
    hullLength = data.gt_sg_settings.length;
else
    try
        hullLength = mode([data.log.LENGTH]);
    catch
        hullLength = 1.8;
    end
end
hullLength = find([HullLength{:,2}] == hullLength);
data.gt_sg_settings.length = HullLength{hullLength,2};

% Set default constraint limits
h_regress_para.hd_a.min = 0.001; h_regress_para.hd_a.max = 0.007;
h_regress_para.hd_b.min = 0.004; h_regress_para.hd_b.max = 0.02;
h_regress_para.hd_c.min = 1.0e-6; h_regress_para.hd_c.max = 3.0e-5;
h_regress_para.volmax.min = 5.0e4; h_regress_para.volmax.max = 5.5e4;
h_regress_para.abs_compress.min = 1.0e-06; h_regress_para.abs_compress.max = 3.0e-5;
h_regress_para.therm_expan.min = 5.0e-05; h_regress_para.therm_expan.max = 5.0e-4;

%% Populate tabs
h_tabs = uitabgroup(h_gui,'Position',[0 0 1 1]);
h_tab(1) = uitab(h_tabs,'Title','1. Data Selection');
h_tab(2) = uitab(h_tabs,'Title','2. Parameter Regression');
h_tab(3) = uitab(h_tabs,'Title','3. Visualisation');
h_tab(4) = uitab(h_tabs,'Title','4. History');
h_tab(5) = uitab(h_tabs,'Title','5. Help');

%% Tab 1 - Data Selection
uicontrol(h_tab(1),'Style','pushbutton','String','Save settings',...
    'Units','normalized','Position',[0.87 0.9 0.12 0.07],...
    'CallBack',@button_data_save);
uicontrol(h_tab(1),'Style','pushbutton','String','Plot selection',...
    'Units','normalized','Position',[0.74 0.9 0.12 0.07],...
    'CallBack',@button_data_plot);

% Apogee and surface cutoff
uicontrol(h_tab(1),'Style','text','String','Surface cutoff (m)',...
    'Units','normalized','Position',[0.01 0.94 0.10 0.03]);
uicontrol(h_tab(1),'Style','text','String','Apogee cutoff (m)',...
    'Units','normalized','Position',[0.01 0.90 0.10 0.03]);
h_data_cutoff_surface = uicontrol(h_tab(1),'Style','edit','String',FRS.surface_cutoff,...
    'Units','normalized','Position',[0.12 0.94 0.06 0.03]);
h_data_cutoff_apogee = uicontrol(h_tab(1),'Style','edit','String',FRS.apogee_cutoff,...
    'Units','normalized','Position',[0.12 0.90 0.06 0.03]);

% Roll, pitch, stall and post motor movement cutoff
uicontrol(h_tab(1),'Style','text','String','Roll max (deg)',...
    'Units','normalized','Position',[0.19 0.94 0.10 0.03]);
uicontrol(h_tab(1),'Style','text','String','Pitch max (deg)',...
    'Units','normalized','Position',[0.19 0.90 0.10 0.03]);
h_data_cutoff_rollmax = uicontrol(h_tab(1),'Style','edit','String',FRS.roll_max,...
    'Units','normalized','Position',[0.30 0.94 0.06 0.03]);
h_data_cutoff_pitchmax = uicontrol(h_tab(1),'Style','edit','String',FRS.pitch_max,...
    'Units','normalized','Position',[0.30 0.90 0.06 0.03]);
uicontrol(h_tab(1),'Style','text','String','Stall angle (deg)',...
    'Units','normalized','Position',[0.37 0.94 0.10 0.03]);
uicontrol(h_tab(1),'Style','text','String','Motor trim (samples)',...
    'Units','normalized','Position',[0.37 0.90 0.10 0.03]);
h_data_cutoff_stallmax = uicontrol(h_tab(1),'Style','edit','String',FRS.stall_max,...
    'Units','normalized','Position',[0.48 0.94 0.06 0.03]);
h_data_cutoff_motor = uicontrol(h_tab(1),'Style','edit','String',FRS.motor_trim,...
    'Units','normalized','Position',[0.48 0.90 0.06 0.03]);
uicontrol(h_tab(1),'Style','text','String','Min dP/dt (m/s)',...
    'Units','normalized','Position',[0.55 0.94 0.10 0.03]);
uicontrol(h_tab(1),'Style','text','String','Max accel. (m/s^2)',...
    'Units','normalized','Position',[0.55 0.90 0.10 0.03]);
h_data_cutoff_vertvelmin = uicontrol(h_tab(1),'Style','edit','String',FRS.vertvel_min,...
    'Units','normalized','Position',[0.66 0.94 0.06 0.03]);
h_data_cutoff_vertaccelmax = uicontrol(h_tab(1),'Style','edit','String',FRS.vertaccel_max,...
    'Units','normalized','Position',[0.66 0.90 0.06 0.03]);

% Dive Selection
uicontrol(h_tab(1),'Style','text','String','Dive numbers',...
    'Units','normalized','Position',[0.01 0.86 0.10 0.03]);
h_data_list_dives = uicontrol(h_tab(1),'Style','listbox','String',dives,...
    'Units','normalized','Position',[0.12 0.01 0.06 0.88],...
    'Max',10000,'Min',1,'Value',indicesIn(dives,FRS.dives));
uicontrol(h_tab(1),'Style','pushbutton','String','Select all',...
    'Units','normalized','Position',[0.01 0.05 0.10 0.03],...
    'CallBack',@(x,y) set(h_data_list_dives,'Value',1:numel(dives)));
uicontrol(h_tab(1),'Style','pushbutton','String','Select none',...
    'Units','normalized','Position',[0.01 0.01 0.10 0.03],...
    'CallBack',@(x,y) set(h_data_list_dives,'Value',[]));

% Plots
h_data_plot_pitch = axes('YAxisLocation','right','Box','on',...
    'Units','normalized','Position',[0.2 0.035 0.74 0.14],'Parent',h_tab(1));
set(h_data_plot_pitch.YLabel,'String','Pitch angle (deg)')
hold on;
plot(h_data_plot_pitch,[data.hydrography(dives).time],[data.eng(dives).pitchAng],'LineStyle','none','Marker','.','MarkerSize',1,'Color',[0.8 0.8 0.8]);
h_data_plot_pitch_blue = plot(h_data_plot_pitch,[],[],'LineStyle','none','Marker','.','MarkerSize',1,'Color','b'); datetick(h_data_plot_pitch,'x'); axis(h_data_plot_pitch,'tight');

h_data_plot_roll = axes('YAxisLocation','right','Box','on',...
    'Units','normalized','Position',[0.2 0.21 0.74 0.14],'Parent',h_tab(1));
set(h_data_plot_roll.YLabel,'String','Roll angle (deg)')
hold on;
plot(h_data_plot_roll,[data.hydrography(dives).time],[data.eng(dives).rollAng],'LineStyle','none','Marker','.','MarkerSize',1,'Color',[0.8 0.8 0.8]);
h_data_plot_roll_blue = plot(h_data_plot_roll,[],[],'LineStyle','none','Marker','.','MarkerSize',1,'Color','b'); datetick(h_data_plot_roll,'x'); axis(h_data_plot_roll,'tight');

h_data_plot_total = axes('YAxisLocation','right','Box','on',...
    'Units','normalized','Position',[0.2 0.385 0.74 0.14],'Parent',h_tab(1));
set(h_data_plot_total.YLabel,'String','dP/dt (m/s)')
hold on;
plot(h_data_plot_total,[data.hydrography(dives).time],[data.flight(dives).glide_vert_spd],'LineStyle','none','Marker','.','MarkerSize',1,'Color',[0.8 0.8 0.8]);
h_data_plot_total_blue = plot(h_data_plot_total,[],[],'LineStyle','none','Marker','.','MarkerSize',1,'Color','b'); datetick(h_data_plot_total,'x'); axis(h_data_plot_total,'tight');

h_data_plot_accel = axes('YAxisLocation','right','Box','on',...
    'Units','normalized','Position',[0.2 0.56 0.74 0.14],'Parent',h_tab(1));
set(h_data_plot_accel.YLabel,'String','Acceleration (m/s^2)')
hold on;
plot(h_data_plot_accel,[data.hydrography(dives).time],[data.flight(dives).glide_vert_spd],'LineStyle','none','Marker','.','MarkerSize',1,'Color',[0.8 0.8 0.8]);
h_data_plot_accel_blue = plot(h_data_plot_total,[],[],'LineStyle','none','Marker','.','MarkerSize',1,'Color','b'); datetick(h_data_plot_accel,'x'); axis(h_data_plot_accel,'tight');

h_data_plot_vbd = axes('YAxisLocation','right','Box','on',...
    'Units','normalized','Position',[0.2 0.735 0.74 0.14],'Parent',h_tab(1));
set(h_data_plot_vbd.YLabel,'String','C\_VBD')
hold on;
plot(h_data_plot_vbd,dives,[data.log(dives).C_VBD],'LineStyle','none','Marker','*','Color',[0.8 0.8 0.8]);
h_data_plot_vbd_blue = plot(h_data_plot_vbd,[],[],'LineStyle','none','Marker','*','Color','b'); axis(h_data_plot_vbd,'tight');

% Populate with current selection
func_data_makeplots;

% Callbacks
    function button_data_save(button_data,event_data)
        button_data_plot(button_data,event_data)
        data.gt_sg_settings.flight_regression_settings = FRS;
        set(h_tabs,'SelectedTab',h_tab(2));
        
        % Initialise second tab
        button_regress_regresscheckbox(0,0);
        button_regress_plot_velocities(0,0);
    end

    function button_data_plot(~,~)
        FRS.surface_cutoff = str2double(get(h_data_cutoff_surface,'String'));
        FRS.apogee_cutoff = str2double(get(h_data_cutoff_apogee,'String'));
        FRS.roll_max = str2double(get(h_data_cutoff_rollmax ,'String'));
        FRS.pitch_max = str2double(get(h_data_cutoff_pitchmax,'String'));
        FRS.stall_max = str2double(get(h_data_cutoff_stallmax,'String'));
        FRS.motor_trim = str2double(get(h_data_cutoff_motor,'String'));
        FRS.vertvel_min = str2double(get(h_data_cutoff_vertvelmin,'String'));
        FRS.vertaccel_max = str2double(get(h_data_cutoff_vertaccelmax,'String'));
        FRS.dives = dives(get(h_data_list_dives,'Value'));
        
        func_data_makeplots;
    end

    function func_data_makeplots
        delete(h_data_plot_total_blue);
        delete(h_data_plot_accel_blue);
        delete(h_data_plot_vbd_blue);
        delete(h_data_plot_pitch_blue);
        delete(h_data_plot_roll_blue);
        
        for dstep = 1:numel(data.hydrography);
            FRS.use_regress{dstep} = false(size(data.hydrography(dstep).depth));
            if any(FRS.dives == dstep)
                GC_phase = data.eng(dstep).GC_phase;
                if FRS.motor_trim < 0
                    GC_phase(:) = 6;
                    motor_trim = 1;
                else
                    motor_trim = FRS.motor_trim;
                end
                FRS.use_regress{dstep}(...
                    data.hydrography(dstep).depth > FRS.surface_cutoff ...
                    & data.hydrography(dstep).depth < data.hydrography(dstep).depth(data.hydrography(dstep).max_pressure_index)-FRS.apogee_cutoff ...
                    & abs(data.eng(dstep).rollAng) < FRS.roll_max ...
                    & abs(data.eng(dstep).pitchAng) < FRS.pitch_max ...
                    & abs(data.eng(dstep).pitchAng) > FRS.stall_max ...
                    & abs(data.flight(dstep).glide_vert_spd) > FRS.vertvel_min ...
                    & abs(gt_sg_sub_diff({data.flight(dstep).glide_vert_spd,data.hydrography(dstep).time.*24.*60.*60})) < FRS.vertaccel_max ...
                    & convn(GC_phase == 6,ones(1,motor_trim*2+1)/(motor_trim*2+1),'same') == 1 ...
                    ) = true;
                
            end
        end
        
        ind = [FRS.use_regress{dives}];
        time = [data.hydrography(dives).time];
        vert_spd = [data.flight(dives).glide_vert_spd];
        pitchAng = [data.eng(dives).pitchAng];
        rollAng = [data.eng(dives).rollAng];
        accel = gt_sg_sub_diff({[data.flight.glide_vert_spd],[data.hydrography.time]*24*60*60});
        
        uicontrol(h_tab(1),'Style','text','String',{'Number of Points:',num2str(sum(ind)),'of',num2str(numel(time))},...
            'Units','normalized','Position',[0.01 0.4 0.10 0.1],...
            'CallBack',@(x,y) set(h_data_list_dives,'Value',1:numel(dives)));
        
        h_data_plot_pitch_blue = plot(h_data_plot_pitch,time(ind),pitchAng(ind),'LineStyle','none','Marker','.','MarkerSize',1,'Color','b'); datetick(h_data_plot_pitch,'x'); axis(h_data_plot_pitch,'tight');
        h_data_plot_roll_blue = plot(h_data_plot_roll,time(ind),rollAng(ind),'LineStyle','none','Marker','.','MarkerSize',1,'Color','b'); datetick(h_data_plot_roll,'x'); axis(h_data_plot_roll,'tight');
        h_data_plot_total_blue = plot(h_data_plot_total,time(ind),vert_spd(ind),'LineStyle','none','Marker','.','MarkerSize',1,'Color','b'); datetick(h_data_plot_total,'x'); axis(h_data_plot_total,'tight');
        h_data_plot_accel_blue = plot(h_data_plot_accel,time(ind),accel(ind),'LineStyle','none','Marker','.','MarkerSize',1,'Color','b'); datetick(h_data_plot_total,'x'); axis(h_data_plot_total,'tight');
        h_data_plot_vbd_blue = plot(h_data_plot_vbd,FRS.dives,[data.log(FRS.dives).C_VBD],'LineStyle','none','Marker','*','Color','b'); axis(h_data_plot_vbd,'tight');
    end

%% Tab 2 - Parameter Regression
uicontrol(h_tab(2),'Style','pushbutton','String','View results',...
    'Units','normalized','Position',[0.87 0.91 0.12 0.07],...
    'CallBack',@button_regress_plot_velocities);
uicontrol(h_tab(2),'Style','pushbutton','String','Store results',...
    'Units','normalized','Position',[0.87 0.83 0.12 0.07],...
    'CallBack',@button_regress_store);
uicontrol(h_tab(2),'Style','pushbutton','String','Undo',...
    'Units','normalized','Position',[0.87 0.75 0.12 0.07],...
    'CallBack',@button_regress_undo);
uicontrol(h_tab(2),'Style','pushbutton','String','Restore defaults',...
    'Units','normalized','Position',[0.87 0.67 0.12 0.07],...
    'CallBack',@button_regress_reset);

% Panels
h_regress_buoy = uipanel('Parent',h_tab(2),'Title','Buoyancy parameters',...
    'Units','normalized','Position',[0.01 0.01 0.485 .49]);
h_regress_hydro = uipanel('Parent',h_tab(2),'Title','Hydrodynamic parameters',...
    'Units','normalized','Position',[0.01 0.505 0.485 .485]);
h_regress_cost = uibuttongroup('Parent',h_tab(2),'Title','Cost function',...
    'Units','normalized','Position',[0.50 0.505 0.18 .485]);
h_regress_opt = uipanel('Parent',h_tab(2),'Title','Settings',...
    'Units','normalized','Position',[0.685 0.505 0.18 .485]);
h_regress_iter = axes('YAxisLocation','right','Box','on',...
    'Units','normalized','Position',[0.50 0.1 0.21 0.39],'Parent',h_tab(2)); hold on;
h_regress_vis = axes('YAxisLocation','right','Box','on',...
    'Units','normalized','Position',[0.76 0.1 0.21 0.39],'Parent',h_tab(2)); ...
    hold on; ylabel(h_regress_vis,'Depth (m)'); xlabel(h_regress_vis,{'Red - dP/dt, Blue - buoyancy,','Green - up/downwelling'});

% Settings
uicontrol(h_regress_opt,'Style','text','String','Fairing type:',...
    'Units','normalized','Position',[0.01 0.93 0.48 0.04]);
h_regress_opt_length = uicontrol(h_regress_opt,'Style','popup','String',HullLength(:,1),'Value',hullLength,...
    'Min',1,'Max',1,'Units','normalized','Position',[0.51 0.93 0.48 0.04],...
    'CallBack',...
    @(x,y) eval(['data.gt_sg_settings.length = HullLength{get(h_regress_opt_length,''Value''),2}'])...
    );
uicontrol(h_regress_opt,'Style','text','String','Constrain:',...
    'Units','normalized','Position',[0.01 0.82 0.48 0.04]); constrainRegression = {'Off','On'};
h_regress_opt_constrain = uicontrol(h_regress_opt,'Style','popup','String',constrainRegression,'Value',1,...
    'Min',1,'Max',1,'Units','normalized','Position',[0.51 0.82 0.48 0.04],...
    'CallBack',@button_regress_constrain);
uicontrol(h_regress_opt,'Style','text','String','Max iterations: ',...
    'Units','normalized','Position',[0.01 0.71 0.48 0.04]);
h_regress_opt_maxiter = uicontrol(h_regress_opt,'Style','edit','String',100,...
    'Units','normalized','Position',[0.50 0.68 0.48 0.1]);
uicontrol(h_regress_opt,'Style','text','String','Function tolerance:',...
    'Units','normalized','Position',[0.01 0.60 0.48 0.04]);
h_regress_opt_tolfun = uicontrol(h_regress_opt,'Style','edit','String',1.0e-8,...
    'Units','normalized','Position',[0.50 0.57 0.48 0.1]);
h_regress_button_process = uicontrol(h_regress_opt,'Style','pushbutton','String','Regress parameters',...
    'Units','normalized','Position',[0.17 0.41 0.74 0.14],...
    'CallBack',@button_regress_process);
h_regress_button_auto = uicontrol(h_regress_opt,'Style','pushbutton','String','Attempt auto regress',...
    'Units','normalized','Position',[0.17 0.13 0.74 0.14],...
    'CallBack',@button_regress_auto,'Enable','on');
uicontrol(h_regress_opt,'Style','text','String',{'Number of','rangefinding passes:'},...
    'Units','normalized','Position',[0.01 0.01 0.48 0.1]);
h_regress_opt_autopass = uicontrol(h_regress_opt,'Style','edit','String','3',...
    'Units','normalized','Position',[0.50 0.01 0.48 0.1],'Enable','on');

% Cost function options
h_regress_cost_radio(1) = uicontrol(h_regress_cost,'Style','radiobutton','String',...
    'RMSD of up & down W(H2O) (hydro)',...
    'Units','normalized','Position',[0.01 0.9 0.98 0.1]);
h_regress_cost_radio(2) = uicontrol(h_regress_cost,'Style','radiobutton','String',...
    'Minimise W(H2O) (buoyancy)',...
    'Units','normalized','Position',[0.01 0.7 0.98 0.1]);
% h_regress_cost_radio(3) = uicontrol(h_regress_cost,'Style','radiobutton','String',...
%     '|W(H_2O)_u_p| + |W(H_2O)_d_o_w_n| -> 0 (coarse - all)',...
%     'Units','normalized','Position',[0.01 0.5 0.98 0.1],'Enable','on');
% h_regress_cost_radio(4) = uicontrol(h_regress_cost,'Style','radiobutton','String',...
%     'Something funky ?',...
%     'Units','normalized','Position',[0.01 0.3 0.98 0.1],'Enable','off');
% h_regress_cost_radio(5) = uicontrol(h_regress_cost,'Style','radiobutton','String',...
%     'Cost fun 1 * cost fun 2 ?? Something nice and overarching?',...
%     'Units','normalized','Position',[0.01 0.1 0.98 0.1],'Enable','off');

% Hydro and buoy panels
h_regress_para.hd_a.panel = uipanel('Parent',h_regress_hydro,'Title','HD_A',...
    'Units','normalized','Position',[0.01 0.04+2*(1-0.06)/3 0.98 (1-0.06)/3]);
h_regress_para.hd_b.panel = uipanel('Parent',h_regress_hydro,'Title','HD_B',...
    'Units','normalized','Position',[0.01 0.03+(1-0.06)/3 0.98 (1-0.06)/3]);
h_regress_para.hd_c.panel = uipanel('Parent',h_regress_hydro,'Title','HD_C',...
    'Units','normalized','Position',[0.01 0.02 0.98 (1-0.06)/3]);
h_regress_para.volmax.panel = uipanel('Parent',h_regress_buoy,'Title','Volmax',...
    'Units','normalized','Position',[0.01 0.04+2*(1-0.06)/3 0.98 (1-0.06)/3]);
h_regress_para.abs_compress.panel = uipanel('Parent',h_regress_buoy,'Title','Absolute compression',...
    'Units','normalized','Position',[0.01 0.03+(1-0.06)/3 0.98 (1-0.06)/3]);
h_regress_para.therm_expan.panel = uipanel('Parent',h_regress_buoy,'Title','Thermal expansion',...
    'Units','normalized','Position',[0.01 0.02 0.98 (1-0.06)/3]);

% Populate parameter panels
parameters = {'hd_a','hd_b','hd_c','volmax','abs_compress','therm_expan'};
for pstep = 1:numel(parameters)
    
    h_regress_para.(parameters{pstep}).regress = uicontrol(h_regress_para.(parameters{pstep}).panel,'Style','checkbox','String','Regress?',...
        'Units','normalized','Position',[0.01 0.8 0.28 0.2],...
        'CallBack',@button_regress_regresscheckbox);
    
    h_regress_para.(parameters{pstep}).initial = uicontrol(h_regress_para.(parameters{pstep}).panel,'Style','edit','String',data.gt_sg_settings.(parameters{pstep}),...
        'Enable',constrainRegression{get(h_regress_opt_constrain,'Value')}, 'Units','normalized','Position',[0.51 0.8 0.48 0.2]);
    uicontrol(h_regress_para.(parameters{pstep}).panel,'Style','text','String','Initial guess: ',...
        'Units','normalized','Position',[0.30 0.8 0.2 0.2]);
    
    h_regress_para.(parameters{pstep}).minval = uicontrol(h_regress_para.(parameters{pstep}).panel,'Style','edit','String',h_regress_para.(parameters{pstep}).min,...
        'Enable',constrainRegression{get(h_regress_opt_constrain,'Value')}, 'Units','normalized','Position',[0.41 0.05 0.23 0.2]);
    uicontrol(h_regress_para.(parameters{pstep}).panel,'Style','text','String','Minimum: ',...
        'Units','normalized','Position',[0.30 0.05 0.11 0.2]);
    
    h_regress_para.(parameters{pstep}).maxval = uicontrol(h_regress_para.(parameters{pstep}).panel,'Style','edit','String',h_regress_para.(parameters{pstep}).max,...
        'Enable',constrainRegression{get(h_regress_opt_constrain,'Value')}, 'Units','normalized','Position',[0.76 0.05 0.23 0.2]);
    uicontrol(h_regress_para.(parameters{pstep}).panel,'Style','text','String','Maximum: ',...
        'Units','normalized','Position',[0.65 0.05 0.11 0.2]);
    
    undo.(parameters{pstep}) = data.gt_sg_settings.(parameters{pstep});
end

% Callbacks
    function button_regress_process(~,~)
        
        stop_regression = false;
        set(h_regress_button_process,'String','Stop regression',...
            'CallBack',@button_regress_stop); drawnow;
        
        for pstep = 1:numel(parameters)
            undo.(parameters{pstep}) = data.gt_sg_settings.(parameters{pstep});
            settings.(parameters{pstep}) = str2double(get(h_regress_para.(parameters{pstep}).initial,'String'));
            param_min.(parameters{pstep}) = str2double(get(h_regress_para.(parameters{pstep}).minval,'String'));
            param_max.(parameters{pstep}) = str2double(get(h_regress_para.(parameters{pstep}).maxval,'String'));
        end
        % 1) make all inputs same order of magnitude
        param_coeffs = get_coeff([settings.hd_a,settings.hd_b,settings.hd_c,settings.volmax,settings.abs_compress,settings.therm_expan]);
        param_exp = get_exp([settings.hd_a,settings.hd_b,settings.hd_c,settings.volmax,settings.abs_compress,settings.therm_expan]);
        clear settings % IMPORTANT SO AS NOT TO CROSS CONTAMINATE REGRESSION SUBFUNCTION
        % 2) set options
        options = optimset('Display','iter',...
            'OutputFcn',@button_regress_plot_iteration,...
            'MaxIter',str2double(get(h_regress_opt_maxiter,'String')),...
            'TolFun',str2double(get(h_regress_opt_tolfun,'String')));
        % 3) initiate iteration plot
        delete(h_regress_iter.Children)
        set(h_regress_iter,...
            'XLim',[0 str2double(get(h_regress_opt_maxiter,'String'))],...
            'XLimMode','manual');
        h_regress_iter_plot = plot(h_regress_iter,...
            [0:str2double(get(h_regress_opt_maxiter,'String'))],...
            nan(str2double(get(h_regress_opt_maxiter,'String'))+1,1),...
            'LineStyle','none','Marker','d');
        % 4) select appropriate cost function
        switch find([h_regress_cost_radio.Value])
            case 1
                h_regress_scorefun = @button_regress_score_minuddiff;
            case 2
                h_regress_scorefun = @button_regress_score_minmeanwh2o;
            case 3
                %h_regress_scorefun = @button_regress_score_minall;
            case 4
            case 5
                %h_regress_scorefun = @
        end
        % 5) do regression AND reapply exponents
        if strcmpi(constrainRegression{get(h_regress_opt_constrain,'Value')},'on')
            for pstep = 1:numel(parameters)
                param_min.(parameters{pstep}) = str2double(get(h_regress_para.(parameters{pstep}).minval,'String'));
                param_max.(parameters{pstep}) = str2double(get(h_regress_para.(parameters{pstep}).maxval,'String'));
            end
            param_min = [param_min.hd_a,param_min.hd_b,param_min.hd_c,param_min.volmax,param_min.abs_compress,param_min.therm_expan] ./ (10.^param_exp);
            param_max = [param_max.hd_a,param_max.hd_b,param_max.hd_c,param_max.volmax,param_max.abs_compress,param_max.therm_expan] ./ (10.^param_exp);
            param_regressed = real(fmincon(@regress_parameters,param_coeffs,[],[],[],[],param_min,param_max,[],options) .* (10.^param_exp));
        else
            param_regressed = real(fminsearch(@regress_parameters,param_coeffs,options) .* (10.^param_exp));
        end
        % 6) joyfully announce it to the world by replacing settings and
        % variables in GUI.
        
        for pstep = 1:numel(parameters)
            if get(h_regress_para.(parameters{pstep}).regress,'Value') == 1
                data.gt_sg_settings.(parameters{pstep}) = param_regressed(pstep);
            end
            set(h_regress_para.(parameters{pstep}).initial,'String',data.gt_sg_settings.(parameters{pstep}));
        end
        
        button_regress_displaydefaults
        button_regress_plot_velocities(0,0)
        
        data.gt_sg_settings.regress_history(end+1,:) = {datestr(now),data.gt_sg_settings.hd_a,data.gt_sg_settings.hd_b,data.gt_sg_settings.hd_c,data.gt_sg_settings.volmax,data.gt_sg_settings.abs_compress,data.gt_sg_settings.therm_expan,' '};
        set(h_his_tables,'Data',data.gt_sg_settings.regress_history);
        
        set(h_regress_button_process,'String','Regress parameters',...
            'CallBack',@button_regress_process); drawnow;
        
        function score = regress_parameters(param)
            % Define settings
            params = param .* (10.^param_exp);
            settings.hd_a = params(1);
            settings.hd_b = params(2);
            settings.hd_c = params(3);
            settings.volmax = params(4);
            settings.abs_compress = params(5);
            settings.therm_expan = params(6);
            settings.hull_length = data.gt_sg_settings.length;
            % Reset variables we don't want to regress
            for pstep = 1:numel(parameters)
                if get(h_regress_para.(parameters{pstep}).regress,'Value') == 0
                    settings.(parameters{pstep}) = data.gt_sg_settings.(parameters{pstep});
                end
            end
            % Run the model
            run_flightmodel(settings)
            
            % SCORE USING THE INDEX TO REMOVE BAD DATA
            score = h_regress_scorefun();
        end
        
        function stop = button_regress_plot_iteration(~,optimValues,~)
            stop = stop_regression;
            h_regress_iter_plot.YData(h_regress_iter_plot.XData == optimValues.iteration) = optimValues.fval;
            drawnow;
        end
        
        function button_regress_stop(~,~)
            stop_regression = true;
        end
        
    end

    function button_regress_auto(~,~)
        for istep = 1:str2double(get(h_regress_opt_autopass,'String'))
            set(h_regress_cost_radio(1),'Value',1);
            set(h_regress_para.hd_a.regress,'Value',1);
            set(h_regress_para.hd_b.regress,'Value',1);
            set(h_regress_para.hd_c.regress,'Value',1);
            set(h_regress_para.volmax.regress,'Value',0);
            set(h_regress_para.abs_compress.regress,'Value',0);
            set(h_regress_para.therm_expan.regress,'Value',0);
            drawnow;
            button_regress_process; drawnow;
            
            
            set(h_regress_cost_radio(2),'Value',1);
            set(h_regress_para.hd_a.regress,'Value',0);
            set(h_regress_para.hd_b.regress,'Value',0);
            set(h_regress_para.hd_c.regress,'Value',0);
            set(h_regress_para.volmax.regress,'Value',1);
            set(h_regress_para.abs_compress.regress,'Value',1);
            set(h_regress_para.therm_expan.regress,'Value',1);
            drawnow;
            button_regress_process; drawnow;
        end
        
        
        set(h_regress_cost_radio(2),'Value',1);
        set(h_regress_para.hd_a.regress,'Value',1);
        set(h_regress_para.hd_b.regress,'Value',1);
        set(h_regress_para.hd_c.regress,'Value',1);
        set(h_regress_para.volmax.regress,'Value',1);
        set(h_regress_para.abs_compress.regress,'Value',1);
        set(h_regress_para.therm_expan.regress,'Value',1);
        drawnow;
        button_regress_process; drawnow;
        
        %Change button to stop option
    end

    function button_regress_store(~,~)
        for pstep = 1:numel(parameters)
            undo.(parameters{pstep}) = data.gt_sg_settings.(parameters{pstep});
            data.gt_sg_settings.(parameters{pstep}) = str2double(get(h_regress_para.(parameters{pstep}).initial,'String'));
        end
        button_regress_displaydefaults
    end

    function button_regress_undo(~,~)
        for pstep = 1:numel(parameters)
            set(h_regress_para.(parameters{pstep}).initial,'String',undo.(parameters{pstep}));
        end
        button_regress_store(0,0)
    end

    function button_regress_reset(~,~)
        for pstep = 1:numel(parameters)
            undo.(parameters{pstep}) = str2double(get(h_regress_para.(parameters{pstep}).initial,'String'));
            set(h_regress_para.(parameters{pstep}).initial,'String',data.gt_sg_settings.(parameters{pstep}));
            set(h_regress_para.(parameters{pstep}).maxval,'String',h_regress_para.(parameters{pstep}).max);
            set(h_regress_para.(parameters{pstep}).minval,'String',h_regress_para.(parameters{pstep}).min);
        end
    end

    function button_regress_regresscheckbox(~,~)
        for fstep = 1:numel(parameters)
            if get(h_regress_para.(parameters{fstep}).regress,'Value') == 1
                set(h_regress_para.(parameters{fstep}).initial,'Enable','on');
                button_regress_constrain(0,0);
            else
                set(h_regress_para.(parameters{fstep}).initial,'Enable','off');
                set(h_regress_para.(parameters{fstep}).maxval,'Enable','off');
                set(h_regress_para.(parameters{fstep}).minval,'Enable','off');
            end
        end
        
        button_regress_displaydefaults;
    end

    function button_regress_constrain(~,~)
        for fstep = 1:numel(parameters)
            if get(h_regress_para.(parameters{fstep}).regress,'Value') == 1
                set(h_regress_para.(parameters{fstep}).maxval,'Enable',constrainRegression{get(h_regress_opt_constrain,'Value')});
                set(h_regress_para.(parameters{fstep}).minval,'Enable',constrainRegression{get(h_regress_opt_constrain,'Value')});
            end
        end
    end

    function button_regress_displaydefaults
        uicontrol(h_tab(2),'Style','text','String',...
            {['HD_A: ',num2str(data.gt_sg_settings.hd_a)],...
            ['HD_B: ',num2str(data.gt_sg_settings.hd_b)],...
            ['HD_C: ',num2str(data.gt_sg_settings.hd_c)],...
            ['V_MAX: ',num2str(data.gt_sg_settings.volmax)],...
            ['A_CMP: ',num2str(data.gt_sg_settings.abs_compress)],...
            ['T_EXP: ',num2str(data.gt_sg_settings.therm_expan)],...
            ['Hull Length: ',num2str(data.gt_sg_settings.length)]},...
            'Units','normalized','Position',[0.87 0.505 0.12 0.16]);
    end

    function button_regress_plot_velocities(~,~)
        for gstep = 1:numel(parameters)
            settings.(parameters{gstep}) = str2double(get(h_regress_para.(parameters{gstep}).initial,'String'));
        end
        settings.hull_length = data.gt_sg_settings.length;
        run_flightmodel(settings)
        
        glide_vert_spd = [data.flight(dives).glide_vert_spd];
        model_vert_spd = [data.flight(dives).model_vert_spd];
        cast_sign = sign([cast_dir{dives}]);
        depth = [data.hydrography(dives).depth];
        ind = [FRS.use_regress{dives}];
        
        mean_w_H2O = bindata2(...
            glide_vert_spd(ind) - model_vert_spd(ind),... ( = w_H20 )
            cast_sign(ind),depth(ind),... % (x and y positions)
            [-2 0 2],[0:4:1000]); % (x and y bins)
        mean_glide = bindata2(...
            glide_vert_spd(ind),...
            cast_sign(ind),depth(ind),... % (x and y positions)
            [-2 0 2],[0:4:1000]); % (x and y bins)
        mean_model = bindata2(...
            model_vert_spd(ind),...
            cast_sign(ind),depth(ind),... % (x and y positions)
            [-2 0 2],[0:4:1000]); % (x and y bins)
        
        % Plot prelim state for future comparison
        % Switch previous to gray, and delete even older one
        persistent h_regress_vis_glide_old h_regress_vis_model_old h_regress_vis_wh2o_old h_regress_vis_zero_old...
            h_regress_vis_glide h_regress_vis_model h_regress_vis_wh2o h_regress_vis_zero
        try
            delete(h_regress_vis_glide_old);
            delete(h_regress_vis_model_old);
            delete(h_regress_vis_wh2o_old);
            delete(h_regress_vis_zero_old);
            h_regress_vis_glide_old = h_regress_vis_glide;
            h_regress_vis_model_old = h_regress_vis_model;
            h_regress_vis_wh2o_old = h_regress_vis_wh2o;
            h_regress_vis_zero_old = h_regress_vis_zero;
            set(h_regress_vis_glide_old,'Color',[0.8 0.8 0.8],'LineStyle','--');
            set(h_regress_vis_model_old,'Color',[0.8 0.8 0.8],'LineStyle','--');
            set(h_regress_vis_wh2o_old,'Color',[0.8 0.8 0.8],'LineStyle','--');
        end
        h_regress_vis_wh2o = plot(h_regress_vis,[mean_w_H2O(:,1);mean_w_H2O(end:-1:1,2)],[2:4:1000,1000:-4:2],'-g');
        h_regress_vis_glide = plot(h_regress_vis,[mean_glide(:,1);mean_glide(end:-1:1,2)],[2:4:1000,1000:-4:2],'-b');
        h_regress_vis_model = plot(h_regress_vis,[mean_model(:,1);mean_model(end:-1:1,2)],[2:4:1000,1000:-4:2],'-r');
        h_regress_vis_zero = plot(h_regress_vis,[0 0],get(h_regress_vis,'YLim'),'-k');
        axis(h_regress_vis,'ij','tight');
    end

    function score = button_regress_score_minmeanwh2o
        % Score function: root mean square difference of the glide and
        % model vertical speeds. This should effectively minimise w_H20 and
        % hopefully also correct slanted w_H20 linked to thermal expansion
        % in high thermal gradient regions.
        w_H2O = abs([data.flight(dives).glide_vert_spd]-[data.flight(dives).model_vert_spd]);
        score = nanmean(w_H2O([FRS.use_regress{dives}]));
        % TODO: grid in time to avoid aliasaing of vertical sampling of
        % vertically moving features??
    end

    function score = button_regress_score_minuddiff
        % Score function: root mean square difference of the up and down
        % estimates of vertical velocity across all dives. This requires
        % gridding of all the estimates and then an RMSD.
        % 1. Grid all up/down estimates of w_H20
        glide_vert_spd = [data.flight(dives).glide_vert_spd];
        model_vert_spd = [data.flight(dives).model_vert_spd];
        cast_Dir = [cast_dir{dives}];
        depth = [data.hydrography(dives).depth];
        ind = [FRS.use_regress{dives}];
        w_H2O = bindata2(...
            glide_vert_spd(ind) - model_vert_spd(ind),... ( = w_H20 )
            cast_Dir(ind),depth(ind),... % (x and y positions)
            [-max(dives):max(dives)],[0:4:1000]); % (x and y bins)
        % w_H2O = downcasts on the left, up casts on the right, with dives
        % mirrored, so "fold" the matrix in the middle and substract the
        % two sides.
        w_H2O = w_H2O(:,max(dives):-1:1) - w_H2O(:,max(dives)+1:end);
        % TODO: Do we want to weight bins by number of samples in each bin
        % (could be useful for severely skewed profile depth distributions)
        score = nanmean(abs(w_H2O(:)));
    end

%% Tab 3 - Visualisation
h_vis_refresh = uicontrol(h_tab(3),'Style','pushbutton','String','Refresh',...
    'Units','normalized','Position',[0.86 0.94 0.13 0.05],...
    'CallBack',@button_vis_refresh);

h_vis_plot.wh2o = axes('YAxisLocation','right','Box','on',...
    'Units','normalized','Position',[0.08 0.57 0.44 .4],'Parent',h_tab(3));
h_vis_plot.wh2odiff = axes('YAxisLocation','right','Box','on',...
    'Units','normalized','Position',[0.08 0.08 0.44 .4],'Parent',h_tab(3));
h_vis_plot.N2 = axes('YAxisLocation','right','Box','on',...
    'Units','normalized','Position',[0.65 0.08 0.34 .85],'Parent',h_tab(3));

vars.abs_salinity = bindata2(...
    [data.hydrography(dives).abs_salinity],...
    [cast_dir{dives}],[data.hydrography(dives).depth],... % (x and y positions)
    [-max(dives):max(dives)],[10:4:990]); % (x and y bins)
vars.cons_temp = bindata2(...
    [data.hydrography(dives).cons_temp],...
    [cast_dir{dives}],[data.hydrography(dives).depth],... % (x and y positions)
    [-max(dives):max(dives)],[10:4:990]); % (x and y bins)
vars.pressure = bindata2(...
    [data.hydrography(dives).pressure],...
    [cast_dir{dives}],[data.hydrography(dives).depth],... % (x and y positions)
    [-max(dives):max(dives)],[10:4:990]); % (x and y bins)

X = [-max(dives):max(dives)];
for kstep=1:numel([-max(dives):max(dives)])
    try
        [vars.N2(:,kstep),~] = gsw_Nsquared(...
            vars.abs_salinity(:,kstep),...
            vars.cons_temp(:,kstep),...
            vars.pressure(:,kstep),...
            data.gps_postdive(abs(X(kstep)),1));
    end
end

    function button_vis_refresh(~,~)
        % Plot difference of up and down casts
        % 1. Grid all up/down estimates of w_H20
        w_H2O = bindata2(...
            [data.flight(dives).glide_vert_spd] - [data.flight(dives).model_vert_spd],... ( = w_H20 )
            abs([cast_dir{dives}]),[data.hydrography(dives).depth],... % (x and y positions)
            [0:max(dives)],[10:4:990]); % (x and y bins)
        w_H2O_percast = bindata2(...
            [data.flight(dives).glide_vert_spd] - [data.flight(dives).model_vert_spd],... ( = w_H20 )
            [cast_dir{dives}],[data.hydrography(dives).depth],... % (x and y positions)
            [-max(dives):max(dives)],[10:4:990]); % (x and y bins)
        % w_H2O = downcasts on the left, up casts on the right, with dives
        % mirrored, so "fold" the matrix in the middle and substract the
        % two sides.
        w_H2O_diff = w_H2O_percast(:,max(dives):-1:1) - w_H2O_percast(:,max(dives)+1:end);
        
        h_vis_plot.wh2o_plot = pcolor(h_vis_plot.wh2o,[1:max(dives)],[10:4:989],real(w_H2O));
        set(h_vis_plot.wh2o,'YDir','reverse')
        xlabel(h_vis_plot.wh2o,'Dive Number'); ylabel(h_vis_plot.wh2o,'Depth (m)');
        set(h_vis_plot.wh2o_plot,'EdgeColor','none')
        h_vis_plot.wh2o_cbar = colorbar(h_vis_plot.wh2o); ylabel(h_vis_plot.wh2o_cbar,{'Estimated W_H_2_O'});
        caxis(h_vis_plot.wh2o,[-1 1]*0.05);
        
        h_vis_plot.wh2odiff_plot = pcolor(h_vis_plot.wh2odiff,[1:max(dives)],[10:4:989],real(w_H2O_diff));
        set(h_vis_plot.wh2odiff,'YDir','reverse')
        xlabel(h_vis_plot.wh2odiff,'Dive Number'); ylabel(h_vis_plot.wh2odiff,'Depth (m)');
        set(h_vis_plot.wh2odiff_plot,'EdgeColor','none')
        h_vis_plot.wh2odiff_cbar = colorbar(h_vis_plot.wh2odiff); ylabel(h_vis_plot.wh2odiff_cbar,{'Absolute difference between up- and downcast estimates of W_H_2_O','Ideally, it should appear as random noise dissimilar T&S structures or dive profile shape.'});
        caxis(h_vis_plot.wh2odiff,[-0.05 0.05]);
        
        w_H2O_N2 = w_H2O_percast(1:end-1,:) + diff(w_H2O_percast);
        h_vis_plot.N2_plot = plot(h_vis_plot.N2,...
            vars.N2(:),...
            w_H2O_N2(:).^2,...
            '-k','Marker','.','LineStyle','none');
        xlabel(h_vis_plot.N2,'N^2 (s^-^2)'); ylabel(h_vis_plot.N2,'W_H_2_O^2 (cm^2 s^-^2)');
        set(h_vis_plot.N2,'XLim',[0 5e-4],'YLim',[0 0.003]);
    end

% W vs Temp, Depth, Time
% Up vs Down
% Histograms like in Frajka Williams

%% Tab 4 - History
h_param_save = uicontrol(h_tab(4),'Style','pushbutton','String','Save notes',...
    'Units','normalized','Position',[0.86 0.94 0.13 0.05],...
    'CallBack',@button_his_save);

columnname = {'Date','hd_a','hd_b','hd_c','volmax','abs_compress','therm_expan','Notes'};
columnformat = {'char','char','char','char','char','char','char','char'};

if ~isfield(data.gt_sg_settings,'regress_history')
    data.gt_sg_settings.regress_history = {datestr(now),data.gt_sg_settings.hd_a,data.gt_sg_settings.hd_b,data.gt_sg_settings.hd_c,data.gt_sg_settings.volmax,data.gt_sg_settings.abs_compress,data.gt_sg_settings.therm_expan,' '};
end

h_his_tables = uitable(h_tab(4),'Data', data.gt_sg_settings.regress_history,...
    'ColumnName', columnname,...
    'ColumnFormat', columnformat,...
    'ColumnEditable', [false false false false false false false true],...
    'ColumnWidth',{200 150 150 150 150 150 150 200},...
    'RowName',[],...
    'Units','normalized','Position',[0.05 0.05 0.9 0.85],...
    'CellEditCallback',@button_his_save);

    function button_his_save(~,~)
        data.gt_sg_settings.regress_history = get(h_his_tables,'Data');
    end

%% Tab 5 - HEEEEEEEEELP ME !
h_auto_process = uicontrol(h_tab(5),'Style','pushbutton','String','Automatic regression',...
    'Units','normalized','Position',[0.86 0.94 0.13 0.05],...
    'CallBack',@button_regress_auto);

uicontrol('Parent',h_tab(5),'Style','text','String',HelpInfo,...
    'Units','normalized','Position',[0.05 0.05 0.5 0.9])

%% Subfunctions
    function out = indicesIn(list,sublist)
        [~, out, ~] = intersect(list,sublist);
    end

% Create get_exp handle to make all regressed coefficients of similar
% units (works better).
get_exp = @(x) floor(log10(x));
get_coeff = @(x) x.*10.^-get_exp(x);

%% Flight model scripts
%%%%%%%%%%%%%% START OF FLIGHT MODEL ROUTINES %%%%%%%%%%%%%%
    function run_flightmodel(settings)
        for istep = intersect([data.eng.dive],[data.log.dive])
            % settings. should have volmax,hd_a,hd_b,hd_c,abs_compress,therm_expan,hull_length
            % Determine glider buoyancy based on volmax and in situ density
            % Calculate volume
            vbd = data.eng(istep).vbdCC; % cm^3
            vol0 = settings.volmax + (data.log(istep).C_VBD - data.gt_sg_settings.vbd_min_cnts)/data.gt_sg_settings.vbd_cnts_per_cc; % cm^3
            vol = (vol0 + vbd) .* ...
                exp(-settings.abs_compress * data.hydrography(istep).pressure + settings.therm_expan * (data.hydrography(istep).temp - data.gt_sg_settings.temp_ref)); % cm^3
            
            % Calculate buoyancy of the glider
            % Buoyancy units
            %   kg + (kg.m-3 * cm^3 * 10^-6)
            % = kg + (kg.m-3 * m^3)
            % = kg
            data.flight(istep).buoyancy =  -data.gt_sg_settings.mass + ( data.hydrography(istep).rho .* vol * 1.e-6 ); % kg
            
            % Calculate glide slope and speed
            [data.flight(istep).model_spd, data.flight(istep).model_slope]... % m.s-1, degrees
                = flightvec(...
                data.flight(istep).buoyancy,... % kg
                data.eng(istep).pitchAng,... % degrees
                settings.hull_length,... % Hull length in meters, excluding antenna
                settings.hd_a,... % rad^-1
                settings.hd_b,... % m^1/4 . kg^1/4 . s^-1/2
                settings.hd_c,... % rad^-2
                data.gt_sg_settings.rho0,...
                istep); % kg.m-3
            
            % Determine vertical and horizontal speed
            data.flight(istep).model_horz_spd = data.flight(istep).model_spd.*cos(data.flight(istep).model_slope*pi/180);
            data.flight(istep).model_vert_spd = data.flight(istep).model_spd.*sin(data.flight(istep).model_slope*pi/180);
        end
    end
    function [ umag, thdeg ] = flightvec( bu, ph, xl, a, b, c, rho0, istep)
        % ********************************************************************************
        %	function flightvec0(bu, ph, xl, a, b, c, rho0 )
        %		Solves unaccelerated flight equations iteratively
        %		for speed magnitude umag and glideangle thdeg
        % ********************************************************************************
        
        gravity = gsw_grav(data.gps_postdive(istep,1),...
            data.hydrography(istep).pressure); % m.s-2
        
        th = (pi/4)*sign(bu);	% initial slope set to 1 or -1
        
        buoyforce = gravity .* bu; % kg.m.s-2
        
        q = ( sign(bu).*buoyforce./(xl*xl*b) ).^(4/3); 	% dynamic pressure for vertical flight
        alpha = 0.;
        
        tol = 0.001;
        j = 0;
        q_old = zeros(size(bu));
        param = ones(size(bu));
        valid = find( bu ~= 0 & sign(bu).*sign(ph) > 0 );
        umag = zeros(size(bu));
        thdeg = zeros(size(bu));
        
        while(~isempty( find( abs( (q(valid)-q_old(valid))./q(valid) ) > tol )) & j <= 15);
            
            q_old = q;
            param_inv = a*a*tan(th).*tan(th).*q.^(0.25)/(4*b*c);
            
            valid = find( param_inv > 1 & sign(bu).*sign(ph) > 0);	% valid solutions for param < 1
            
            param(valid) = 4*b*c./(a*a*tan(th(valid)).*tan(th(valid)).*q(valid).^(0.25));
            q(valid) = ( buoyforce(valid).*sin(th(valid))./(2*xl*xl*b*q(valid).^(-0.25)) ).*  ...
                (1 + sqrt(1-param(valid)));
            alpha(valid) = ( -a*tan(th(valid))./(2*c) ).*(1 - sqrt(1-param(valid)));
            thdeg(valid) = ph(valid) - alpha(valid);
            
            stall = find( param_inv <= 1 | sign(bu).*sign(ph) < 0 );	% stall if param >= 1
            
            % TODO: REMOVE ME
            stall = [];
            %
            
            q(stall) = 0.;
            thdeg(stall) = 0.;
            
            th = thdeg*pi/180;
            
            j = j+1;
            % TODO: Should this be replaced with actual density.....??
            umag = sqrt( 2*q/rho0 );
        end
    end
%%%%%%%%%%%%%% END OF FLIGHT MODEL ROUTINES %%%%%%%%%%%%%%

%% WAIT FOR GUI CLOSE
uiwait(h_gui)
end

function [ym,yb] = bindata2(y,x1,x2,x1rg,x2rg)
%function [ym,yb] = bindata2(y,x1,x2,x1rg,x2rg)
%Computes:
%ym(ii,jj) = mean(y(x1>=x1rg(ii) & x1 < x1rg(ii+1) & x2>=x2rg(jj) & x2 < x2rg(jj+1))
%for every ii, jj
%If a bin is empty it returns nan for that bin
%using a fast algorithm which uses no looping
%Also returns yb, the approximation of y using binning (useful for r^2
%calculations). Example:
%
%x = randn(500,2);
%y = sum(x.^2,2) + randn(500,1);
%xrg = linspace(-3,3,10)';
%[ym,yb] = bindata2(y,x(:,1),x(:,2),xrg,xrg);
%subplot(1,2,1);plot3(x(:,1),x(:,2),y,'.');
%subplot(1,2,2);h = imagesc(xrg,xrg,ym);
%set(h,'AlphaData',~isnan(ym)); box off;
%
%By Patrick Mineault
%Refs: http://xcorr.net/?p=3326
%      http://www-pord.ucsd.edu/~matlab/bin.htm
% Modified by Bastien Y. Queste to account for NaNs and to permute the matrix at the end.

good = ~isnan(x1+x2+y);

[~,whichedge1] = histc(x1(good),x1rg(:)');
[~,whichedge2] = histc(x2(good),x2rg(:)');

bins1 = min(max(whichedge1,1),length(x1rg)-1);
bins2 = min(max(whichedge2,1),length(x2rg)-1);

bins = (bins2-1)*(length(x1rg)-1)+bins1;

xpos = ones(size(bins,1),1);
ns = sparse(bins,xpos,1,(length(x1rg)-1)*(length(x2rg)-1),1);
ysum = sparse(bins,xpos,y(good),(length(x1rg)-1)*(length(x2rg)-1),1);
ym = full(ysum)./(full(ns));
yb = ym(bins);
ym = reshape(ym,length(x1rg)-1,length(x2rg)-1)';
end