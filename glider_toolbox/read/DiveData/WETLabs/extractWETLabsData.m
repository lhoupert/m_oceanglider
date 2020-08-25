% *************************************************************************
% The purpose of this function is to extract WETLabs Sensor Data
% from the structure containing the data from the engineering file.  This
% function currently only processes the backscatter, chlorophyll and CDOM
% values.  Temperature is currently not utilized for plotting.
%
% The Scattering, CDOM and Chlorophyll fields are used from the calData
% parameter to perform scaling on the CDOM, Scattering and Chlorophyll
% fields in the engineering data structure.
%
% The resultant WETLabsData structure is a nested data structure referenced
% by sensor Id in the first dimension and signal id in the nested
% dimension.  The child entries of the structure are:
%  refData - Signal reference data
%  sigData - Measured signal value
%  depth - Depth measurement for each array entry(ie. refData, sigData)
%  iwn - Indexes of points for downcast
%  iwp - Indexes of points for upcast
%  scaled - Scaled data values produced from calData and sigData
%  units - String description of the output units
%
% A warning will be displayed in the console if calibration constants are
% not available for a specific sensor.
%
% @param eng Engineering data file
% @param mp Number of rows in the dataset
% @param descriptors List of string headers for relevant columns
% @param calData WETLabs calibration data set.
% *************************************************************************
function [WETLabsData] = extractWETLabsData(eng, mp, descriptors, calData)
    WETLabsData = [];

    for i=1:length(descriptors)
        % PARSE OUT THE SENSOR ID, SIGNAL ID AND SIGNAL TYPE FROM THE DESCRIPTOR STRING IN THE ENGINEERING DATA COLUMN HEADER
        [sensorId signalId signalType] = parseSensorInformation(descriptors(i));

        if (isempty(sensorId) == 0) && (isempty(signalType) == 0)
            % IF THE SG_CALIB_CONSTANTS FILE IS MISSING WET Labs CALIBRATION DATA, WE NEED TO PREVENT THE SOFTWARE FROM
            % CRASHING DUE TO ANY ATTEMPTS TO INDEX INTO AN EMPTY ARRAY. THEREFORE, VERIFY WHETHER OR NOT DATA EXISTS. IF THE ARRAY
            % HAS DATA, THEN INDEX INTO THE REQUIRED DATA. IF NOT, THEN SET THE "POINTER" TO NULL. THIS WILL PREVENT THE ACT OF
            % INDEXING INTO AN EMPTY ARRAY FROM CAUSING A CRASH.
            if isempty(calData) == 1
                pCalData = [];
            else
                %pCalData = calData(sensorId);
                if strcmpi(sensorId, 'wlbbfl2') && isfield(calData, 'wlbbfl2')
                    pCalData = calData.wlbbfl2;
                elseif strcmpi(sensorId, 'wlbb2fl') && isfield(calData, 'wlbb2fl')
                    pCalData = calData.wlbb2fl;
                elseif strcmpi(sensorId, 'wlbb3') && isfield(calData, 'wlbb3')
                    pCalData = calData.wlbb3;
                elseif strcmpi(sensorId, 'wlfl3') && isfield(calData, 'wlfl3')
                    pCalData = calData.wlfl3;
                elseif isfield(calData, 'wl')
                    % Old approach - puck type not included in calibration
                    % fields.
                    pCalData = calData.wl;
                else    
                    pCalData = [];
                end
            end
            
            if strcmpi(signalType, 'sig')
                if (isempty(signalId) == 1)
                    [signalId] = generateSignalID(WETLabsData, sensorId, signalId, 'Scattering');
                end
                %[WETLabsData{sensorId}.Scattering{signalId}] = handleWLData(eng, mp, pCalData, 'Scatter', '', sensorId, signalId,...
                [WETLabsData.(sensorId).Scattering{signalId}] = handleWLData(eng, mp, pCalData, 'Scatter', '', sensorId, signalId,...
                                                                            '\beta(\theta_c) (m^-^1 sr^-^1)', '\beta(\theta_c) (A/D Counts)', descriptors(i)); %#ok<AGROW>

            elseif (strcmpi(signalType, 'BB1sig') || strcmpi(signalType, 'BB2sig') || strcmpi(signalType, 'BB3sig')  || ...
                    strcmpi(signalType, 'BBsig'))
                if (isempty(signalId) == 1)
                    [signalId] = generateSignalID(WETLabsData, sensorId, signalId, 'Scattering');
                end
                %[WETLabsData{sensorId}.Scattering{signalId}] = handleWLData(eng, mp, pCalData, 'Scatter', '', sensorId, signalId,...
                [WETLabsData.(sensorId).Scattering{signalId}] = handleWLData(eng, mp, pCalData, 'Scatter', '', sensorId, signalId,...
                                                                            '\beta(\theta_c) (m^-^1 sr^-^1)', '\beta(\theta_c) (A/D Counts)', descriptors(i)); %#ok<AGROW>

            elseif (strcmpi(signalType, 'FL1sig') || strcmpi(signalType, 'FL2sig') || strcmpi(signalType, 'FL3sig') || ...
                    strcmpi(signalType, 'Chlsig') || strcmpi(signalType, 'Cdomsig'))
                % Use corresponding 'ref' field (FL1ref, FL2ref, FL3ref, Chlref, Cdomref) to
                % determine the type of fluorescence data.
                [flType] = getFluorType(eng, mp, descriptors(i), pCalData); %calData);
                if (strcmpi(flType, 'Chl'))
                    if (isempty(signalId) == 1)
                        [signalId] = generateSignalID(WETLabsData, sensorId, signalId, 'Chlorophyll');
                    end

                    %[WETLabsData{sensorId}.Chlorophyll{signalId}] = handleWLData(eng, mp, pCalData, 'Chl', 'Chl', sensorId,...
                    [WETLabsData.(sensorId).Chlorophyll{signalId}] = handleWLData(eng, mp, pCalData, 'Chl', 'Chl', sensorId,...
                                                                    signalId, '{\mu}g/l', '(A/D Counts)', descriptors(i));
                elseif (strcmpi(flType, 'Cdom'))
                    if (isempty(signalId) == 1)
                        [signalId] = generateSignalID(WETLabsData, sensorId, signalId, 'CDOM');
                    end

                    %[WETLabsData{sensorId}.CDOM{signalId}] = handleWLData(eng, mp, pCalData, 'CDOM', 'Cdom', sensorId, ...
                    [WETLabsData.(sensorId).CDOM{signalId}] = handleWLData(eng, mp, pCalData, 'CDOM', 'Cdom', sensorId, ...
                                                            signalId, 'ppb', '(A/D Counts)', descriptors(i));
                elseif (strcmpi(flType, 'Phycoerythrin'))
                    if (isempty(signalId) == 1)
                        [signalId] = generateSignalID(WETLabsData, sensorId, signalId, 'Phycoerythrin');
                    end
                
                    %[WETLabsData{sensorId}.Phycoerythrin{signalId}] = handleWLData(eng, mp, pCalData, 'Phycoerythrin', 'Phycoerythrin',...
                    [WETLabsData.(sensorId).Phycoerythrin{signalId}] = handleWLData(eng, mp, pCalData, 'Phycoerythrin', 'Phycoerythrin',...
                                                                sensorId, signalId, 'ppb', '(A/D Counts)', descriptors(i));
                else
                    display(sprintf('Unknown fluorescence type <%s>.', flType));
                end
                
            end
        end
    end

%     WETLabsData = WETLabsData(~cellfun('isempty',WETLabsData));
end

function [flType] = getFluorType(eng, mp, descriptor, calData)

    strSig = [];
    strRef = [];
    name = 'none';
    calField = 'UNKNOWN';
    flType = 'Unknown';
    
    % Find specified fluorometer signal name (wlXXXX_FL1sig, wlXXXX_Chlsig, etc.) in eng structure
    list = strcmpi(fieldnames(eng), descriptor);
    idx = find(list > 0);
    if isempty(idx) == 0
        strSig = descriptor;
        % Use signal name to generate reference name (FL1ref, FL2ref, FL3ref, Chlref, Cdomref)
        strRef = regexprep(strSig, 'sig', 'ref', 'preservecase');
        % Extract the reference data from eng.
        refData = fillchk(eng, strRef, mp);
        if ~isempty(refData)
            % First value tells us the wavelength of fluorescence data
            fl = refData(1);
            % Map wavelength to appropriate string.
            % See http://www.wetlabs.com/eco-triplet, Specifications tab.
            if (isfield(calData, 'Chlorophyll') && isfield(calData.Chlorophyll, 'wavelength') && fl == calData.Chlorophyll.wavelength)
                flType = 'Chl';             % 695
            elseif (isfield(calData, 'CDOM') && isfield(calData.CDOM, 'wavelength') && fl == calData.CDOM.wavelength)
                flType = 'Cdom';            % 460
            elseif (isfield(calData, 'Phycoerythrin') && isfield(calData.Phycoerythrin, 'wavelength') && fl == calData.Phycoerythrin.wavelength)
                flType = 'Phycoerythrin';   % 570 - warning, same as Rhodamine
            elseif (isfield(calData, 'Uranine') && isfield(calData.Uranine, 'wavelength') && fl == calData.Uranine.wavelength)
                flType = 'Uranine';         % 530
            elseif (isfield(calData, 'Rhodamine') && isfield(calData.Rhodamine, 'wavelength') && fl == calData.Rhodamine.wavelength)
                flType = 'Rhodamine';       % 570 - warning, same as Phycoerythrin
            elseif (isfield(calData, 'Phycocyanin') && isfield(calData.Phycocyanin, 'wavelength') && fl == calData.Phycocyanin.wavelength)
                flType = 'Phycocyanin';     % 680
            else
                % the wavelength field was probably not included in
                % sg_calib_constants.m when the .nc file was created.
                
                % output a message
                display(sprintf('Error - unknown WETLabs fluorescence wavelength %d', fl));
                          
                % try to recover using some well-known values
                if (fl == 695)
                    flType = 'Chl';
                elseif (fl == 460)
                    flType = 'Cdom';
                elseif (fl == 570)
                    flType = 'Phycoerythrin';
                elseif (fl == 530)
                    flType = 'Uranine';
                %elseif (fl == 570)    % 570 - same as Phycoerythrin
                %    flType = 'Rhodamine';
                elseif (fl == 680)
                    flType = 'Phycocyanin';
                else
                    flType = 'unknown';
                end
            end

        end
    end
end

% *************************************************************************
% This is an internal function to this source file.
%
% This function performs the scaling operation on WETLabs triplet puck
% signals using the provided counts, darkCount and scaleFactor.
% *************************************************************************
function [scaled] = scaleWETlabsData(counts, darkCounts, scaleFactor)
    scaled = scaleFactor * (counts - darkCounts);
end

% *************************************************************************
% This is an internal function to this source file.
%
% This function is responsible for extracting the signal and reference data
% from the specified eng file using sigType, sensorId and signalId.  A
% default description, specified in description, will also be stored in the
% units field of the structure.
%
% Data is pulled from the engineering data structure by accessing the
% fields 'wl<sensor Id>_<signal type><ref AND sig><signalId>.  For example
% if the data for signal two from sensor 1 which is a backscatter is
% accessed the following elements would be retrieved and used for
% processing: wl1_bb2ref AND wl1bb2sig
%
% The calibration structure is built using signal specific structure names.
%  In order to generalize the caller of the function specifies the top
%  level cal structure for the sensor and also the member name to retrieve
%  the signal specific data from.  To do this the calField name from the
%  calData element is retrieved and subsequently the scaleFactor and
%  darkCounts for the specific signal are retrieved.
%
% @param eng Engineering data structure containing data
% *************************************************************************
function [entry] = handleWLData(eng, mp, calData, calField, sigType, sensorId, signalId, scaledDescription, rawDescription, descriptor)

    strSig = [];
    strRef = [];
    signalCalData = [];
    entry = [];

    % ------------------------------------------------------------------------------------------------------------------------------
    % Retrieve the reference column of data
    % ------------------------------------------------------------------------------------------------------------------------------

    list = strcmpi(fieldnames(eng), descriptor);
    idx = find(list > 0);
    if isempty(idx) == 0
        strSig = descriptor;
    end

    % THIS LINE CREATES THE NAME OF THE REFERENCE FIELD BY REPLACING THE SIGNAL FIELD NAME WITH THE REFERENCE FIELD NAME
    if isempty(strSig) == 0
        strRef = regexprep(strSig, 'sig', 'ref', 'preservecase');
    end

    % EXTRACT THE DATA FROM THE ENGINEERING FILE STRUCTURE
    refData = fillchk(eng, strRef, mp);
    sigData = fillchk(eng, strSig, mp);

    if ~isempty(sigData) && ~isempty(refData) 
        % FIND ANY EMPTY VALUES REPRESENTED BY NAN'S
        idx = find(~isnan(sigData));

        % PLACE THE WET LABS DATA INTO THE WET LABS STRUCTURE
        entry = createWLEntry(idx, eng, refData, sigData, rawDescription);

        % LOCATE THE CALIBRATION DATA FOR THE SENSOR IN THE CALIBRATION DATA STRUCTURE
        signalCalData = [];
        if ~isempty(calData)
            calFieldNames = fieldnames(calData);
            idx = regexpi(calFieldNames, calField);
            idx = find(not(cellfun('isempty', idx)));

            % IF WE HAVE MORE THAN ONE SENSOR OF A PARTICULAR TYPE (USUALLY BACKSCATTER), EXTRACT OUT THE CALIBRATION DATA FOR THE
            % PROPER WAVELENGTH
            if (strcmpi(calField, 'scatter') == 1)
                % Use the first value in the refData array to generate the
                % calibration field prefix, e.g. Scatter700.
                idx = regexpi(calFieldNames, sprintf('%s%3.3d', calField, refData(1)));
            else
                idx = regexpi(calFieldNames, calField);
            end
            idx = find(not(cellfun('isempty', idx)));

            if ~isempty(idx)
                if isfield(calData, (char(calFieldNames(idx))))
                    signalCalData = calData.(char(calFieldNames(idx)));
                end
            else
                %display(sprintf('Error - WET Labs calibration constants missing for %s %s%3.3d ?', sensorId, calField, refData(1)));
                display(sprintf('Error - WET Labs calibration constants missing for %s, %s%03d ?', sensorId, calField, refData(1)));
                %display(sprintf('(update sg_calib_constants.m on base station and reprocess the dive)'));
            end
        end
    else
        display(sprintf('Unable to extract WET Labs data from engineering file.'));
    end
        
    % ATTEMPT CONVERSION OF THE RAW SENSOR VALUES INTO ENGINEERING UNITS. IF THE ATTEMPT FAILS, REPORT THE ERROR AND CONTINUE.
    if ~isempty(signalCalData)
        [entry.scaled] = scaleWETlabsData(entry.counts,signalCalData.scaleFactor,signalCalData.darkCounts);
        entry.units = scaledDescription;
    else
        % jfaust - we already output a good error msg above, so don't
        % confuse the user
        %display(sprintf('WETLabs Calibration Data Error For %u.%u.  darkCounts or scaleFactor calibration constants are missing', sensorId, signalId));
        %display(sprintf('WETLabs Calibration Data Error For %s.%s.  darkCounts or scaleFactor calibration constants are missing', sensorId, char(strSig)));
    end
end

% *************************************************************************
% This is an internal function to this source file.
%
% This function parses the sensor id, signal id and signal type from the
% specified descriptor string.
% *************************************************************************
function [sensorId, signalId, signalType] = parseSensorInformation(descriptor)
    signalType = [];
    signalId   = [];

    % Wetlabs field names are constructed using _'s between the relevant
    % pieces of data so we tokenize them using the _ to get to the actual
    % field values
    [prefix remain] = strtok(descriptor, '_');

    % The WETLabs fields start with the prefix value from the .cnf file.
    % This gives us the puck type (wlbbfl2, wlbb2fl, wlbb3, wlfl3).
    sensorId = char(prefix);
    
    % Look for the wet labs field names.  These can be bbsig, bbref,
    % Chlreg, etc. as indicated by the below regular expression
    %
    % jf - this is old - probably not needed
    [tokens, ~] = regexpi(char(remain), '(Chlref|Chlsig|Cdomsig|Cdomref|BBsig|BBref|temp)(\d+)', 'tokens', 'match');

    if isempty(tokens) == 0
        signalType = tokens{:}{1};
        signalId = str2double(tokens{:}{2});
    else
        % current expected values are:
        %  BB1ref, BB1sig
        %  BB2ref, BB2sig
        %  BB3ref, BB3sig
        %  FL1ref, FL1sig
        %  FL2ref, FL2sig
        %  FL3ref, FL3sig
        % Anything else is for backward compatibility 
        [tokens, ~] = regexpi(char(remain), '(Chlref|Chlsig|Cdomsig|Cdomref|BBref|BBsig|BB1ref|BB1sig|BB2ref|BB2sig|FL1ref|FL1sig|FL2ref|FL2sig|FL3ref|FL3sig)', 'tokens', 'match');
        if isempty(tokens) == 0
            signalType = char(tokens{:});
        else
            [tokens, ~] = regexpi(char(remain), '(\d+)(sig|ref|temp)', 'tokens', 'match');
            if isempty(tokens) == 0
                signalType = tokens{:}{2};
            end
        end
                
    end
end

% *************************************************************************
% This is an internal function to this source file.
%
% This function packages the specified engineering data into a valid
% WETLabs data structure.
% *************************************************************************
function entry = createWLEntry(idx, eng, refData, sigData, description)
    entry.refData = refData(idx);
    entry.counts  = sigData(idx);
    entry.depth   = eng.depth(idx);
    entry.scaled  = [];
    entry.units   = description;

    entry.iwn     = find(eng.pitchCtl(idx) < 0);
    entry.iwp     = find(eng.pitchCtl(idx) >= 0);
%     %------------------------------------------------------------------------
%     % DIFFERENTIATE BETWEEN DOWNCAST DATA AND UPCAST DATA
%     %------------------------------------------------------------------------
%     % FIND INDEXES FOR DOWNCAST/DESECENT DATA POINTS
%     entry.iwn = 1:find(eng.depth==max(eng.depth));
%     
%     % FIND INDEXES FOR UPCAST/ASCENT DATA POINTS
%     entry.iwp = find(eng.depth==max(eng.depth)):length(eng.depth);
% 
%     % 1. SEARCH FOR NEGATIVE PITCH CONTROL VALUES AFTER MAX DEPTH IS IDENTIFIED
%     % 2. REMOVE ANY DATA ON THE UPCAST WHERE THE PITCHCTL IS NEGATIVE. THIS IS A CARTE BLANCHE APPROACH FOR NOW. IT WILL REMOVE ANY
%     %    DATA ON THE UPCAST WITH A NEGATIVE PITCH CONTROL VALUE ASSOCIATED WITH IT.
%     %
%     % NOTE: THIS IS ONLY AN ISSUE WITH ALI SENSORS SINCE THEY LOG DATA INDEPENDENT OF THE DATA LOGGED IN THE ENGINEERING FILE.
%     if ~isempty(entry.iwp)
%         pitchPos = find(eng.pitchCtl(entry.iwp)<0);
%         if ~isempty(pitchPos)
%             entry.iwp(pitchPos) = [];
%         end
%     end
end

% **********************************************************************************************************************************
% Function: validateSignalID
%
% Purpose: the purpose of this function is to generate a validate signal ID when the signal ID is missing due to omission from the
%          column header; in particular when dealing with older data (i.e., pre-66.07.14).
% **********************************************************************************************************************************
function [signalId] = generateSignalID(WETLabsData, sensorId, signalId, strField)
    % jfaust - i believe the signalId is simply an array index for the
    % different sets of backscatter data for a given puck (e.g., 470 nm vs 700 nm).
    if isempty(WETLabsData) || isfield(WETLabsData, sensorId) == 0 || isfield(WETLabsData.(sensorId), strField) == 0
        signalId = 1;
    else
        signalId = length(WETLabsData.(sensorId).(strField)) + 1;
    end
end
