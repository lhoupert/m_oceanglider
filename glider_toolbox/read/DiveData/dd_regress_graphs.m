%==========================================================================
% COPYRIGHT
% © 2014 Kongsberg Underwater Technology Inc. All rights reserved. No part 
% of this work covered by the copyright hereon may be reproduced or 
% otherwise copied without prior permission from Kongsberg Underwater 
% Technology, Inc.
%==========================================================================
%
% Description: Launch regression plots
%

% regress_dives = { selected_dives };

existingPath = pwd;
cd(settings.pathname);

% Saves as it gets cleared during execution of regress_nnn scripts
noDivesToRegress = regress_dives; 

try

    count = sum(settings.plotSelectionState(15:17));
    progressPercent = .70;
    progressIncrement = .20/count;

    if count > 0
        hWaitb = waitbar(progressPercent, 'Completed computing regression plots.');       
    end
    
    % --------------------------------------------------------------------------------
    % Plot 15 - Regression - Pitch (regression plot)
    % --------------------------------------------------------------------------------
    if settings.plotSelectionState(15) == 1
        waitbar(progressPercent, hWaitb, 'Computing regression for pitch data...');
        regress_dives = noDivesToRegress;
        regress_pitch;
        set(gcf, 'Name','Pitch Regression Plot (Attitude)');
        set(gcf, 'UserData', 'Plot_RegressPitch');       
        progressPercent = progressPercent+progressIncrement;
    end


    % --------------------------------------------------------------------------------
    % Plot 16 - Regression - Roll (regression plot)
    % --------------------------------------------------------------------------------
    if settings.plotSelectionState(16) == 1
        waitbar(progressPercent, hWaitb, 'Computing regression for roll data...');
        regress_dives = noDivesToRegress;
        regress_roll;

        try
            % HACK: Might be possible gcf-1 does not exist in some circumstances
            hFig1 = gcf-1;
            sFig1 = get(hFig1);
            set(hFig1,'Name', 'Roll Regression Plot (Attitude)'); 
            set(hFig1, 'UserData', 'Plot_RegressRoll_Attitude');
        catch
        end
        
        set(gcf, 'Name', 'Roll Regression plot (Turn Rate)');
        set(gcf, 'UserData', 'Plot_RegressRoll_Turn_Rate');
        progressPercent = progressPercent+progressIncrement;
    end


    % --------------------------------------------------------------------------------
    % Plot 17 - Regression - VBD (regressionplot)
    % --------------------------------------------------------------------------------

    % Only accessible via the Pregress VBD button for now.
    
    % if settings.plotSelectionState(17) == 1
    %     % Supress prompt in the following script (set to some invalid value, this 
    %     % configuration is not used by KUTI)
    %     % mass_comp = ?; 
    % 
    %     waitbar(progressPercent, hWaitb, 'Computing regression for VBD data...');
    %     regress_dives = noDivesToRegress;
    % 
    %     % KUTI modified version of regress_vbd.m :       
    %     [vbdBias, w_rms_final] = dd_regress_vbd(handles, reference_c_vbd, [diveNumberList], settings);
    % 
    %     set(gcf, 'Name','VBD Regression Plot (turn rate)');
    %     set(gcf, 'UserData', 'Plot_RegressVBD');
    %     progressPercent = progressPercent+progressIncrement;
    % end

    if count > 0
        waitbar(progressPercent, hWaitb, 'Completed computing regression plots.');      
    end
catch
end

% Restore path
cd(existingPath);


