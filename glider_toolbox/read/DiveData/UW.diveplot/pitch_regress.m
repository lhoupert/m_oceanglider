% 
% Copyright (c) 2006-2013 by University of Washington.  All rights reserved. Confidential
%
%
function [AA] = pitch_regress(ip,pitch,pitch_control,buoy,vbd,roll_control,pitch_intended,vol_comp)
	% reduce the arrays to valid points
	pitch = pitch(ip); % observed pitch response to the following factors:
	pitch_control = pitch_control(ip); % pitch mass position wrt center
	buoy = buoy(ip); % impact of seawater density, vbd system and compressee, if any
	vbd = vbd(ip); % oil shift during pumps
	roll_control = roll_control(ip); % pitch/roll coupling post-TTI mass-shifter
	pitch_intended = pitch_intended(ip); % general attack angle
	vol_comp = vol_comp(ip); % volume change of compressee based on T and P

	% avoid rank deficiences
	% check for single pitch_intended
	n_pitch_desired = length(find( abs(diff(unique(abs(pitch_intended)))) > 1)); % use 1 degree
	n_pitch_desired = n_pitch_desired + 1; % how many actual different pitch_desired there are, unless pitch_intended is empty
	if (n_pitch_desired == 1)
		if (isempty(find(vol_comp,1)))
			% the last terms can be dropped; dodge rank deficiency warning, which can scare the natives
			if (isempty(find(roll_control,1)))
				XX = [ones(size(pitch_control)) pitch_control buoy vbd];
				AA = XX\pitch;
				AA = [AA; 0; 0; 0];
			else
				XX = [ones(size(pitch_control)) pitch_control buoy vbd roll_control];
				AA = XX\pitch;
				AA = [AA; 0; 0];
			end
			return;
		end
		pitch_intended = zeros(size(pitch_control)); % effectively drop this term
	end
	XX = [ones(size(pitch_control)) pitch_control buoy vbd roll_control pitch_intended vol_comp];
	AA = XX\pitch;
    %DEBUG  AA = single(XX)\single(pitch); % experiment to test onboard possibility
