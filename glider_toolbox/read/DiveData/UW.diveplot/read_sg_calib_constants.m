function [calib_consts] = read_sg_calib_constants()
  % read calibration constants and return as a structure ala python
    sg_calib_constants
	% Assign no local variables before this point!!
	% Thus the only variables who() sees are from sg_calib_constants
    x__vars = who;
    for x__i = 1:size(x__vars,1)
      x__vname = x__vars{x__i};
      eval(sprintf('calib_consts.%s = %s;',x__vname,x__vname));
    end
    calib_consts.fileinfo = dir('sg_calib_constants.m');
    