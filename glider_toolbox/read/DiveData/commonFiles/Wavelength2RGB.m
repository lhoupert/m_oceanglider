% ****************************************************************************************************
% Converts the given wavelength in the visible spectrum into red, green and blue components
% ****************************************************************************************************
function [R G B] = Wavelength2RGB(wavelength)
    if((wavelength >= 380.0) && (wavelength <= 440.0))
        R = -1.0*(wavelength-440.0)/(440.0-380.0);
        G = 0.0;
        B = 1.0;
    elseif ((wavelength >= 440.0) && (wavelength <= 490.0))
        R = 0.0;
        G = (wavelength-440.0)/(490.0-440.0);
        B = 1.0;
    elseif ((wavelength >= 490.0) && (wavelength <= 510.0))
        R = 0.0;
        G = 1.0;
        B = -1.0*(wavelength-510.0)/(510.0-490.0);
    elseif((wavelength >= 510.0) && (wavelength <= 580.0))
        R = (wavelength-510.0)/(580.0-510.0);
        G = 1.0;
        B = 0.0;
    elseif ((wavelength >= 580.0) && (wavelength <= 645.0))
        R = 1.0;
        G = -1.0*(wavelength-645.0)/(645.0-580.0);
        B = 0.0;
    elseif ((wavelength >= 645.0) && (wavelength <= 780.0))
        R = 1.0;
        G = 0.0;
        B = 0.0;
    else
        R = 0.0;
        G = 0.0;
        B = 0.0;
    end
end
