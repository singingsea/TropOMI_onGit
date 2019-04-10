function [doy] = utc2doy(day,month,year)
% UTC2DOY Converts day, month and year to the day of year.
%   [doy] = utc2doy(day,month,year) finds the day of year (1-365) from the 
%   given day, month and year.

mjd = utc2mjd(day,month,year);

[doy] = mjd2doy(mjd);
