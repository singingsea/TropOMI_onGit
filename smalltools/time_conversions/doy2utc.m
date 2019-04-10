function [month,day] = doy2utc(doy,year)
% DOY2UTC Converts day of year to MM and DD.
%   [month,day] = doy2utc(doy,year) finds the month and day for the given 
%   day of year (1-365) and year.

mjd = doy2mjd(doy,year);

[year_tmp,month,day] = mjd2utc(mjd);
