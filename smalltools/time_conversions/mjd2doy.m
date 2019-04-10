function [doy,year] = mjd2doy(mjd)
% MJD2DOY Converts modified julian date to day of year.
%   [doy,year] = mjd2doy(mjd) finds the day of year (1-365) for the given 
%   modified julian date.

[year,month,day] = mjd2utc(mjd);

doy = utc2mjd(day,month,year)-utc2mjd(ones(size(year)),ones(size(year)),year)+1;
