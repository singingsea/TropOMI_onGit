function [mjd] = doy2mjd(doy,year)
% DOY2MJD Converts day of year to MJD.
%   [mjd] = doy2mjd(doy,year) finds the mjd for the given day of year
%   (1-365) and year.

mjd = utc2mjd(ones(size(doy)),ones(size(doy)),year)+doy-1;
