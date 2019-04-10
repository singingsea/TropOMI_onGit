function [mjd] = jd2mjd(jd)
% JD2MJD Convert Julian Date to Modified Julian Date
%   [mjd] = jd2mjd(jd) converts the given Julian Date to the Modified
%   Julian Date.

%----------------------------------------------------------------
% Craig Haley 09-01-07
%----------------------------------------------------------------

mjd = jd-2400000.5;
