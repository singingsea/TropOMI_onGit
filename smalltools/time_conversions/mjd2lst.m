function [lst] = mjd2lst(mjd,longitude)
% MJD2LST Convert modified julian date to local solar time.
%   [lst] = mjd2lst(mjd,longitude) converts the given modified julian date
%   to local solar time (in hours) at the given east longitude (degrees).

%----------------------------------------------------------------
% Craig Haley 07/03/02
%----------------------------------------------------------------

[year,month,day,hour,minute,secs,ticks] = mjd2utc(mjd);

lst = mod(hour+minute/60+(secs+ticks)/60/60+longitude/15.0,24.0);
