function [mjd] = tai932mjd(tai)
% TAI2MJD Converts international atomic time to modified julian date.
%   [mjd] = tai2mjd(tai) converts the given international atomic time (TAI)
%   to modified julian date.

%----------------------------------------------------------------
% Chris McLinden, modified Craig Haley's original 06-02-07
%----------------------------------------------------------------

%get approximate mjd for determining correction
mjd = jd2mjd(tai/(24*60*60)+mjd2jd(utc2mjd(1,1,1993)));

%load in the correction
[correction] = load_leapsecs(mjd) - load_leapsecs(utc2mjd(1,1,1993));

%make the correction
secs = tai-correction;

%convert to days
jd = secs/(24*60*60);

%TAI93 seconds 0 is defined at 01 Jan 1993
jd = jd+mjd2jd(utc2mjd(1,1,1993));

%get the modified julian date
mjd = jd2mjd(jd);
