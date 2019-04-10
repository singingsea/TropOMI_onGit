function [year,month,day,hour,minute,secs,ticks] = mjd2utc(mjd)
% MJD2UTC Converts modified julian date to universal time coordinated.
%   [year,month,day,hour,minute,secs,ticks] = mjd2utc(mjd) finds the civil
%   calendar gate for the given modified julian date. Julian calendar is 
%   used up to 1582 October 4, Gregorian calendar is used from 1582
%   October 15 onwards. Follows Algorithm given in "Astronomy on the
%   Personal Computer"  O. Montenbruck and T. Pfleger.

%  19-03-05 CSH vectorized with simple for loop
%  02-10-05 CSH vectorized

mjd=double(mjd);

dayfrac = mjd-floor(mjd);           % Get the fraction of UTC day.
jd = mjd+2400000.5;                 % Get the Julian Date
jd0 = floor(jd+0.5);                % and add a half.

c = jd0+1524;                       % Julian
pts = find(jd0 >= 2299161);
b = fix(((jd0(pts)-1867216.25)/36524.25));
c(pts) = jd0+(b-fix(b/4))+1525.0;   % Gregorian
d = fix(((c-122.1)/365.25));
e = 365.0*d+fix(d/4);
f = fix(((c-e)/30.6001));

day = fix((c-e+0.5)-fix(30.6001*f));
month = fix((f-1-12*fix(f/14)));
year = fix((d-4715-fix((7+month)/10)));

if (nargout == 0) || (nargout > 3)
    second = dayfrac*86400;
    hour = fix(second/3600);
    second = second-3600*hour;
    minute = fix(second/60);
    second = second-60*minute;
    secs = fix(second);
    ticks = second-secs;
end

if nargout == 0
    year = [num2str(year(1)) '-' num2str(month(1),'%02d') '-' num2str(day(1),'%02d') ' '...
        num2str(hour(1)) ':' num2str(minute(1),'%02d') ':' num2str(secs(1)+ticks(1),'%05.3f')];
end
