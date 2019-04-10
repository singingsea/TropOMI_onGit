function [mjd] = utc2mjd(day,month,year,hour,minute,secs,ticks)
% UTC2MJD Convert universal time coordinated to modified julian date.
%   function [mjd] = utc2mjd(day,month,year,hour,minute,secs,ticks) finds 
%   the modified julian date from the given civil calendar date. The
%   Julian calendar is used up to 1582 October 4 and Gregorian calendar 
%   is used from 1582 October 15 onwards. Follows Algorithm given in 
%   "Astronomy on the Personal Computer"  O. Montenbruck and T. Pfleger.

% 19-03-05 CSH vectorized with simple for loop
% 07-09-05 CSH modified to not require definition of hh,mm,ss,ticks
% 02-10-05 CSH vectorized

switch nargin
    case {1,2}
        error('Day, Month and Year required.')
    case 3
        hour = zeros(size(day));
        minute = zeros(size(day));
        secs = zeros(size(day));
        ticks = zeros(size(day));
    case 4
        minute = zeros(size(day));
        secs = zeros(size(day));
        ticks = zeros(size(day));
    case 5
        secs = zeros(size(day));
        ticks = zeros(size(day));
    case 6
        ticks = zeros(size(day));
    case 7
        %nothing needs to be done
    otherwise
        error('Too many input arguments.');
end

year(year <= 40) = year(year <= 40)+2000;
year((year > 40) & (year < 100)) = year((year > 40) & (year < 100))+1900;

y = fix(year);

pts = find(month <= 2);
y(pts) = y(pts)-1;
month(pts) = month(pts)+12;

B = fix(y/400)-fix(y/100)+fix(y/4);
pts = find((year < 1582) | ((year == 1582) & ((month < 10) | ((month == 10) & (day < 15)))));
B(pts) = fix(-2+fix((y(pts)+4716)/4)-1179);

A = 365.0*y-679004.0;
mjd = A+B+(fix(30.6001*(month+1)))+day+(hour/24.0)+minute/1440.0+(secs+ticks)/86400.0;
