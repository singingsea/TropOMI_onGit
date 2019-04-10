function [secs] = load_leapsecs(mjd)
% LOAD_LEAPSECS Number of leap seconds.
%   function [secs] = load_leapsecs(mjd) reads in the number of leap
%   seconds for the given date.  The data are from the International Earth
%   Rotation Service (IERS).

%----------------------------------------------------------------
% Craig Haley 19/06/03
%  02-11-03 CSH made small change to remove for loop
%  29-06-04 CSH fixed for non-vectorized utc2mjd
%  02-10-05 CSH vectorized
%----------------------------------------------------------------

%assume the file is in the path
filename = 'leapseconds.txt';

%open the file
fid = fopen(filename,'rt');
if (fid == -1)
    error(['Can''t read file: ',filename]);
end

%read in the header info
jnk = fgetl(fid);
jnk = fgetl(fid);

%read in the table
[data] = fscanf(fid,'%lf',[4 inf]);
data = data';
fclose(fid);

%find the correct value
leap_secs = data(:,4);
leap_mjd = utc2mjd(data(:,3),data(:,2),data(:,1));  %at 00:00:00
for i=1:length(mjd)
    pts = find(leap_mjd <= mjd(i));
    secs = leap_secs(pts(length(pts)));
end
