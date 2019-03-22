function [hh,mm,ss] = fracDay2time(fracday)
hh = floor(fracday*24);
mm = floor((fracday-hh/24)*24*60);
ss = floor((fracday-hh/24 - mm/24/60)*24*60*60);
end