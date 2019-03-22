function f = date2fracyear(year,month,day,hour,minute,second)
jd = juliandate(year,month,day,hour,minute,second);
f = year  + jd/(365/leapyear(year));
end