clc
clear
format long g

alpha = [.3820e-7   .1490e-7  -.1790e-6   .0000];
beta = [.1430e6   .0000  -.3280e6   .1130e6];

% Time (UTC)
year = 2000;
month = 1;
day = 1;
hour = 20;
MIN = 45;
sec = 0;  
elev = 20;                        % Elevation (deg)
azimuth = 210;                    % Azimuth   (deg)
fi = 40;                          % Latitude (deg)
lambda = 260;                     % Longitude (deg)

% Compute Time-of-Week in seconds (UTC)
[week,tow] = cal2gpstime(year,month,day,hour,MIN,sec)

% Ionospheric L1 corr.
dIon1 = klobuchar(fi,lambda,elev,azimuth,tow,alpha,beta)

