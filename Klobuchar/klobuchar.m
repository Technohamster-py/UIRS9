% ==================================================================
% function for computing an Ionospheric range correction for the *
% GPS L1 frequency from the parameters broadcasted in the GPS      *
% Navigation Message.                                              *
% ==================================================================
% References:                                                      *
% Klobuchar, J.A., (1996) "Ionosphercic Effects on GPS", in        *
%   Parkinson, Spilker (ed), "Global Positioning System Theory and *
%   Applications, pp.513-514.                                      *
% ICD-GPS-200, Rev. C, (1997), pp. 125-128                         *
% NATO, (1991), "Technical Characteristics of the NAVSTAR GPS",    *
%   pp. A-6-31   -   A-6-33                                        *
% ==================================================================
% Inputs:                                                          *
%   fi            : Geodetic latitude of receiver          (deg)   *
%   lambda        : Geodetic longitude of receiver         (deg)   *
%   elev          : Elevation angle of satellite           (deg)   *
%   azimuth       : Geodetic azimuth of satellite          (deg)   *
%   tow           : Time of Week                           (sec)   *
%   alfa(4)       : The coefficients of a cubic equation           *
%                   representing the amplitude of the vertical     *
%                   dalay (4 coefficients - 8 bits each)           *
%   beta(4)       : The coefficients of a cubic equation           *
%                   representing the period of the model           *
%                   (4 coefficients - 8 bits each)                 *
% Output:                                                          *
%   dIon1         : Ionospheric slant range correction for         *
%                   the L1 frequency                       (metre) *
% ==================================================================
function dIon1 = klobuchar(fi,lambda,elev,azimuth,tow,alfa,beta)

c        =  2.99792458e8;             % speed of light
deg2semi =  1./180.;                  % degees to semisircles
semi2rad =  pi;                       % semisircles to radians
deg2rad  =  pi/180.;                  % degrees to radians

a = azimuth*deg2rad;                  % asimuth in radians
e = elev*deg2semi;                    % elevation angle in
                                      % semicircles

psi = 0.0137 / (e+0.11) - 0.022;      % Earth Centered angle

lat_i = fi*deg2semi + psi*cos(a);     % Subionospheric lat
if (lat_i > 0.416)
    lat_i = 0.416;
  elseif(lat_i < -0.416)
      lat_i = -0.416;
end

                                      % Subionospheric long
long_i = lambda*deg2semi + (psi*sin(a)/cos(lat_i*semi2rad));

                                      % Geomagnetic latitude
lat_m = lat_i + 0.064*cos((long_i-1.617)*semi2rad);

t = 4.32e4*long_i + tow;
t = mod(t,86400.);                    % Seconds of day
if (t > 86400.)
    t = t - 86400.;
end
if (t < 0.)
    t = t + 86400.;
end

sF = 1. + 16.*(0.53-e)^3;             % Slant factor

                                      % Period of model
PER = beta(1) + beta(2)*lat_m + beta(3)*lat_m^2 +beta(4)*lat_m^3;

if (PER < 72000.)
    PER = 72000.;
end

x = 2.*pi*(t-50400.)/PER ;            % Phase of the model
                                      % (Max at 14.00 =
                                      % 50400 sec local time)

                                      % Amplitud of the model
AMP = alfa(1) + alfa(2)*lat_m + alfa(3)*lat_m^2 +alfa(4)*lat_m^3;
if(AMP < 0.)
    AMP = 0.;
end

                                      % Ionospheric corr.
if(abs(x) > 1.57)
    dIon1 = sF * (5.e-9);
else
    dIon1 = sF * (5.e-9 + AMP*(1. - x*x/2. + x*x*x*x/24.));
end

dIon1 = c * dIon1;

