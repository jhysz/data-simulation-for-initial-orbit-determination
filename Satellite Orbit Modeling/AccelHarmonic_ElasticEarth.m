%--------------------------------------------------------------------------
%
% AccelHarmonic: Computes the acceleration due to the harmonic gravity
%                field of the central body (Solid Earth Tides and Ocean
%                Tides effects are considered)
%
% Inputs:
%   r_Sun       Geocentric equatorial position (in [m]) referred to the
%               mean equator and equinox of J2000 (EME2000, ICRF)
%   r_Moon      Geocentric equatorial position (in [m]) referred to the
%               mean equator and equinox of J2000 (EME2000, ICRF)
%   Mjd_UTC     Modified Julian Date of UTC
%   r           Satellite position vector in the inertial system
%   E           Transformation matrix to body-fixed system
%   Cnm,Snm     Spherical harmonics coefficients (normalized)
%   n_max       Maximum degree
%   m_max       Maximum order (m_max<=n_max; m_max=0 for zonals, only)
%   UT1_UTC     UT1-UTC time difference [s]
%   TT_UTC      TT-UTC time difference [s]
%   x_pole      X-Pole [rad]
%   y_pole      Y-Pole [rad]
%
% Output:
%   a           Acceleration (a=d^2r/dt^2)
%
% Last modified:   2020/03/24   Meysam Mahooti
% 
%--------------------------------------------------------------------------
function [a] = AccelHarmonic_ElasticEarth(Mjd_UTC,r_Sun,r_Moon,r,E,UT1_UTC,TT_UTC,x_pole,y_pole)

global Cnm Snm AuxParam const

r_ref = 6378.1363e3;   % Earth's radius [m]; GGM03S
gm    = 398600.4415e9; % [m^3/s^2]; GGM03S

C = Cnm;
S = Snm;

r_Moon = E*r_Moon;
[lM, phiM, rM] = CalcPolarAngles(r_Moon);
r_Sun = E*r_Sun;
[lS, phiS, rS] = CalcPolarAngles(r_Sun);

Mjd_TT = Mjd_UTC + TT_UTC/86400;
Mjd_UT1 = Mjd_UTC + UT1_UTC/86400;

T  = (Mjd_TT-const.MJD_J2000)/36525;
T2 = T*T;
T3 = T2*T;
rev = 360*3600;  % arcsec/revolution

if (AuxParam.SolidEarthTides)
    % Effect of Solid Earth Tides (elastic Earth)
    % For dC21 and dS21
    % The coefficients we choose are in-phase(ip) amplitudes and out-of-phase amplitudes of the
    % corrections for frequency dependence, and multipliers of the Delaunay variables
    % Refer to Table 6.5a in IERS2010
    coeff0 = [...
	%  l   l'  F   D   Om  Amp(R) Amp(I)
       2,  0,  2,  0,  2,  -0.1,    0;    
       0,  0,  2,  2,  2,  -0.1,    0;    
       1,  0,  2,  0,  1,  -0.1,    0;    
       1,  0,  2,  0,  2,  -0.7,    0.1;  
      -1,  0,  2,  2,  2,  -0.1,    0;    
       0,  0,  2,  0,  1,  -1.3,    0.1;  
       0,  0,  2,  0,  2,  -6.8,    0.6;  
       0,  0,  0,  2,  0,   0.1,    0;    
       1,  0,  2, -2,  2,   0.1,    0;    
      -1,  0,  2,  0,  1,   0.1,    0;    
      -1,  0,  2,  0,  2,   0.4,    0;    
       1,  0,  0,  0,  0,   1.3,   -0.1;  
       1,  0,  0,  0,  1,   0.3,    0;    
      -1,  0,  0,  2,  0,   0.3,    0;    
      -1,  0,  0,  2,  1,   0.1,    0;    
       0,  1,  2, -2,  2,  -1.9,    0.1;  
       0,  0,  2, -2,  1,   0.5,    0;    
       0,  0,  2, -2,  2,  -43.4,   2.9;
       0, -1,  2, -2,  2,   0.6,    0;    
       0,  1,  0,  0,  0,   1.6,   -0.1;  
      -2,  0,  2,  0,  1,   0.1,    0;    
       0,  0,  0,  0, -2,   0.1,    0;    
       0,  0,  0,  0, -1,  -8.8,    0.5;  
       0,  0,  0,  0,  0,   470.9, -30.2; 
       0,  0,  0,  0,  1,   68.1,  -4.6;  
       0,  0,  0,  0,  2,  -1.6,    0.1;  
      -1,  0,  0,  1,  0,   0.1,    0;    
       0, -1,  0,  0, -1,  -0.1,    0;    
       0, -1,  0,  0,  0,  -20.6,  -0.3;  
       0,  1, -2,  2, -2,   0.3,    0;    
       0, -1,  0,  0,  1,  -0.3,    0;    
      -2,  0,  0,  2,  0,  -0.2,    0;    
      -2,  0,  0,  2,  1,  -0.1,    0;    
       0,  0, -2,  2, -2,  -5.0,    0.3;  
       0,  0, -2,  2, -1,   0.2,    0;    
       0, -1, -2,  2, -2,  -0.2,    0;    
       1,  0,  0, -2,  0,  -0.5,    0;    
       1,  0,  0, -2,  1,  -0.1,    0;    
      -1,  0,  0,  0, -1,   0.1,    0;    
      -1,  0,  0,  0,  0,  -2.1,    0.1;  
      -1,  0,  0,  0,  1,  -0.4,    0;    
       0,  0,  0, -2,  0,  -0.2,    0;    
      -2,  0,  0,  0,  0,  -0.1,    0;    
       0,  0, -2,  0, -2,  -0.6,    0;    
       0,  0, -2,  0, -1,  -0.4,    0;    
       0,  0, -2,  0,  0,  -0.1,    0;    
      -1,  0, -2,  0, -2,  -0.1,    0;    
      -1,  0, -2,  0, -1,  -0.1,    0;    
    ];
    % For dC22 and dS22
    % Refer to Table 6.5c in IERS2010
    coeff2 = [...
    % l  l' F  D  Om   Amp
      1, 0, 2, 0, 2,  -0.3;
      0, 0, 2, 0, 2,  -1.2;
    ];
    
    % STEP1 CORRECTIONS
    [lgM, dlgM] = Legendre(4,4,phiM);
    [lgS, dlgS] = Legendre(4,4,phiS);
    dCnm20 = (0.29525/5)*( (const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,1)...
           + (const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,1) );
    dCnm21 = (0.29470/5)*( (const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,2)*cos(lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,2)*cos(lS) );
    dSnm21 = (0.29470/5)*( (const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,2)*(sin(lM))...
           + (const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,2)*(sin(lS)) );
    dCnm22 = (0.29801/5)*( (const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,3)*cos(2*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,3)*cos(2*lS) );
    dSnm22 = (0.29801/5)*( (const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,3)*(sin(2*lM))...
           + (const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,3)*(sin(2*lS)) );
    dCnm30 = (0.093/7)*( (const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,1)...
           + (const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,1) );
    dCnm31 = (0.093/7)*( (const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,2)*cos(lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,2)*cos(lS) );
    dSnm31 = (0.093/7)*( (const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,2)*sin(lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,2)*sin(lS) );
    dCnm32 = (0.093/7)*( (const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,3)*cos(2*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,3)*cos(2*lS) );
    dSnm32 = (0.093/7)*( (const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,3)*sin(2*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,3)*sin(2*lS) );
    dCnm33 = (0.094/7)*( (const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,4)*cos(3*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,4)*cos(3*lS) );
    dSnm33 = (0.094/7)*( (const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,4)*sin(3*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,4)*sin(3*lS) );
    dCnm40 = (-0.00087/5)*( (const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(5,1)...
           + (const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(5,1) );
    dCnm41 = (-0.00079/5)*( (const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(5,2)*cos(lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(5,2)*cos(lS) );
    dSnm41 = (-0.00079/5)*( (const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(5,2)*sin(lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(5,2)*sin(lS) );
    dCnm42 = (-0.00057/5)*( (const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(5,3)*cos(2*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(5,3)*cos(2*lS) );
    dSnm42 = (-0.00057/5)*( (const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(5,3)*sin(2*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(5,3)*sin(2*lS) );
    
    % STEP2 CORRECTIONS
    theta_g = gmst(Mjd_UT1);
    dC21 = 0;
    dS21 = 0;
    for i=1:48
        dC21 = dC21 + 1e-12*coeff0(i,6)*sin(theta_g+pi);
        dS21 = dS21 + 1e-12*coeff0(i,6)*cos(theta_g+pi);
    end
    dCnm21 = dCnm21 + dC21;
    dSnm21 = dSnm21 + dS21;
    
    dC22 = 0;
    dS22 = 0;
    for i=1:2
		dC22 = dC22 + 1e-12*coeff2(i,6)*sin(theta_g+pi);
        dS22 = dS22 + 1e-12*coeff2(i,6)*cos(theta_g+pi);
    end
    dCnm22 = dCnm22 + dC22;
    dSnm22 = dSnm22 + dS22;
    
    % Treatment of the Permanent Tide (elastic Earth)
    dC20 = 4.4228e-8*(-0.31460)*0.29525;
    dCnm20 = dCnm20 - dC20;
    
    % Effect of Solid Earth Pole Tide (elastic Earth)
    dC21 = -1.290e-9*(x_pole);
    dS21 = 1.290e-9*(y_pole);
    dCnm21 = dCnm21 + dC21;
    dSnm21 = dSnm21 + dS21;
    
    C(3,1) = C(3,1) + dCnm20;
    C(3,2) = C(3,2) + dCnm21;
    C(3,3) = C(3,3) + dCnm22;
    S(3,2) = S(3,2) + dSnm21;
    S(3,3) = S(3,3) + dSnm22;
    
    C(4,1) = C(4,1) + dCnm30;
    C(4,2) = C(4,2) + dCnm31;
    C(4,3) = C(4,3) + dCnm32;
    C(4,4) = C(4,4) + dCnm33;
    S(4,2) = S(4,2) + dSnm31;
    S(4,3) = S(4,3) + dSnm32;
    S(4,4) = S(4,4) + dSnm33;
    
    C(5,1) = C(5,1) + dCnm40;
    C(5,2) = C(5,2) + dCnm41;
    C(5,3) = C(5,3) + dCnm42;
    S(5,2) = S(5,2) + dSnm41;
    S(5,3) = S(5,3) + dSnm42;    
end

if (AuxParam.OceanTides)
    % Ocean Tides
    [lgM, dlgM] = Legendre(6,6,phiM);
    [lgS, dlgS] = Legendre(6,6,phiS);
    
    dCnm20 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.3075)/5*( (const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,1)...
           + (const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,1) );
    dCnm21 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.3075)/5*( (const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,2)*cos(lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,2)*cos(lS) );
    dSnm21 = -0.3075/5*( (const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,2)*sin(lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,2)*sin(lS) );
    dCnm22 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.3075)/5*( (const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,3)*cos(2*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,3)*cos(2*lS) );
    dSnm22 = -0.3075/5*( (const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,3)*sin(2*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,3)*sin(2*lS) );
    dCnm30 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.195)/7*( (const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,1)...
           + (const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,1) );
    dCnm31 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.195)/7*( (const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,2)*cos(lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,2)*cos(lS) );
    dSnm31 = -0.195/7*( (const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,2)*sin(lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,2)*sin(lS) );
    dCnm32 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.195)/7*( (const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,3)*cos(2*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,3)*cos(2*lS) );
    dSnm32 = -0.195/7*( (const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,3)*sin(2*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,3)*sin(2*lS) );
    dCnm33 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.195)/7*( (const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,4)*cos(3*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,4)*cos(3*lS) );
    dSnm33 = -0.195/7*( (const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,4)*sin(3*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,4)*sin(3*lS) );
    dCnm40 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.132)/9*( (const.GM_Moon/gm)*((r_ref/rM)^5)*lgM(5,1)...
           + (const.GM_Sun/gm)*((r_ref/rS)^5)*lgS(5,1) );
    dCnm41 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.132)/9*( (const.GM_Moon/gm)*((r_ref/rM)^5)*lgM(5,2)*cos(lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^5)*lgS(5,2)*cos(lS) );
    dSnm41 = -0.132/9*( (const.GM_Moon/gm)*((r_ref/rM)^5)*lgM(5,2)*sin(lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^5)*lgS(5,2)*sin(lS) );
    dCnm42 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.132)/9*( (const.GM_Moon/gm)*((r_ref/rM)^5)*lgM(5,3)*cos(2*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^5)*lgS(5,3)*cos(2*lS) );
    dSnm42 = -0.132/9*( (const.GM_Moon/gm)*((r_ref/rM)^5)*lgM(5,3)*sin(2*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^5)*lgS(5,3)*sin(2*lS) );
    dCnm43 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.132)/9*( (const.GM_Moon/gm)*((r_ref/rM)^5)*lgM(5,4)*cos(3*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^5)*lgS(5,4)*cos(3*lS) );
    dSnm43 = -0.132/9*( (const.GM_Moon/gm)*((r_ref/rM)^5)*lgM(5,4)*sin(3*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^5)*lgS(5,4)*sin(3*lS) );
    dCnm44 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.132)/9*( (const.GM_Moon/gm)*((r_ref/rM)^5)*lgM(5,5)*cos(4*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^5)*lgS(5,5)*cos(4*lS) );
    dSnm44 = -0.132/9*( (const.GM_Moon/gm)*((r_ref/rM)^5)*lgM(5,5)*sin(4*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^5)*lgS(5,5)*sin(4*lS) );
    dCnm50 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (const.GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,1)...
           + (const.GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,1) );
    dCnm51 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (const.GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,2)*cos(lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,2)*cos(lS) );
    dSnm51 = -0.1032/9*( (const.GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,2)*sin(lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,2)*sin(lS) );
    dCnm52 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (const.GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,3)*cos(2*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,3)*cos(2*lS) );
    dSnm52 = -0.1032/9*( (const.GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,3)*sin(2*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,3)*sin(2*lS) );
    dCnm53 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (const.GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,4)*cos(3*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,4)*cos(3*lS) );
    dSnm53 = -0.1032/9*( (const.GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,4)*sin(3*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,4)*sin(3*lS) );
    dCnm54 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (const.GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,5)*cos(4*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,5)*cos(4*lS) );
    dSnm54 = -0.1032/9*( (const.GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,5)*sin(4*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,5)*sin(4*lS) );
    dCnm55 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (const.GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,6)*cos(5*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,6)*cos(5*lS) );
    dSnm55 = -0.1032/9*( (const.GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,6)*sin(5*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,6)*sin(5*lS) );
    dCnm60 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (const.GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,1)...
           + (const.GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,1) );
    dCnm61 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (const.GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,2)*cos(lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,2)*cos(lS) );
    dSnm61 = -0.0892/9*( (const.GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,2)*sin(lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,2)*sin(lS) );
    dCnm62 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (const.GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,3)*cos(2*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,3)*cos(2*lS) );
    dSnm62 = -0.0892/9*( (const.GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,3)*sin(2*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,3)*sin(2*lS) );
    dCnm63 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (const.GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,4)*cos(3*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,4)*cos(3*lS) );
    dSnm63 = -0.0892/9*( (const.GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,4)*sin(3*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,4)*sin(3*lS) );
    dCnm64 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (const.GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,5)*cos(4*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,5)*cos(4*lS) );
    dSnm64 = -0.0892/9*( (const.GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,5)*sin(4*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,5)*sin(4*lS) );
    dCnm65 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (const.GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,6)*cos(5*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,6)*cos(5*lS) );
    dSnm65 = -0.0892/9*( (const.GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,6)*sin(5*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,6)*sin(5*lS) );
    dCnm66 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (const.GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,7)*cos(6*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,7)*cos(6*lS) );
    dSnm66 = -0.0892/9*( (const.GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,7)*sin(6*lM)...
           + (const.GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,7)*sin(6*lS) );
    
    C(3,1) = C(3,1) + dCnm20;
    C(3,2) = C(3,2) + dCnm21;
    C(3,3) = C(3,3) + dCnm22;
    S(3,2) = S(3,2) + dSnm21;
    S(3,3) = S(3,3) + dSnm22;
    
    C(4,1) = C(4,1) + dCnm30;
    C(4,2) = C(4,2) + dCnm31;
    C(4,3) = C(4,3) + dCnm32;
    C(4,4) = C(4,4) + dCnm33;
    S(4,2) = S(4,2) + dSnm31;
    S(4,3) = S(4,3) + dSnm32;
    S(4,4) = S(4,4) + dSnm33;
    
    C(5,1) = C(5,1) + dCnm40;
    C(5,2) = C(5,2) + dCnm41;
    C(5,3) = C(5,3) + dCnm42;
    C(5,4) = C(5,4) + dCnm43;
    C(5,5) = C(5,5) + dCnm44;
    S(5,2) = S(5,2) + dSnm41;
    S(5,3) = S(5,3) + dSnm42;
    S(5,4) = S(5,4) + dSnm43;
    S(5,5) = S(5,5) + dSnm44;
    
    C(6,1) = C(6,1) + dCnm50;
    C(6,2) = C(6,2) + dCnm51;
    C(6,3) = C(6,3) + dCnm52;
    C(6,4) = C(6,4) + dCnm53;
    C(6,5) = C(6,5) + dCnm54;
    C(6,6) = C(6,6) + dCnm55;
    S(6,2) = S(6,2) + dSnm51;
    S(6,3) = S(6,3) + dSnm52;
    S(6,4) = S(6,4) + dSnm53;
    S(6,5) = S(6,5) + dSnm54;
    S(6,6) = S(6,6) + dSnm55;
    
    C(7,1) = C(7,1) + dCnm60;
    C(7,2) = C(7,2) + dCnm61;
    C(7,3) = C(7,3) + dCnm62;
    C(7,4) = C(7,4) + dCnm63;
    C(7,5) = C(7,5) + dCnm64;
    C(7,6) = C(7,6) + dCnm65;
    C(7,7) = C(7,7) + dCnm66;
    S(7,2) = S(7,2) + dSnm61;
    S(7,3) = S(7,3) + dSnm62;
    S(7,4) = S(7,4) + dSnm63;
    S(7,5) = S(7,5) + dSnm64;
    S(7,6) = S(7,6) + dSnm65;
    S(7,7) = S(7,7) + dSnm66;    
end

% Body-fixed position 
r_bf = E * r;

% Auxiliary quantities
d = norm(r_bf);                     % distance
latgc = asin(r_bf(3)/d);
lon = atan2(r_bf(2),r_bf(1));

[pnm, dpnm] = Legendre(AuxParam.n,AuxParam.m,latgc);

dUdr = 0;
dUdlatgc = 0;
dUdlon = 0;
q3 = 0; q2 = q3; q1 = q2;
for n=0:AuxParam.n
    b1 = (-gm/d^2)*(r_ref/d)^n*(n+1);
    b2 =  (gm/d)*(r_ref/d)^n;
    b3 =  (gm/d)*(r_ref/d)^n;
    for m=0:AuxParam.m
        q1 = q1 + pnm(n+1,m+1)*(C(n+1,m+1)*cos(m*lon)+S(n+1,m+1)*sin(m*lon));
        q2 = q2 + dpnm(n+1,m+1)*(C(n+1,m+1)*cos(m*lon)+S(n+1,m+1)*sin(m*lon));
        q3 = q3 + m*pnm(n+1,m+1)*(S(n+1,m+1)*cos(m*lon)-C(n+1,m+1)*sin(m*lon));
    end
    dUdr     = dUdr     + q1*b1;
    dUdlatgc = dUdlatgc + q2*b2;
    dUdlon   = dUdlon   + q3*b3;
    q3 = 0; q2 = q3; q1 = q2;
end

% Body-fixed acceleration
r2xy = r_bf(1)^2+r_bf(2)^2;

ax = (1/d*dUdr-r_bf(3)/(d^2*sqrt(r2xy))*dUdlatgc)*r_bf(1)-(1/r2xy*dUdlon)*r_bf(2);
ay = (1/d*dUdr-r_bf(3)/(d^2*sqrt(r2xy))*dUdlatgc)*r_bf(2)+(1/r2xy*dUdlon)*r_bf(1);
az =  1/d*dUdr*r_bf(3)+sqrt(r2xy)/d^2*dUdlatgc;

a_bf = [ax ay az]';

% Inertial acceleration 
a = E'*a_bf;

