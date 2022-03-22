%--------------------------------------------------------------------------
%
% gast: Greenwich Apparent Sidereal Time
%
% Input:
%   Mjd_UT1   Modified Julian Date UT1
%
% Output:
%   gstime    GAST in [rad]
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------
function gstime = gast(Mjd_UT1)

global const

gstime = mod(gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), const.pi2);

