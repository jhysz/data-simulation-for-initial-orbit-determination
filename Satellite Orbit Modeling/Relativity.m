%--------------------------------------------------------------------------
%
% Relativisty: Computes the perturbational acceleration due to relativistic
%              effects
%
% Inputs:
%   r           Satellite position vector
%   v           Satellite velocity vector
% 
% Output:
%   a    		Acceleration (a=d^2r/dt^2)
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------
function a = Relativity(r, v)

global const

% Relative position vector of satellite w.r.t. point mass 
r_Sat = norm(r);
v_Sat = norm(v);

% Acceleration 
a = const.GM_Earth/(const.c_light^2*r_Sat^3)*((4*const.GM_Earth/r_Sat-v_Sat^2)*r+4*dot(r,v)*v);

