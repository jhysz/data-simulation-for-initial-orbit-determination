%--------------------------------------------------------------------------
%
% AccelSolrad: Computes the acceleration due to solar radiation pressure
%              assuming the spacecraft surface normal to the Sun direction
%
% Inputs:
%   r           Spacecraft position vector 
%   r_Sun       Sun position vector 
%   Area        Cross-section 
%   mass        Spacecraft mass
%   CR          Solar radiation pressure coefficient
%   P0          Solar radiation pressure at 1 AU 
%   AU          Length of one Astronomical Unit
%
% Output:
%   a    		Acceleration (a=d^2r/dt^2)
%
% Notes:
%   r, r_sun, Area, mass, P0 and AU must be given in consistent units,
%   e.g. m, m^2, kg and N/m^2. 
%
% Last modified:   2018/01/27   M. Mahooti
% 
%--------------------------------------------------------------------------
function a = AccelSolrad(r,r_Sun,Area,mass,Cr,P0,AU)

nu = Cylindrical(r, r_Sun);

% Relative position vector of spacecraft w.r.t. Sun
d = r - r_Sun;

% Acceleration
a = nu*Cr*(Area/mass)*P0*(AU*AU) * d/(norm(d)^3);

