%--------------------------------------------------------------------------
%
% nrlmsise00: Computes the atmospheric density for the nrlmsise00 model
%
% Inputs:
%   Mjd_UTC     Modified Julian Date UTC
%   r_ECEF      Satellite position vector in ECEF coordinate system [m]
%   UT1_UTC     UT1-UTC time difference [s]
%   TT_UTC      TT-UTC time difference [s]
% 
% Output:
%   dens        Density [kg/m^3]
%
% Last modified:   2018/02/11   M. Mahooti
% 
%--------------------------------------------------------------------------
function dens = nrlmsise00(Mjd_UTC, r_ECEF, UT1_UTC, TT_UTC)

global const swdata
global pt pd ps pdl ptl pma sam ptm pdm pavgm
global gsurf re dd dm04 dm16 dm28 dm32 dm40 dm01 dm14
global meso_tn1 meso_tn2 meso_tn3 meso_tgn1 meso_tgn2 meso_tgn3
global dfa plg ctloc stloc c2tloc s2tloc s3tloc c3tloc apdf apt
nrlmsise00_coeff

gsurf = 0;
re = 6367.088132098377173;
dd = 0;
dm04 = 0;
dm16 = 0;
dm28 = 0;
dm32 = 0;
dm40 = 0;
dm01 = 0;
dm14 = 0;
meso_tn1 = zeros(5,1);
meso_tn2 = zeros(4,1);
meso_tn3 = zeros(5,1);
meso_tgn1 = zeros(2,1);
meso_tgn2 = zeros(2,1);
meso_tgn3 = zeros(2,1);
dfa = 0;
plg = zeros(4,9);
ctloc = 0;
stloc = 0;
c2tloc = 0;
s2tloc = 0;
s3tloc = 0;
c3tloc = 0;
apdf = 0;
apt = zeros(4,1);

flags = struct('switches',zeros(24,1),'sw',zeros(24,1),'swc',zeros(24,1));
ap = struct('a',zeros(7,1));
input = struct('year',0,'doy',0,'sec',0,'alt',0,'g_lat',0,'g_long',0,'lst',0,'f107A',0,'f107',0,'ap',0,'ap_a',ap);

flags.switches(1)=0;
for i=2:24
    flags.switches(i)=1;
end
flags.switches(10)=-1;

[year, mon, day, hour, minute, sec] = invjday(Mjd_UTC+2400000.5);
days = finddays(year, mon, day, hour, minute, sec);
input.doy = floor(days);
input.year = 0; % without effect
input.sec = hour*3600+minute*60+sec; % seconds in day (UT)

[lon, lat, height] = Geodetic(r_ECEF);
input.alt = height/1000;
input.g_lat = lat*const.Deg;
input.g_long = lon*const.Deg;

Mjd_UT1 = Mjd_UTC + UT1_UTC/86400;
lst = lon + gast(Mjd_UT1);
lst = mod(lst,const.pi2);
lst = (lst*24)/(const.pi2); % hours
input.lst=lst;

i = find((year==swdata(1,:)) & (mon==swdata(2,:)) & (day==swdata(3,:)),1,'first');
sw = swdata(:,i);
input.ap = sw(23);        % Arithmetic average of the 8 Ap indices for the day
input.ap_a.a(1) = sw(23); % daily AP
input.ap_a.a(2) = sw(15); % 3 hr AP index for current time
sw_1 = swdata(:,i-1);
% Define Solar Flux Values
input.ap_a.a(3) = sw_1(22); % 3 hr AP index for 3 hrs before current time
input.ap_a.a(4) = sw_1(21); % 3 hr AP index for 6 hrs before current time
input.ap_a.a(5) = sw_1(20); % 3 hr AP index for 9 hrs before current time
sum  = sw_1(19)+sw_1(18)+sw_1(17)+sw_1(16)+sw_1(15);
sw_2 = swdata(:,i-2);
sum = sum+sw_2(22)+sw_2(21)+sw_2(20);
input.ap_a.a(6) = sum/8; % Average of eight 3 hr AP indicies from 12 to 33
                         % hrs prior to current time
sw_3 = swdata(:,i-3);
sum = sw_2(19)+sw_2(18)+sw_2(17)+sw_2(16)+sw_2(15)+sw_3(22)+sw_3(21)+sw_3(20);
input.ap_a.a(7) = sum/8; % Average of eight 3 hr AP indicies from 36 to 57
                         % hrs prior to current time 

input.f107 = sw(31);     % solar radio noise flux (jansky)
input.f107A = sw(33);    % 162-day average F10 (jansky)

output = gtd7d(input, flags);
dens = 1e3*output.d(6);

end


% ------------------------------------------------------------------- %
% ------------------------------ TSELEC ----------------------------- %
% ------------------------------------------------------------------- %
function flags = tselec(flags)

for i=1:24
    if (i~=10)
        if (flags.switches(i)==1)
            flags.sw(i)=1;
        else
            flags.sw(i)=0;
        end        
        if (flags.switches(i)>0)
            flags.swc(i)=1;
        else
            flags.swc(i)=0;
        end
    else
        flags.sw(i)=flags.switches(i);
        flags.swc(i)=flags.switches(i);
    end
end

end

% ------------------------------------------------------------------- %
% ------------------------------ GLATF ------------------------------ %
% ------------------------------------------------------------------- %
function [gv, reff] = glatf(lat)

dgtr = 1.74533E-2;
c2 = cos(2*dgtr*lat);
gv = 980.616 * (1 - 0.0026373 * c2);
reff = 2 * (gv) / (3.085462E-6 + 2.27E-9 * c2) * 1E-5;

end

% ------------------------------------------------------------------- %
% ------------------------------ CCOR ------------------------------- %
% ------------------------------------------------------------------- %
function out = ccor(alt, r, h1, zh)

%        CHEMISTRY/DISSOCIATION CORRECTION FOR MSIS MODELS
%        ALT - altitude
%        R - target ratio
%        H1 - transition scale length
%        ZH - altitude of 1/2 R
e = (alt - zh) / h1;
if (e>70)
    out = exp(0);
    return
end
if (e<-70)
    out = exp(r);
    return
end
ex = exp(e);
e = r / (1 + ex);
out = exp(e);

end

% ------------------------------------------------------------------- %
% ------------------------------ CCOR ------------------------------- %
% ------------------------------------------------------------------- %
function out = ccor2(alt, r, h1, zh, h2)
%        CHEMISTRY/DISSOCIATION CORRECTION FOR MSIS MODELS
%        ALT - altitude
%        R - target ratio
%        H1 - transition scale length
%        ZH - altitude of 1/2 R
%        H2 - transition scale length #2 ?
e1 = (alt - zh) / h1;
e2 = (alt - zh) / h2;
if ((e1 > 70) || (e2 > 70))
    out = exp(0);
    return
end
if ((e1 < -70) && (e2 < -70))
    out = exp(r);
    return
end
ex1 = exp(e1);
ex2 = exp(e2);
ccor2v = r / (1 + 0.5 * (ex1 + ex2));
out = exp(ccor2v);

end

% ------------------------------------------------------------------- %
% ------------------------------- SCALH ----------------------------- %
% ------------------------------------------------------------------- %
function g = scalh(alt, xm, temp)

global gsurf re

rgas = 831.4;
g = gsurf / ((1 + alt/re)^2);
g = rgas * temp / (g * xm);

end

% ------------------------------------------------------------------- %
% -------------------------------- DNET ----------------------------- %
% ------------------------------------------------------------------- %
function out = dnet (dd, dm, zhm, xmm, xm)
%       TURBOPAUSE CORRECTION FOR MSIS MODELS
%       Root mean density
%       DD - diffusive density
%       DM - full mixed density
%       ZHM - transition scale length
%       XMM - full mixed molecular weight
%       XM  - species molecular weight
%       DNET - combined density
a  = zhm / (xmm-xm);
if (~((dm>0) && (dd>0)))
    fprintf('dnet log error %e %e %e\n',dm,dd,xm);
    if ((dd==0) && (dm==0))
        dd=1;
    end
    if (dm==0)
        out = dd;
        return
    end
    if (dd==0)
        out = dm;
        return
    end
end
ylog = a * log(dm/dd);
if (ylog<-10)
    out = dd;
    return
end
if (ylog>10)
    out = dm;
    return
end
a = dd*(1 + exp(ylog))^(1/a);
out = a;

end


% ------------------------------------------------------------------- %
% ------------------------------- SPLINI ---------------------------- %
% ------------------------------------------------------------------- %
function y = splini(xa, ya, y2a, n, x)
%      INTEGRATE CUBIC SPLINE FUNCTION FROM XA(1) TO X
%      XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
%      Y2A: ARRAY OF SECOND DERIVATIVES
%      N: SIZE OF ARRAYS XA,YA,Y2A
%      X: ABSCISSA ENDPOINT FOR INTEGRATION
%      Y: OUTPUT VALUE
yi=0;
klo=0;
khi=1;
while ((x>xa(klo+1)) && (khi<n))
    xx=x;
    if (khi<(n-1))
        if (x<xa(khi+1))
            xx=x;
        else
            xx=xa(khi+1);
        end
    end
    h = xa(khi+1) - xa(klo+1);
    a = (xa(khi+1) - xx)/h;
    b = (xx - xa(klo+1))/h;
    a2 = a*a;
    b2 = b*b;
    yi = yi + ((1 - a2) * ya(klo+1) / 2 + b2 * ya(khi+1) / 2 + ((-(1+a2*a2)/4 + a2/2) * y2a(klo+1) + (b2*b2/4 - b2/2) * y2a(khi+1)) * h * h / 6) * h;
    klo= klo+1;
    khi= khi+1;
end
y = yi;

end


% ------------------------------------------------------------------- %
% ------------------------------- SPLINT ---------------------------- %
% ------------------------------------------------------------------- %
function y = splint (xa, ya, y2a, n, x)
%      CALCULATE CUBIC SPLINE INTERP VALUE
%      ADAPTED FROM NUMERICAL RECIPES BY PRESS ET AL.
%      XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
%      Y2A: ARRAY OF SECOND DERIVATIVES
%      N: SIZE OF ARRAYS XA,YA,Y2A
%      X: ABSCISSA FOR INTERPOLATION
%      Y: OUTPUT VALUE
klo=0;
khi=n-1;
while ((khi-klo)>1)
    k=int8((khi+klo)/2);
    if (xa(k+1)>x)
        khi=k;
    else
        klo=k;
    end
end
h = xa(khi+1) - xa(klo+1);
if (h==0)
    fprintf('bad XA input to splint');
end
a = (xa(khi+1) - x)/h;
b = (x - xa(klo+1))/h;
yi = a * ya(klo+1) + b * ya(khi+1) + ((a*a*a - a) * y2a(klo+1) + (b*b*b - b) * y2a(khi+1)) * h * h/6;
y = yi;

end


% ------------------------------------------------------------------- %
% ------------------------------- SPLINE ---------------------------- %
% ------------------------------------------------------------------- %
function y2 = spline (x, y, n, yp1, ypn)
%      CALCULATE 2ND DERIVATIVES OF CUBIC SPLINE INTERP FUNCTION
%      ADAPTED FROM NUMERICAL RECIPES BY PRESS ET AL
%      X,Y: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
%      N: SIZE OF ARRAYS X,Y
%      YP1,YPN: SPECIFIED DERIVATIVES AT X(0) AND X(N-1); VALUES
%               >= 1E30 SIGNAL SIGNAL SECOND DERIVATIVE ZERO
%      Y2: OUTPUT ARRAY OF SECOND DERIVATIVES
u = zeros(5,1);
y2 = zeros(5,1);

if (yp1>0.99E30)
    y2(1)=0;
    u(1)=0;
else
    y2(1)=-0.5;
    u(1)=(3/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1);
end
for i=2:n-1
	sig = (x(i)-x(i-1))/(x(i+1) - x(i-1));
	p = sig * y2(i-1) + 2;
	y2(i) = (sig - 1) / p;
	u(i) = (6 * ((y(i+1) - y(i))/(x(i+1) - x(i)) -(y(i) - y(i-1)) / (x(i) - x(i-1)))/(x(i+1) - x(i-1)) - sig * u(i-1))/p;
end
if (ypn>0.99E30)
    qn = 0;
    un = 0;
else
    qn = 0.5;
    un = (3 / (x(n) - x(n-1))) * (ypn - (y(n) - y(n-1))/(x(n) - x(n-1)));
end
y2(n) = (un - qn * u(n-1)) / (qn * y2(n-1) + 1);
for k=n-1:-1:1
    y2(k) = y2(k) * y2(k+1) + u(k);
end

end


% ------------------------------------------------------------------- %
% ------------------------------- DENSM ----------------------------- %
% ------------------------------------------------------------------- %
function [densm_tmp, tz] = densm (alt, d0, xm, tz, mn3, zn3, tn3, tgn3, mn2, zn2, tn2, tgn2)

global gsurf re

zeta = inline('(zz-zl)*(6367.088132098377173+zl)/(6367.088132098377173+zz)','zz','zl');

% Calculate Temperature and Density Profiles for lower atmos.
xs = zeros(10,1);
ys = zeros(10,1);
rgas = 831.4;
densm_tmp = d0;

if (alt>zn2(1))
    if (xm==0)
        densm_tmp = tz;
        return            
    else
        densm_tmp = d0;
        return
    end
end

% STRATOSPHERE/MESOSPHERE TEMPERATURE
if (alt>zn2(mn2))
    z = alt;
else
    z = zn2(mn2);
end
mn=mn2;
z1=zn2(1);
z2=zn2(mn);
t1=tn2(1);
t2=tn2(mn);
zg = zeta(z, z1);
zgdif = zeta(z2, z1);

% set up spline nodes
for k=1:mn
    xs(k)=zeta(zn2(k),z1)/zgdif;
    ys(k)=1 / tn2(k);
end

yd1 = -tgn2(1) / (t1*t1) * zgdif;
yd2 = -tgn2(2) / (t2*t2) * zgdif * (((re+z2)/(re+z1))^2);

% calculate spline coefficients
y2out = spline (xs, ys, mn, yd1, yd2);
x = zg/zgdif;
y = splint (xs, ys, y2out, mn, x);

% temperature at altitude
tz = 1 / y;
if (xm~=0)
    % calaculate stratosphere / mesospehere density
    glb = gsurf / ((1 + z1/re)^2);
    gamm = xm * glb * zgdif / rgas;
    
    % Integrate temperature profile
    yi = splini(xs, ys, y2out, mn, x);
    expl = gamm*yi;
    if (expl>50)
        expl = 50;
    end
    % Density at altitude
    densm_tmp = densm_tmp * (t1 / tz) * exp(-expl);
end

if (alt>zn3(1))
    if (xm==0)
        densm_tmp = tz;
        return        
    else
        return
    end
end

% troposhere / stratosphere temperature
z = alt;
mn = mn3;
z1 = zn3(1);
z2 = zn3(mn);
t1 = tn3(1);
t2 = tn3(mn);
zg = zeta(z,z1);
zgdif = zeta(z2,z1);

% set up spline nodes
for k=1:mn
    xs(k) = zeta(zn3(k),z1) / zgdif;
    ys(k) = 1 / tn3(k);
end
yd1 = -tgn3(1) / (t1*t1) * zgdif;
yd2 = -tgn3(2) / (t2*t2) * zgdif * (((re+z2)/(re+z1))^2);

% calculate spline coefficients
y2out =	spline (xs, ys, mn, yd1, yd2);
x = zg/zgdif;
y =	splint (xs, ys, y2out, mn, x);

% temperature at altitude
tz = 1 / y;
if (xm~=0)
    % calaculate tropospheric / stratosphere density
    glb = gsurf / ((1 + z1/re)^2);
    gamm = xm * glb * zgdif / rgas;
    
    % Integrate temperature profile
    yi = splini(xs, ys, y2out, mn, x);
    expl = gamm*yi;
    if (expl>50)
        expl = 50;
    end
    % Density at altitude
    densm_tmp = densm_tmp * (t1 / tz) * exp(-expl);
end
if (xm==0)
    densm_tmp = tz;
    return
else
    return
end

end


% ------------------------------------------------------------------- %
% ------------------------------- DENSU ----------------------------- %
% ------------------------------------------------------------------- %
function [densu_temp, tz] = densu (alt, dlb, tinf, tlb, xm, alpha, tz, zlb, s2, mn1, zn1, tn1, tgn1)
%      Calculate Temperature and Density Profiles for MSIS models
%      New lower thermo polynomial

global gsurf re

rgas = 831.4;
xs = zeros(5,1);
ys = zeros(5,1);
y2out = zeros(5,1);

% joining altitudes of Bates and spline
za=zn1(1);
if (alt>za)
    z=alt;
else
    z=za;
end
% geopotential altitude difference from ZLB
zg2 = zeta(z, zlb);

% Bates temperature
tt = tinf - (tinf - tlb) * exp(-s2*zg2);
ta = tt;
tz = tt;
densu_temp = tz;

if (alt<za)
    % calculate temperature below ZA
	% temperature gradient at ZA from Bates profile
	dta = (tinf - ta) * s2 * ((re+zlb)/(re+za))^2;
	tgn1(1)=dta;
	tn1(1)=ta;
	if (alt>zn1(mn1))
        z=alt;
    else
        z=zn1(mn1);
    end
	mn=mn1;
	z1=zn1(1);
	z2=zn1(mn);
	t1=tn1(1);
	t2=tn1(mn);
	% geopotental difference from z1
	zg = zeta (z, z1);
	zgdif = zeta(z2, z1);
    % set up spline nodes
	for k=1:mn
		xs(k) = zeta(zn1(k), z1) / zgdif;
		ys(k) = 1 / tn1(k);
    end
	% end node derivatives
	yd1 = -tgn1(1) / (t1*t1) * zgdif;
	yd2 = -tgn1(2) / (t2*t2) * zgdif * ((re+z2)/(re+z1))^2;
	% calculate spline coefficients
	y2out = spline (xs, ys, mn, yd1, yd2);
    x = zg / zgdif;
	y = splint (xs, ys, y2out, mn, x);
    % temperature at altitude
	tz = 1 / y;
	densu_temp = tz;
end
if (xm==0)
	return
end
% calculate density above za
glb = gsurf / (1 + zlb/re)^2;
gamma = xm * glb / (s2 * rgas * tinf);
expl = exp(-s2 * gamma * zg2);
if (expl>50)
	expl=50;
end
if (tt<=0)
	expl=50;
end
% density at altitude
densa = dlb * (tlb/tt)^(1+alpha+gamma) * expl;
densu_temp=densa;
if (alt>=za)
	return
end
% calculate density below za
glb = gsurf / (1 + z1/re)^2;
gamm = xm * glb * zgdif / rgas;

% integrate spline temperatures
yi = splini (xs, ys, y2out, mn, x);
expl = gamm * yi;
if (expl>50)
	expl=50;
end
if (tz<=0)
	expl=50;
end

% density at altitude
densu_temp = densu_temp * (t1 / tz)^(1 + alpha) * exp(-expl);

end

% ------------------------------------------------------------------- %
% ------------------------------- GLOBE7 ---------------------------- %
% ------------------------------------------------------------------- %
function tinf = globe7(p, input, flags)

global dfa plg ctloc stloc c2tloc s2tloc s3tloc c3tloc apdf apt

% CALCULATE G(L) FUNCTION 
% Upper Thermosphere Parameters
t = zeros(15,1);
sr = 7.2722E-5;
dgtr = 1.74533E-2;
dr = 1.72142E-2;
hr = 0.2618;

tloc = input.lst;

for j=1:14
    t(j) = 0;
end

% calculate legendre polynomials
c = sin(input.g_lat * dgtr);
s = cos(input.g_lat * dgtr);
c2 = c*c;
c4 = c2*c2;
s2 = s*s;

plg(1,2) = c;
plg(1,3) = 0.5*(3*c2 -1);
plg(1,4) = 0.5*(5*c*c2-3*c);
plg(1,5) = (35*c4 - 30*c2 + 3)/8;
plg(1,6) = (63*c2*c2*c - 70*c2*c + 15*c)/8;
plg(1,7) = (11*c*plg(1,6) - 5*plg(1,5))/6;
% plg(1,8) = (13*c*plg(1,7) - 6*plg(1,6))/7;
plg(2,2) = s;
plg(2,3) = 3*c*s;
plg(2,4) = 1.5*(5*c2-1)*s;
plg(2,5) = 2.5*(7*c2*c-3*c)*s;
plg(2,6) = 1.875*(21*c4 - 14*c2 +1)*s;
plg(2,7) = (11*c*plg(2,6)-6*plg(2,5))/5;
% plg(2,8) = (13*c*plg(2,7)-7*plg(2,6))/6;
% plg(2,9) = (15*c*plg(2,8)-8*plg(2,7))/7;
plg(3,3) = 3*s2;
plg(3,4) = 15*s2*c;
plg(3,5) = 7.5*(7*c2 -1)*s2;
plg(3,6) = 3*c*plg(3,5)-2*plg(3,4);
plg(3,7) =(11*c*plg(3,6)-7*plg(3,5))/4;
plg(3,8) =(13*c*plg(3,7)-8*plg(3,6))/5;
plg(4,4) = 15*s2*s;
plg(4,5) = 105*s2*s*c; 
plg(4,6) =(9*c*plg(4,5)-7.*plg(4,4))/2;
plg(4,7) =(11*c*plg(4,6)-8.*plg(4,5))/3;

if (~(((flags.sw(8)==0)&&(flags.sw(9)==0))&&(flags.sw(15)==0)))
    stloc = sin(hr*tloc);
    ctloc = cos(hr*tloc);
    s2tloc = sin(2*hr*tloc);
    c2tloc = cos(2*hr*tloc);
    s3tloc = sin(3*hr*tloc);
    c3tloc = cos(3*hr*tloc);
end

cd32 = cos(dr*(input.doy-p(32)));
cd18 = cos(2*dr*(input.doy-p(18)));
cd14 = cos(dr*(input.doy-p(14)));
cd39 = cos(2*dr*(input.doy-p(39)));

% F10.7 EFFECT
df = input.f107 - input.f107A;
dfa = input.f107A - 150;
t(1) =  p(20)*df*(1+p(60)*dfa) + p(21)*df*df + p(22)*dfa + p(30)*(dfa^2);
f1 = 1 + (p(48)*dfa +p(20)*df+p(21)*df*df)*flags.swc(2);
f2 = 1 + (p(50)*dfa+p(20)*df+p(21)*df*df)*flags.swc(2);

% TIME INDEPENDENT
t(2) = (p(2)*plg(1,3)+ p(3)*plg(1,5)+p(23)*plg(1,7)) + ...
       (p(15)*plg(1,3))*dfa*flags.swc(2) +p(27)*plg(1,2);

% SYMMETRICAL ANNUAL
t(3) = p(19)*cd32;

% SYMMETRICAL SEMIANNUAL
t(4) = (p(16)+p(17)*plg(1,3))*cd18;

% ASYMMETRICAL ANNUAL
t(5) = f1*(p(10)*plg(1,2)+p(11)*plg(1,4))*cd14;

% ASYMMETRICAL SEMIANNUAL
t(6) = p(38)*plg(1,2)*cd39;

% DIURNAL
if (flags.sw(8))
    t71 = (p(12)*plg(2,3))*cd14*flags.swc(6);
    t72 = (p(13)*plg(2,3))*cd14*flags.swc(6);
    t(7) = f2*((p(4)*plg(2,2) + p(5)*plg(2,4) + p(28)*plg(2,6) + t71) * ...
    ctloc + (p(7)*plg(2,2) + p(8)*plg(2,4) + p(29)*plg(2,6) ...
    + t72)*stloc);
end

% SEMIDIURNAL
if (flags.sw(9))
    t81 = (p(24)*plg(3,4)+p(36)*plg(3,6))*cd14*flags.swc(6);
    t82 = (p(34)*plg(3,4)+p(37)*plg(3,6))*cd14*flags.swc(6);
    t(8) = f2*((p(6)*plg(3,3)+ p(42)*plg(3,5) + t81)*c2tloc +(p(9)*plg(3,3) + p(43)*plg(3,5) + t82)*s2tloc);
end

% TERDIURNAL
if (flags.sw(15))
    t(14) = f2 * ((p(40)*plg(4,4)+(p(94)*plg(4,5)+p(47)*plg(4,7))*cd14*flags.swc(6))* s3tloc +(p(41)*plg(4,4)+(p(95)*plg(4,5)+p(49)*plg(4,7))*cd14*flags.swc(6))* c3tloc);
end

% magnetic activity based on daily ap
if (flags.sw(10)==-1)
    ap = input.ap_a;
    if (p(52)~=0)
        exp1 = exp(-10800*sqrt(p(52)*p(52))/(1+p(139)*(45-sqrt(input.g_lat*input.g_lat))));
        if (exp1>0.99999)
            exp1=0.99999;
        end
        if (p(25)<1E-4)
            p(25)=1E-4;
        end
        apt(1)=sg0(exp1,p,ap.a);
        % apt(2)=sg2(exp1,p,ap.a);
        % apt(3)=sg0(exp2,p,ap.a);
        % apt(4)=sg2(exp2,p,ap.a);
        
        if (flags.sw(10))
            t(9) = apt(1)*(p(51)+p(97)*plg(1,3)+p(55)*plg(1,5)+ ...
                (p(126)*plg(1,2)+p(127)*plg(1,4)+p(128)*plg(1,6))*cd14*flags.swc(6)+ ...
                (p(129)*plg(2,2)+p(130)*plg(2,4)+p(131)*plg(2,6))*flags.swc(8)* ...
                cos(hr*(tloc-p(132))));
        end
    end
else
    apd = input.ap-4;
    p44 = p(44);
    p45 = p(45);
    if (p44<0)
        p44 = 1E-5;
    end
    apdf = apd + (p45-1)*(apd + (exp(-p44 * apd) - 1)/p44);
    if (flags.sw(10))
        t(9)=apdf*(p(33)+p(46)*plg(1,3)+p(35)*plg(1,5)+ ...
            (p(101)*plg(1,2)+p(102)*plg(1,4)+p(103)*plg(1,6))*cd14*flags.swc(6)+ ...
            (p(122)*plg(2,2)+p(123)*plg(2,4)+p(124)*plg(2,6))*flags.swc(8)* ...
            cos(hr*(tloc-p(125))));
    end
end

if ((flags.sw(11))&&(input.g_long>-1000))
    % longitudinal
    if (flags.sw(12))
        t(11) = (1 + p(81)*dfa*flags.swc(2))* ...
            ((p(65)*plg(2,3)+p(66)*plg(2,5)+p(67)*plg(2,7)...
            +p(104)*plg(2,2)+p(105)*plg(2,4)+p(106)*plg(2,6)...
            +flags.swc(6)*(p(110)*plg(2,2)+p(111)*plg(2,4)+p(112)*plg(2,6))*cd14)* ...
            cos(dgtr*input.g_long) ...
            +(p(91)*plg(2,3)+p(92)*plg(2,5)+p(93)*plg(2,7)...
            +p(107)*plg(2,2)+p(108)*plg(2,4)+p(109)*plg(2,6)...
            +flags.swc(6)*(p(113)*plg(2,2)+p(114)*plg(2,4)+p(115)*plg(2,6))*cd14)* ...
            sin(dgtr*input.g_long));
    end
    
    % ut and mixed ut, longitude
    if (flags.sw(13))
        t(12)=(1+p(96)*plg(1,2))*(1+p(82)*dfa*flags.swc(2))*...
            (1+p(120)*plg(1,2)*flags.swc(6)*cd14)*...
            ((p(69)*plg(1,2)+p(70)*plg(1,4)+p(71)*plg(1,6))*...
            cos(sr*(input.sec-p(72))));
        t(12)= t(12) + flags.swc(12)*...
            (p(77)*plg(3,4)+p(78)*plg(3,6)+p(79)*plg(3,8))*...
            cos(sr*(input.sec-p(80))+2*dgtr*input.g_long)*(1+p(138)*dfa*flags.swc(2));
    end
    
    % ut, longitude magnetic activity
    if (flags.sw(14))
        if (flags.sw(10)==-1)
            if (p(52))
                t(13)=apt(1)*flags.swc(12)*(1.+p(133)*plg(1,2))*...
                    ((p(53)*plg(2,3)+p(99)*plg(2,5)+p(68)*plg(2,7))*...
                    cos(dgtr*(input.g_long-p(98))))...
                    +apt(1)*flags.swc(12)*flags.swc(6)*...
                    (p(134)*plg(2,2)+p(135)*plg(2,4)+p(136)*plg(2,6))*...
                    cd14*cos(dgtr*(input.g_long-p(137))) ...
                    +apt(1)*flags.swc(13)* ...
                    (p(56)*plg(1,2)+p(57)*plg(1,4)+p(58)*plg(1,6))*...
                    cos(sr*(input.sec-p(59)));
            end
        else
            t(13) = apdf*flags.swc(12)*(1+p(121)*plg(1,2))*...
                ((p(61)*plg(2,3)+p(62)*plg(2,5)+p(63)*plg(2,7))*...
                cos(dgtr*(input.g_long-p(64))))...
                +apdf*flags.swc(12)*flags.swc(6)* ...
                (p(116)*plg(2,2)+p(117)*plg(2,4)+p(118)*plg(2,6))* ...
                cd14*cos(dgtr*(input.g_long-p(119))) ...
                + apdf*flags.swc(13)* ...
                (p(84)*plg(1,2)+p(85)*plg(1,4)+p(86)*plg(1,6))* ...
                cos(sr*(input.sec-p(76)));
        end
    end
end

% parms not used: 82, 89, 99, 139-149
tinf = p(31);
for i=1:14
    tinf = tinf + abs(flags.sw(i+1))*t(i);
end

end


% ------------------------------------------------------------------- %
% ------------------------------- GLOB7S ---------------------------- %
% ------------------------------------------------------------------- %
function tt = glob7s(p,input,flags)
% VERSION OF GLOBE FOR LOWER ATMOSPHERE 10/26/99

global dfa plg ctloc stloc c2tloc s2tloc s3tloc c3tloc apdf apt

pset = 2;
t = zeros(14,1);
dr = 1.72142E-2;
dgtr = 1.74533E-2;
% confirm parameter set
if (p(100)==0)
    p(100)=pset;
end
if (p(100)~=pset)
    fprintf('Wrong parameter set for glob7s\n');
    tt = -1;
    return
end
for j=1:14
    t(j)=0;
end
cd32 = cos(dr*(input.doy-p(32)));
cd18 = cos(2*dr*(input.doy-p(18)));
cd14 = cos(dr*(input.doy-p(14)));
cd39 = cos(2*dr*(input.doy-p(39)));

% F10.7
t(1) = p(22)*dfa;

% time independent
t(2)=p(2)*plg(1,3) + p(3)*plg(1,5) + p(23)*plg(1,7) + p(27)*plg(1,2) + p(15)*plg(1,4) + p(60)*plg(1,6);

% SYMMETRICAL ANNUAL
t(3)=(p(19)+p(48)*plg(1,3)+p(30)*plg(1,5))*cd32;

% SYMMETRICAL SEMIANNUAL
t(4)=(p(16)+p(17)*plg(1,3)+p(31)*plg(1,5))*cd18;

% ASYMMETRICAL ANNUAL
t(5)=(p(10)*plg(1,2)+p(11)*plg(1,4)+p(21)*plg(1,6))*cd14;

% ASYMMETRICAL SEMIANNUAL
t(6)=(p(38)*plg(1,2))*cd39;

% DIURNAL
if (flags.sw(8))
    t71 = p(12)*plg(2,3)*cd14*flags.swc(6);
    t72 = p(13)*plg(2,3)*cd14*flags.swc(6);
    t(7) = ((p(4)*plg(2,2) + p(5)*plg(2,4) + t71) * ctloc + (p(7)*plg(2,2) + p(8)*plg(2,4) + t72) * stloc) ;
end

% SEMIDIURNAL
if (flags.sw(9))
    t81 = (p(24)*plg(3,4)+p(36)*plg(3,6))*cd14*flags.swc(6);
    t82 = (p(34)*plg(3,4)+p(37)*plg(3,6))*cd14*flags.swc(6);
    t(8) = ((p(6)*plg(3,3) + p(42)*plg(3,5) + t81) * c2tloc + (p(9)*plg(3,3) + p(43)*plg(3,5) + t82) * s2tloc);
end

% TERDIURNAL
if (flags.sw(15))
    t(14) = p(40) * plg(4,4) * s3tloc + p(41) * plg(4,4) * c3tloc;
end

% MAGNETIC ACTIVITY
if (flags.sw(10))
    if (flags.sw(10)==1)
        t(9) = apdf * (p(33) + p(46) * plg(1,3) * flags.swc(3));
    end
    if (flags.sw(10)==-1)
        t(9)=(p(51)*apt(1) + p(97)*plg(1,3) * apt(1)*flags.swc(3));
    end
end

% LONGITUDINAL
if (~((flags.sw(11)==0) || (flags.sw(12)==0) || (input.g_long<=-1000)))
    t(11) = (1 + plg(1,2)*(p(81)*flags.swc(6)*cos(dr*(input.doy-p(82)))...
    +p(86)*flags.swc(7)*cos(2*dr*(input.doy-p(87))))...
    +p(84)*flags.swc(4)*cos(dr*(input.doy-p(85)))...
    +p(88)*flags.swc(5)*cos(2*dr*(input.doy-p(89))))...
    *((p(65)*plg(2,3)+p(66)*plg(2,5)+p(67)*plg(2,7)...
    +p(75)*plg(2,2)+p(76)*plg(2,4)+p(77)*plg(2,6)...
    )*cos(dgtr*input.g_long)...
    +(p(91)*plg(2,3)+p(92)*plg(2,5)+p(93)*plg(2,7)...
    +p(78)*plg(2,2)+p(79)*plg(2,4)+p(80)*plg(2,6)...
    )*sin(dgtr*input.g_long));
end
tt=0;

for i=1:14
    tt = tt+abs(flags.sw(i+1))*t(i);   
end

end


% ------------------------------------------------------------------- %
% ------------------------------- GTD7 ------------------------------ %
% ------------------------------------------------------------------- %
function output = gtd7(input, flags)

global pt pd ps pdl ptl pma sam ptm pdm pavgm
global gsurf re dd dm04 dm16 dm28 dm32 dm40 dm01 dm14
global meso_tn1 meso_tn2 meso_tn3 meso_tgn1 meso_tgn2 meso_tgn3

tz = 0;

output = struct('d',zeros(9,1),'t',zeros(2));

mn3 = 5;
zn3 = [32.5,20,15,10,0];
mn2 = 4;
zn2 = [72.5,55,45,32.5];
zmix= 62.5;

flags = tselec(flags);

% Latitude variation of gravity (none for sw(3)=0)
xlat=input.g_lat;
if (flags.sw(3)==0)
    xlat=45;
end

[gsurf, re] = glatf(xlat);

xmm = pdm(3,5);

% THERMOSPHERE / MESOSPHERE (above zn2(1))
if (input.alt>zn2(1))
    altt = input.alt;
else
    altt = zn2(1);
end
tmp=input.alt;
input.alt = altt;
soutput = gts7(input, flags);
input.alt = tmp;
if (flags.sw(1)) % metric adjustment
    dm28m=dm28*1E6;
else
    dm28m=dm28;
end
output.t(1)=soutput.t(1);
output.t(2)=soutput.t(2);
if (input.alt>=zn2(1))
    for i=1:9
        output.d(i)=soutput.d(i);
    end
    return
end

% LOWER MESOSPHERE/UPPER STRATOSPHERE (between zn3(1) and zn2(1))
% Temperature at nodes and gradients at end nodes
% Inverse temperature a linear function of spherical harmonics
meso_tgn2(1)=meso_tgn1(2);
meso_tn2(1)=meso_tn1(5);
meso_tn2(2)=pma(1,1)*pavgm(1)/(1-flags.sw(21)*glob7s(pma(1,:), input, flags));
meso_tn2(3)=pma(2,1)*pavgm(2)/(1-flags.sw(21)*glob7s(pma(2,:), input, flags));
meso_tn2(4)=pma(3,1)*pavgm(3)/(1-flags.sw(21)*flags.sw(23)*glob7s(pma(3,:), input, flags));
meso_tgn2(2)=pavgm(9)*pma(10,1)*(1+flags.sw(21)*flags.sw(23)*glob7s(pma(10,:), input, flags))*meso_tn2(4)*meso_tn2(4)/((pma(3,1)*pavgm(3))^2);
meso_tn3(1)=meso_tn2(4);

if (input.alt<zn3(1))
% LOWER STRATOSPHERE AND TROPOSPHERE (below zn3(1))
% Temperature at nodes and gradients at end nodes
% Inverse temperature a linear function of spherical harmonics
	meso_tgn3(1)=meso_tgn2(2);
	meso_tn3(2)=pma(4,1)*pavgm(4)/(1-flags.sw(23)*glob7s(pma(4,:), input, flags));
	meso_tn3(3)=pma(5,1)*pavgm(5)/(1-flags.sw(23)*glob7s(pma(5,:), input, flags));
	meso_tn3(4)=pma(6,1)*pavgm(6)/(1-flags.sw(23)*glob7s(pma(6,:), input, flags));
	meso_tn3(5)=pma(7,1)*pavgm(7)/(1-flags.sw(23)*glob7s(pma(7,:), input, flags));
	meso_tgn3(2)=pma(8,1)*pavgm(8)*(1+flags.sw(23)*glob7s(pma(8,:), input, flags)) *meso_tn3(5)*meso_tn3(5)/((pma(7,1)*pavgm(7))^2);
end

% LINEAR TRANSITION TO FULL MIXING BELOW zn2(1)
dmc = 0;
if (input.alt>zmix)
    dmc = 1 - (zn2(1)-input.alt)/(zn2(1) - zmix);
end
dz28 = soutput.d(3);
	
%*** N2 density ***
dmr=soutput.d(3) / dm28m - 1;
[output.d(3), tz] = densm(input.alt, dm28m, xmm, tz, mn3, zn3, meso_tn3, meso_tgn3, mn2, zn2, meso_tn2, meso_tgn2);
output.d(3) = output.d(3) * (1 + dmr*dmc);

%*** HE density ***
dmr = soutput.d(1) / (dz28 * pdm(1,2)) - 1;
output.d(1) = output.d(3) * pdm(1,2) * (1 + dmr*dmc);

%*** O density ***
output.d(2) = 0;
output.d(9) = 0;

%*** O2 density ***
dmr = soutput.d(4) / (dz28 * pdm(4,2)) - 1;
output.d(4) = output.d(3) * pdm(4,2) * (1 + dmr*dmc);

%*** AR density ***
dmr = soutput.d(5) / (dz28 * pdm(5,2)) - 1;
output.d(5) = output.d(3) * pdm(5,2) * (1 + dmr*dmc);

%*** Hydrogen density ***
output.d(7) = 0;

%*** Atomic nitrogen density ***
output.d(8) = 0;

%*** Total mass density ***
output.d(6) = 1.66E-24 * (4 * output.d(1) + 16 * output.d(2) + 28 * output.d(3) + 32 * output.d(4) + 40 * output.d(5) + output.d(7) + 14 * output.d(8));

if (flags.sw(1))
    output.d(6)=output.d(6)/1000;
end

%*** temperature at altitude ***
[dd, tz] = densm(input.alt, 1, 0, tz, mn3, zn3, meso_tn3, meso_tgn3, mn2, zn2, meso_tn2, meso_tgn2);

output.t(2)=tz;

end


% ------------------------------------------------------------------- %
% ------------------------------- GTD7D ----------------------------- %
% ------------------------------------------------------------------- %
function output = gtd7d(input, flags)

output = gtd7(input, flags);
output.d(6) = 1.66E-24 * (4 * output.d(1) + 16 * output.d(2) + 28 * output.d(3) + 32 * output.d(4) + 40 * output.d(5) + output.d(7) + 14 * output.d(8) + 16 * output.d(9));
if (flags.sw(1))
    output.d(6)=output.d(6)/1000;
end

end


% ------------------------------------------------------------------- %
% -------------------------------- GHP7 ----------------------------- %
% ------------------------------------------------------------------- %
function output = ghp7(input, flags, press)

global gsurf re

bm = 1.3806E-19;
rgas = 831.4;
test = 0.00043;
ltest = 12;
pl = log10(press);
if (pl >= -5)
    if (pl>2.5)
        zi = 18.06 * (3 - pl);
    elseif ((pl>0.075) && (pl<=2.5))
        zi = 14.98 * (3.08 - pl);
    elseif ((pl>-1) && (pl<=0.075))
        zi = 17.80 * (2.72 - pl);
    elseif ((pl>-2) && (pl<=-1))
        zi = 14.28 * (3.64 - pl);
    elseif ((pl>-4) && (pl<=-2))
        zi = 12.72 * (4.32 -pl);
    else
        zi = 25.3 * (0.11 - pl);
    end
    cl = input.g_lat/90;
    cl2 = cl*cl;
    if (input.doy<182)
        cd = (1 - input.doy) / 91.25;
    else
        cd = (input.doy) / 91.25 - 3;
    end
    ca = 0;
    if ((pl > -1.11) && (pl<=-0.23))
        ca = 1;
    end
    if (pl > -0.23)
        ca = (2.79 - pl) / (2.79 + 0.23);
    end
    if ((pl <= -1.11) && (pl>-3))
        ca = (-2.93 - pl)/(-2.93 + 1.11);
    end
    z = zi - 4.87 * cl * cd * ca - 1.64 * cl2 * ca + 0.31 * ca * cl;
else
    z = 22 * ((pl + 4)^2) + 110;
end

% iteration loop
l = 0;
while(1)
    l = l+1;
    input.alt = z;
    output = gtd7(input, flags);
    z = input.alt;
    xn = output.d(1) + output.d(2) + output.d(3) + output.d(4) + output.d(5) + output.d(7) + output.d(8);
    p = bm * xn * output.t(2);
    if (flags.sw(1))
        p = p*1E-6;
    end
    diff = pl - log10(p);
    if (sqrt(diff*diff)<test)
        return
    end
    if (l==ltest)
        fprintf('ERROR: ghp7 not converging for press %e, diff %e',press,diff);
        return
    end
    xm = output.d(6) / xn / 1.66E-24;
    if (flags.sw(1))
        xm = xm * 1E3;
    end
    g = gsurf / ((1 + z/re)^2);
    sh = rgas * output.t(2) / (xm * g);
    
    % new altitude estimate using scale height %
    if (l <  6)
        z = z - sh * diff * 2.302;
    else
        z = z - sh * diff;
    end
end

end


% ------------------------------------------------------------------- %
% ------------------------------- GTS7 ------------------------------ %
% ------------------------------------------------------------------- %
function output = gts7(input, flags)
%     Thermospheric portion of NRLMSISE-00
%     See GTD7 for more extensive comments
%     alt > 72.5 km!

global pt pd ps pdl ptl pma sam ptm pdm pavgm
global gsurf re dd dm04 dm16 dm28 dm32 dm40 dm01 dm14
global meso_tn1 meso_tn2 meso_tn3 meso_tgn1 meso_tgn2 meso_tgn3

tz = 0;

output = struct('d',zeros(9,1),'t',zeros(2));

zn1 = [120, 110, 100, 90, 72.5];
mn1 = 5;
dgtr = 1.74533E-2;
dr = 1.72142E-2;
alpha = [-0.38, 0, 0, 0, 0.17, 0, -0.38, 0, 0];
altl = [200, 300, 160, 250, 240, 450, 320, 450];
za = pdl(2,16);
zn1(1) = za;
for j=1:9 
    output.d(j)=0;
end

% TINF VARIATIONS NOT IMPORTANT BELOW ZA OR ZN1(1)
if (input.alt>zn1(1))
    tinf = ptm(1)*pt(1) * ...
        (1+flags.sw(17)*globe7(pt,input,flags));
else
    tinf = ptm(1)*pt(1);
end
output.t(1)=tinf;

% GRADIENT VARIATIONS NOT IMPORTANT BELOW ZN1(5)
if (input.alt>zn1(5))
    g0 = ptm(4)*ps(1) * ...
        (1+flags.sw(20)*globe7(ps,input,flags));
else
    g0 = ptm(4)*ps(1);
end
tlb = ptm(2) * (1 + flags.sw(18)*globe7(pd(4,:),input,flags))*pd(4,1);
s = g0 / (tinf - tlb);

% Lower thermosphere temp variations not significant for density above 300 km
if (input.alt<300)
    meso_tn1(2)=ptm(7)*ptl(1,1)/(1-flags.sw(19)*glob7s(ptl(1,:), input, flags));
    meso_tn1(3)=ptm(3)*ptl(2,1)/(1-flags.sw(19)*glob7s(ptl(2,:), input, flags));
    meso_tn1(4)=ptm(8)*ptl(3,1)/(1-flags.sw(19)*glob7s(ptl(3,:), input, flags));
    meso_tn1(5)=ptm(5)*ptl(4,1)/(1-flags.sw(19)*flags.sw(21)*glob7s(ptl(4,:), input, flags));
    meso_tgn1(2)=ptm(9)*pma(9,1)*(1+flags.sw(19)*flags.sw(21)*glob7s(pma(9,:), input, flags))*meso_tn1(5)*meso_tn1(5)/((ptm(5)*ptl(4,1))^2);
else
    meso_tn1(2)=ptm(7)*ptl(1,1);
    meso_tn1(3)=ptm(3)*ptl(2,1);
    meso_tn1(4)=ptm(8)*ptl(3,1);
    meso_tn1(5)=ptm(5)*ptl(4,1);
    meso_tgn1(2)=ptm(9)*pma(9,1)*meso_tn1(5)*meso_tn1(5)/((ptm(5)*ptl(4,1))^2);
end

% N2 variation factor at Zlb
g28 = flags.sw(22)*globe7(pd(3,:), input, flags);

% VARIATION OF TURBOPAUSE HEIGHT
zhf = pdl(2,25)*(1+flags.sw(6)*pdl(1,25)*sin(dgtr*input.g_lat)*cos(dr*(input.doy-pt(14))));
output.t(1) = tinf;
xmm = pdm(3,5);
z = input.alt;

%*** N2 DENSITY ***

% Diffusive density at Zlb
db28 = pdm(3,1)*exp(g28)*pd(3,1);
% Diffusive density at Alt
[output.d(3), output.t(2)] = densu(z,db28,tinf,tlb,28,alpha(3),output.t(2),ptm(6),s,mn1,zn1,meso_tn1,meso_tgn1);
dd = output.d(3);
% Turbopause
zh28 = pdm(3,3)*zhf;
zhm28 = pdm(3,4)*pdl(2,6); 
xmd = 28-xmm;
% Mixed density at Zlb
[b28, tz] = densu(zh28,db28,tinf,tlb,xmd,(alpha(3)-1),tz,ptm(6),s,mn1, zn1,meso_tn1,meso_tgn1);
if ((flags.sw(16))&&(z<=altl(3)))
    %  Mixed density at Alt
    [dm28, ~] = densu(z,b28,tinf,tlb,xmm,alpha(3),tz,ptm(6),s,mn1,zn1,meso_tn1,meso_tgn1);
    %  Net density at Alt
    output.d(3) = dnet(output.d(3),dm28,zhm28,xmm,28);
end

%*** HE DENSITY ***

%   Density variation factor at Zlb
g4 = flags.sw(22)*globe7(pd(1,:), input, flags);
%  Diffusive density at Zlb
db04 = pdm(1,1)*exp(g4)*pd(1,1);
%  Diffusive density at Alt
[output.d(1), output.t(2)] = densu(z,db04,tinf,tlb, 4,alpha(1),output.t(2),ptm(6),s,mn1,zn1,meso_tn1,meso_tgn1);
dd = output.d(1);
if ((flags.sw(16)) && (z<altl(1)))
    %  Turbopause
    zh04 = pdm(1,3);
    %  Mixed density at Zlb
    [b04, output.t(2)] = densu(zh04,db04,tinf,tlb,4-xmm,alpha(1)-1.,output.t(2),ptm(6),s,mn1,zn1,meso_tn1,meso_tgn1);
    %  Mixed density at Alt
    [dm04, output.t(2)] = densu(z,b04,tinf,tlb,xmm,0,output.t(2),ptm(6),s,mn1,zn1,meso_tn1,meso_tgn1);
    zhm04 = zhm28;
    %  Net density at Alt
    output.d(1) = dnet(output.d(1),dm04,zhm04,xmm,4);
    %  Correction to specified mixing ratio at ground
    rl = log(b28*pdm(1,2)/b04);
    zc04 = pdm(1,5)*pdl(2,1);
    hc04 = pdm(1,6)*pdl(2,2);
    %  Net density corrected at Alt
    output.d(1) = output.d(1)*ccor(z,rl,hc04,zc04);
end

%*** O DENSITY ***

% Density variation factor at Zlb
g16 = flags.sw(22)*globe7(pd(2,:),input,flags);
%  Diffusive density at Zlb
db16 = pdm(2,1)*exp(g16)*pd(2,1);
% Diffusive density at Alt
[output.d(2), output.t(2)] = densu(z,db16,tinf,tlb, 16,alpha(2),output.t(2),ptm(6),s,mn1, zn1,meso_tn1,meso_tgn1);
dd = output.d(2);
if ((flags.sw(16)) && (z<=altl(2)))
    % Turbopause
    zh16 = pdm(2,3);
    % Mixed density at Zlb
    [b16, output.t(2)] = densu(zh16,db16,tinf,tlb,16-xmm,(alpha(2)-1),output.t(2),ptm(6),s,mn1,zn1,meso_tn1,meso_tgn1);
    %  Mixed density at Alt
    [dm16, output.t(2)] = densu(z,b16,tinf,tlb,xmm,0,output.t(2),ptm(6),s,mn1,zn1,meso_tn1,meso_tgn1);
    zhm16 = zhm28;
    % Net density at Alt
    output.d(2) = dnet(output.d(2),dm16,zhm16,xmm,16);
    rl = pdm(2,2)*pdl(2,17)*(1+flags.sw(2)*pdl(1,24)*(input.f107A-150));
    hc16 = pdm(2,6)*pdl(2,4);
    zc16 = pdm(2,5)*pdl(2,3);
    hc216 = pdm(2,6)*pdl(2,5);
    output.d(2) = output.d(2)*ccor2(z,rl,hc16,zc16,hc216);
    % Chemistry correction
    hcc16 = pdm(2,8)*pdl(2,14);
    zcc16 = pdm(2,7)*pdl(2,13);
    rc16 = pdm(2,4)*pdl(2,15);
    % Net density corrected at Alt
    output.d(2) = output.d(2)*ccor(z,rc16,hcc16,zcc16);
end

%*** O2 DENSITY ***

% Density variation factor at Zlb
g32= flags.sw(22)*globe7(pd(5,:), input, flags);
% Diffusive density at Zlb
db32 = pdm(4,1)*exp(g32)*pd(5,1);
% Diffusive density at Alt
[output.d(4), output.t(2)] = densu(z,db32,tinf,tlb, 32,alpha(4),output.t(2),ptm(6),s,mn1, zn1,meso_tn1,meso_tgn1);
dd = output.d(4);
if (flags.sw(16))
    if (z<=altl(4))
        % Turbopause
        zh32 = pdm(4,3);
        % Mixed density at Zlb
        [b32, output.t(2)] = densu(zh32,db32,tinf,tlb,32-xmm,alpha(4)-1.,output.t(2),ptm(6),s,mn1,zn1,meso_tn1,meso_tgn1);         
        % Mixed density at Alt
        [dm32, output.t(2)] = densu(z,b32,tinf,tlb,xmm,0,output.t(2),ptm(6),s,mn1,zn1,meso_tn1,meso_tgn1);
        zhm32 = zhm28;
        % Net density at Alt
        output.d(4) = dnet(output.d(4),dm32,zhm32,xmm,32);
        % Correction to specified mixing ratio at ground
        rl = log(b28*pdm(4,2)/b32);
        hc32 = pdm(4,6)*pdl(2,8);
        zc32 = pdm(4,5)*pdl(2,7);
        output.d(4) = output.d(4)*ccor(z,rl,hc32,zc32);
    end
    % Correction for general departure from diffusive equilibrium above Zlb
    hcc32 = pdm(4,8)*pdl(2,23);
    hcc232 = pdm(4,8)*pdl(1,23);
    zcc32 = pdm(4,7)*pdl(2,22);
    rc32 = pdm(4,4)*pdl(2,24)*(1+flags.sw(2)*pdl(1,24)*(input.f107A-150));
    % Net density corrected at Alt
    output.d(4) = output.d(4)*ccor2(z,rc32,hcc32,zcc32,hcc232);
end

%*** AR DENSITY ***

% Density variation factor at Zlb
g40 = flags.sw(22)*globe7(pd(6,:),input,flags);
% Diffusive density at Zlb
db40 = pdm(5,1)*exp(g40)*pd(6,1);
% Diffusive density at Alt
[output.d(5), output.t(2)] = densu(z,db40,tinf,tlb, 40,alpha(5),output.t(2),ptm(6),s,mn1,zn1,meso_tn1,meso_tgn1);
dd = output.d(5);
if ((flags.sw(16)) && (z<=altl(5)))
    % Turbopause
    zh40 = pdm(5,3);
    % Mixed density at Zlb
    [b40, output.t(2)] = densu(zh40,db40,tinf,tlb,40-xmm,alpha(5)-1.,output.t(2),ptm(6),s,mn1,zn1,meso_tn1,meso_tgn1);
    % Mixed density at Alt
    [dm40, output.t(2)] = densu(z,b40,tinf,tlb,xmm,0,output.t(2),ptm(6),s,mn1,zn1,meso_tn1,meso_tgn1);
    zhm40 = zhm28;
    % Net density at Alt
    output.d(5) = dnet(output.d(5),dm40,zhm40,xmm,40);
    % Correction to specified mixing ratio at ground
    rl = log(b28*pdm(5,2)/b40);
    hc40 = pdm(5,6)*pdl(2,10);
    zc40 = pdm(5,5)*pdl(2,9);
    % Net density corrected at Alt
    output.d(5) = output.d(5)*ccor(z,rl,hc40,zc40);
end

%*** HYDROGEN DENSITY ***

% Density variation factor at Zlb
g1 = flags.sw(22)*globe7(pd(7,:), input, flags);
% Diffusive density at Zlb
db01 = pdm(6,1)*exp(g1)*pd(7,1);
% Diffusive density at Alt
[output.d(7), output.t(2)] = densu(z,db01,tinf,tlb,1,alpha(7),output.t(2),ptm(6),s,mn1,zn1,meso_tn1,meso_tgn1);
dd = output.d(7);
if ((flags.sw(16)) && (z<=altl(7)))
    % Turbopause
    zh01 = pdm(6,3);
    % Mixed density at Zlb
    [b01, output.t(2)] = densu(zh01,db01,tinf,tlb,1-xmm,alpha(7)-1,output.t(2),ptm(6),s,mn1,zn1,meso_tn1,meso_tgn1);
    % Mixed density at Alt
    [dm01, output.t(2)] = densu(z,b01,tinf,tlb,xmm,0,output.t(2),ptm(6),s,mn1,zn1,meso_tn1,meso_tgn1);
    zhm01 = zhm28;
    % Net density at Alt
    output.d(7) = dnet(output.d(7),dm01,zhm01,xmm,1);
    % Correction to specified mixing ratio at ground
    rl = log(b28*pdm(6,2)*sqrt(pdl(2,18)*pdl(2,18))/b01);
    hc01 = pdm(6,6)*pdl(2,12);
    zc01 = pdm(6,5)*pdl(2,11);
    output.d(7) = output.d(7)*ccor(z,rl,hc01,zc01);
    % Chemistry correction
    hcc01 = pdm(6,8)*pdl(1,20);
    zcc01 = pdm(6,7)*pdl(2,19);
    rc01 = pdm(6,4)*pdl(2,21);
    % Net density corrected at Alt
    output.d(7) = output.d(7)*ccor(z,rc01,hcc01,zcc01);
end

%*** ATOMIC NITROGEN DENSITY ***

% Density variation factor at Zlb
g14 = flags.sw(22)*globe7(pd(8,:),input,flags);
% Diffusive density at Zlb
db14 = pdm(7,1)*exp(g14)*pd(8,1);
% Diffusive density at Alt
[output.d(8), output.t(2)] = densu(z,db14,tinf,tlb,14,alpha(8),output.t(2),ptm(6),s,mn1,zn1,meso_tn1,meso_tgn1);
dd = output.d(8);
if ((flags.sw(16)) && (z<=altl(8)))
    % Turbopause
    zh14 = pdm(7,3);
    % Mixed density at Zlb
    [b14, output.t(2)] = densu(zh14,db14,tinf,tlb,14-xmm,alpha(8)-1,output.t(2),ptm(6),s,mn1,zn1,meso_tn1,meso_tgn1);
    % Mixed density at Alt
    [dm14, output.t(2)] = densu(z,b14,tinf,tlb,xmm,0,output.t(2),ptm(6),s,mn1,zn1,meso_tn1,meso_tgn1);
    zhm14 = zhm28;
    % Net density at Alt
    output.d(8) = dnet(output.d(8),dm14,zhm14,xmm,14);
    % Correction to specified mixing ratio at ground
    rl = log(b28*pdm(7,2)*sqrt(pdl(1,3)*pdl(1,3))/b14);
    hc14 = pdm(7,6)*pdl(1,2);
    zc14 = pdm(7,5)*pdl(1,1);
    output.d(8) = output.d(8)*ccor(z,rl,hc14,zc14);
    % Chemistry correction
    hcc14 = pdm(7,8)*pdl(1,5);
    zcc14 = pdm(7,7)*pdl(1,4);
    rc14 = pdm(7,4)*pdl(1,6);
    % Net density corrected at Alt
    output.d(8) = output.d(8)*ccor(z,rc14,hcc14,zcc14);
end

%*** Anomalous OXYGEN DENSITY ***

g16h = flags.sw(22)*globe7(pd(9,:),input,flags);
db16h = pdm(8,1)*exp(g16h)*pd(9,1);
tho = pdm(8,10)*pdl(1,7);
[dd, output.t(2)] = densu(z,db16h,tho,tho,16,alpha(9),output.t(2),ptm(6),s,mn1, zn1,meso_tn1,meso_tgn1);
zsht = pdm(8,6);
zmho = pdm(8,5);
zsho = scalh(zmho,16,tho);
output.d(9) = dd*exp(-zsht/zsho*(exp(-(z-zmho)/zsht)-1));

% total mass density
output.d(6) = 1.66E-24*(4*output.d(1)+16*output.d(2)+28*output.d(3)+32*output.d(4)+40*output.d(5)+ output.d(7)+14*output.d(8));

% temperature
z = sqrt(input.alt*input.alt);
[~, output.t(2)] = densu(z,1, tinf, tlb, 0, 0,output.t(2),ptm(6), s, mn1, zn1, meso_tn1, meso_tgn1);

if (flags.sw(1))
    for i=1:9
        output.d(i)=output.d(i)*1E6;
    end
    output.d(6)=output.d(6)/1000;
end

end

function out = zeta(zz,zl)
global re
out = ((zz-zl)*(re+zl)/(re+zz));
end

% 3hr Magnetic activity functions
% Eq. A24d
function out = g0(a,p)
out = (a - 4 + (p(26) - 1) * (a - 4 + (exp(-sqrt(p(25)*p(25)) * (a - 4)) - 1) / sqrt(p(25)*p(25))));
end

% Eq. A24c
function out = sumex(ex)
out = (1 + (1 - (ex^19)) / (1 - ex) * (ex^0.5));
end

% Eq. A24a
function out = sg0(ex,p,ap)
out = (g0(ap(2),p)+(g0(ap(3),p)*ex+g0(ap(4),p)*ex*ex+...
       g0(ap(5),p)*(ex^3)+(g0(ap(6),p)*(ex^4)+...
       g0(ap(7),p)*(ex^12))*(1-(ex^8))/(1-ex)))/sumex(ex);
end

