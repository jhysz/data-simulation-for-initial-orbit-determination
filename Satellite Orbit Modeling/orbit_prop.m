function Eph_eci = orbit_prop(obs_k)

format long g
close all
profile on
tic

global const Cnm Snm AuxParam swdata eopdata PC

SAT_Const
load DE430Coeff.mat
PC = DE430Coeff;

% read Earth gravity field coefficients
Cnm = zeros(181,181);
Snm = zeros(181,181);
fid = fopen('GGM03S.txt','r');
for n=0:180
    for m=0:n
        temp = fscanf(fid,'%d %d %f %f %f %f',[6 1]);
        Cnm(n+1,m+1) = temp(3);
        Snm(n+1,m+1) = temp(4);
    end
end
fclose(fid);

% read Earth orientation parameters
fid = fopen('eop19620101.txt','r');
%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s 
%  ----------------------------------------------------------------------------------------------------
eopdata = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
fclose(fid);

% read space weather data
fid = fopen('sw19571001.txt','r');
%  ---------------------------------------------------------------------------------------------------------------------------------
% |                                                                                             Adj     Adj   Adj   Obs   Obs   Obs 
% | yy mm dd BSRN ND Kp Kp Kp Kp Kp Kp Kp Kp Sum Ap  Ap  Ap  Ap  Ap  Ap  Ap  Ap  Avg Cp C9 ISN F10.7 Q Ctr81 Lst81 F10.7 Ctr81 Lst81
%  ---------------------------------------------------------------------------------------------------------------------------------
swdata = fscanf(fid,'%4i %3d %3d %5i %3i %3i %3i %3i %3i %3i %3i %3i %3i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4f %2i %4i %6f %2i %6f %6f %6f %6f %6f',[33 inf]);
fclose(fid);

% model parameters
AuxParam = struct('Mjd_UTC',0,'area_solar',0,'area_drag',0,'mass',0,'Cr',0,...
                  'Cd',0,'n',0,'m',0,'sun',0,'moon',0,'sRad',0,'drag',0,...
                  'planets',0,'SolidEarthTides',0,'OceanTides',0,'Relativity',0);
              
% epoch state (Envisat)
fid = fopen('InitialState.txt','r');

tline = fgetl(fid);
year = str2num(tline(1:4));
mon = str2num(tline(6:7));
day = str2num(tline(9:10));
hour = str2num(tline(12:13));
min = str2num(tline(15:16));
sec = str2num(tline(18:23));
Y0 = zeros(1,6);
for i = 1:6
    tline = fgetl(fid);
    Y0(i) = obs_k(i);
end
tline = fgetl(fid);
AuxParam.area_solar = str2num(tline(49:end));
tline = fgetl(fid);
AuxParam.area_drag = str2num(tline(38:end));
tline = fgetl(fid);
AuxParam.mass = str2num(tline(19:end));
tline = fgetl(fid);
AuxParam.Cr = str2num(tline(5:end));
tline = fgetl(fid);
AuxParam.Cd = str2num(tline(5:end));
fclose(fid);

% epoch
Mjd_UTC = Mjday(year, mon, day, hour, min, sec);
% Y0 = ECEF2ECI(Mjd_UTC, Y0); % ECI 

AuxParam.Mjd_UTC = Mjd_UTC;
AuxParam.n       = 36;
AuxParam.m       = 36;
AuxParam.sun     = 1;
AuxParam.moon    = 1;
AuxParam.planets = 1;
AuxParam.sRad    = 1;
AuxParam.drag    = 1;
AuxParam.SolidEarthTides = 0;
AuxParam.OceanTides = 0;
AuxParam.Relativity = 0;

Mjd0   = Mjd_UTC;
Step   = 2;   % [s]
N_Step = 240; 

% shorten PC, eopdata, swdata, Cnm, and Snm to decrease CPU time
num = fix(N_Step*Step/86400)+2;
JD = Mjd_UTC+2400000.5;
i = find(PC(:,1)<=JD & JD<=PC(:,2),1,'first');
PC = PC(i:i+num,:);
mjd = (floor(Mjd_UTC));
i = find(mjd==eopdata(4,:),1,'first');
eopdata = eopdata(:,i:i+num);
i = find((year==swdata(1,:)) & (mon==swdata(2,:)) & (day==swdata(3,:)),1,'first');
swdata = swdata(:,i-3:i+num);
Cnm = Cnm(1:AuxParam.n+1,1:AuxParam.n+1);
Snm = Snm(1:AuxParam.n+1,1:AuxParam.n+1);

% propagation
Eph = Ephemeris(Y0, N_Step, Step);


Eph_eci = zeros(N_Step+1,6);
fid = fopen('SatelliteStates.txt','w');
for i=1:N_Step+1
    [year,mon,day,hr,min,sec] = invjday(Mjd0+Eph(i,1)/86400+2400000.5);
    fprintf(fid,'  %4d/%2.2d/%2.2d  %2d:%2d:%6.3f',year,mon,day,hr,min,sec);
    fprintf(fid,'  %14.3f%14.3f%14.3f%12.3f%12.3f%12.3f\n',...
            Eph(i,2),Eph(i,3),Eph(i,4),Eph(i,5),Eph(i,6),Eph(i,7));
    Eph_eci(i,:) =  Eph(i,2:7);
end
fclose(fid);

toc

