c = 3e8;            %speed of light
range_max = 180;    %max detection range
range_res = 1;
tm = 5.5*(2*range_max/c); %sweep time

bw = c/(2*range_res);

sweep_slope = bw/tm;
v_max = 180*1000/3600;  %target max velocity
fc = 77e9;          %radar frequency
lambda = c/fc;      %radar wavelength
fs = 72e6;          %sampling rate

fbmx = (2*range_max*bw/tm) + ((2*v_max/c)*fc);

ts = 1/fs;
t = 0:ts:tm;

chirps = 64;        

samples = ceil(tm*fs);  %samples in one chirp

%% target

R0 = 110; %range in meters
V = 36; %radial velocity, m/s
R2 = 50; %range in meters
V2 = 45; %radial velocity, m/s