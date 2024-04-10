function [nexrad_aoi,Ns] = helperSelectAOI(nexrad,az1,az2,rg1,rg2,blind_rg,range_res,num_az,num_bin)
% This function helperSelectAOI is in support of WeatherRadarExample.
% helperSelectAOI is to select an area of interest(AOI) that targets
% occupy. nexrad is a struct object holding NEXRAD moment data. az1 and az2
% are starting and ending azimuth angles specified in degree, rg1 and rg2
% are starting and ending ranges specified in meter. blind_rg is NEXRAD
% blind range, range_res is radar range resolution, num_az is the number of
% azimuth angles in the AOI, num_bin is the number of range bins in the
% AOI. The returned nexrad_aoi is a struct object holding AOI data, Ns is
% the number of available weather targets in space by removing NaN data.

    % Azimuth angle for each radial data in the NEXRAD coordinate.
    az_nexrad = mod(90-nexrad.azimuth,360);  
    % Transform from Cartesian coordinate to NEXRAD coordinate, to get the
    % starting and ending radial indexes and range bin indexes of the AOI.
    [r1,r2,b1,b2,aznum,rgnum] = cart2rad(az_nexrad,az1,az2,rg1,rg2,blind_rg,range_res);
    % Define an area of interest(AOI) in terms of azimuth and range.
    nexrad_aoi.ZH = nexrad.ZH(r1:r2,b1:b2);
    nexrad_aoi.Vr = nexrad.Vr(r1:r2,b1:b2);
    nexrad_aoi.SW = nexrad.SW(r1:r2,b1:b2);
    nexrad_aoi.ZDR = nexrad.ZDR(r1:r2,b1:b2);
    nexrad_aoi.Rhohv = nexrad.Rhohv(r1:r2,b1:b2);
    nexrad_aoi.Phidp = nexrad.Phidp(r1:r2,b1:b2);
    
    nexrad_aoi.r1 = r1;
    nexrad_aoi.r2 = r2;
    nexrad_aoi.b1 = b1;
    nexrad_aoi.b2 = b2;
    nexrad_aoi.aznum = aznum;
    nexrad_aoi.rgnum = rgnum;
    
    % Calculate the number of available weather targets in the 2D space by
    % removing nan data.
    r0 = (num_az-aznum)/2;
    nexrad_aoi.rlow = round(r1-r0);
    nexrad_aoi.rup = round(r2+r0);
    Ns = num_az*num_bin-sum(sum(isnan(nexrad.ZH(r1-r0:r2+r0,:)))); 
    
    
function [r1,r2,b1,b2,aznum,rgnum] = cart2rad(Az0,az1,az2,rg1,rg2,blind_rg,range_res)
    az1_rad = mod(90-az1,360);
    az2_rad = mod(90-az2,360);
    r1 = min(find(abs(Az0-az1_rad)<0.5),find(abs(Az0-az2_rad)<0.5));
    r2 = max(find(abs(Az0-az1_rad)<0.5),find(abs(Az0-az2_rad)<0.5));
    b1 = (rg1-blind_rg)/range_res+1;
    b2 = (rg2-blind_rg)/range_res+1;
    aznum = r2-r1+1;  % Number of azimuth angles in AOI
    rgnum = b2-b1+1;  % Number of range bins in AOI
