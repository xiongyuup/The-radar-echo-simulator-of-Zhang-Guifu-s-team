function meanRCS = helperdBZ2RCS(theta,range,lambda,tau,dBZ,c)
% This function helperdBZ2RCS is in support of WeatherRadarExample.
% helperdBZ2RCS is to transform reflectivity in the resolution volume to
% mean RCS. theta is the radar 3dB beamwidth, range is the distance to
% target, lambda is the radar wavelength, tau is the radar pulse width, dBZ
% is the NEXRAD reflectivity, c is the speed of light. The returned meanRCS
% is the mean RCS in the resolution volume.
Z = 10.^(dBZ/10)*1e-18; 
% reflectivity
nta = pi^5*0.93*Z/(lambda^4); 
% resolution volume size
Volume = range.^2*c*tau*pi*deg2rad(theta)^2/(16*log(2));
% mean RCS
meanRCS = nta*Volume;
 