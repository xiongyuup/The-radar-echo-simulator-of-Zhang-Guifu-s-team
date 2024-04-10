function [Moment_bias,Moment_std] = helperDataQuality(nexrad,moment,range,azimuth,flag)
% This function helperDataQuality is in support of WeatherRadarExample.
% helperDataQuality is to evaluate data quality based on estimated 'moment'
% and 'nexrad' truth field within an area defined by range and azimuth, for
% various radar moments specified by 'flag'. The returned Moment_bias and
% Moment_std are bias and standard deviation of estimated moments compared
% with NEXRAD truth data respectively.
    [rr,aa] = meshgrid(range/1e3,azimuth); 
    [xx,yy] = pol2cart(deg2rad(aa),rr);
    if strcmp(flag,'Z')
        % Error statistics for Z
        diff = moment.ZH-nexrad.ZH; 
        rows = reshape(diff,[],1);
        Moment_bias = mean(rows,'omitnan'); 
        Moment_std = std(rows,'omitnan'); 

        % Plot NEXRAD Z
        figure;
        subplot(121);
        sector = nexrad.ZH;
        sectorplot(xx,yy,sector,flag);
        title('NEXRAD Z_{H}');

        % Plot estimated Z
        subplot(122);
        calibrated = moment.ZH-Moment_bias;
        sectorplot(xx,yy,calibrated,flag);
        title('Simulated Z_{H}')

        % Plot range profile for Z
        figure;
        subplot(121);
        range_profile(range,nexrad.ZH,calibrated,flag)
        title('Z_{H} range profile');

        % Scatter plot for Z
        subplot(122);
        scatter_plot(nexrad.ZH,calibrated,flag)
        title('Scatter plot for Z_{H}');
         
    elseif strcmp(flag,'Vr')  
        % Error statistics for Vr
        diff = moment.Vr-nexrad.Vr; 
        rows = reshape(diff,[],1);
        Moment_bias = mean(rows,'omitnan'); 
        Moment_std = std(rows,'omitnan'); 

        % Plot NEXRAD Vr
        figure;
        subplot(121);
        sector = nexrad.Vr;
        sectorplot(xx,yy,sector,flag);
        title('NEXRAD V_{r}');

        % Plot estimated Vr
        subplot(122);
        calibrated = moment.Vr;
        sectorplot(xx,yy,calibrated,flag);
        title('Simulated V_{r}')

        % Plot range profile for Vr
        figure;
        subplot(121);
        range_profile(range,nexrad.Vr,calibrated,flag)
        title('Range profile for V_{r}');

        % Scatter plot for Vr
        subplot(122);
        scatter_plot(nexrad.Vr,calibrated,flag)
        title('Scatter plot for V_{r}');
        
    elseif strcmp(flag,'sigmav')
        % Error statistics for sigmav
        diff = moment.SW-nexrad.SW; 
        rows = reshape(diff,[],1);
        Moment_bias = mean(rows,'omitnan'); 
        Moment_std = std(rows,'omitnan'); 

        % Plot NEXRAD sigmav
        figure;
        subplot(121);
        sector = nexrad.SW;
        sectorplot(xx,yy,sector,flag);
        title('NEXRAD \sigma_{v}')

        % Plot estimated sigmav
        subplot(122);
        calibrated = moment.SW;
        sectorplot(xx,yy,calibrated,flag);
        title('Simulated \sigma_{v}')
         
        % Plot range profile for sigmav 
        figure;
        subplot(121);
        range_profile(range,nexrad.SW,calibrated,flag)
        title('Range profile for \sigma_{v}');

        % Scatter plot for sigmav
        subplot(122);
        scatter_plot(nexrad.SW,calibrated,flag)
        title('Scatter plot for \sigma_{v}');
        
    elseif strcmp(flag,'ZDR')
        % Error statistics for ZDR
        diff = moment.ZDR-nexrad.ZDR; 
        rows = reshape(diff,[],1);
        Moment_bias = mean(rows,'omitnan'); 
        Moment_std = std(rows,'omitnan'); 

        % Plot NEXRAD ZDR
        figure;
        subplot(121);
        sector = nexrad.ZDR;
        sectorplot(xx,yy,sector,flag);
        title('NEXRAD Z_{DR}')

        % Plot estimated ZDR
        subplot(122);
        calibrated = moment.ZDR-Moment_bias;
        sectorplot(xx,yy,calibrated,flag);
        title('Simulated Z_{DR}');

        % Plot range profile for ZDR
        figure;
        subplot(121);
        range_profile(range,nexrad.ZDR,calibrated,flag)
        title('Z_{DR} range profile');

        % Scatter plot for ZDR
        subplot(122);
        scatter_plot(nexrad.ZDR,calibrated,flag)
        title('Scatter plot for Z_{DR}');
        
        
    elseif strcmp(flag,'Rhohv')   
        % Error statistics for Rhohv
        diff = moment.Rhohv-nexrad.Rhohv; 
        rows = reshape(diff,[],1);
        Moment_bias = mean(rows,'omitnan'); 
        Moment_std = std(rows,'omitnan'); 

        % Plot NEXRAD Rhohv
        figure;
        subplot(121);
        sector = nexrad.Rhohv;
        sectorplot(xx,yy,sector,flag);
        title('NEXRAD \rho_{hv}');

        % Plot estimated Rhohv
        subplot(122);
        calibrated = moment.Rhohv-Moment_bias;
        sectorplot(xx,yy,calibrated,flag);
        title('Simulated \rho_{hv}')

        % Plot range profile for Rhohv
        figure;
        subplot(121);
        range_profile(range,nexrad.Rhohv,calibrated,flag)
        title('\rho_{hv} range profile');

        % Scatter plot for Rhohv
        subplot(122);
        scatter_plot(nexrad.Rhohv,calibrated,flag)
        title('Scatter plot for \rho_{hv}');
        
    elseif strcmp(flag,'Phidp')  
        % Error statistics for Phidp
        diff = moment.Phidp-nexrad.Phidp; 
        rows = reshape(diff,[],1);
        Moment_bias = mean(rows,'omitnan'); 
        Moment_std = std(rows,'omitnan'); 

        % Plot NEXRAD Phidp
        figure;
        subplot(121);
        sector = nexrad.Phidp;
        sectorplot(xx,yy,sector,flag);
        title('NEXRAD \phi_{DP}')

        % Plot estimated Phidp
        subplot(122);
        calibrated = moment.Phidp-Moment_bias;
        sectorplot(xx,yy,calibrated,flag);
        title('Simulated \phi_{DP}')

        % Plot range profile for Phidp 
        figure;
        subplot(121);
        range_profile(range,nexrad.Phidp,calibrated,flag)
        title('\phi_{DP} range profile');

        % Scatter plot for Phidp
        subplot(122);
        scatter_plot(nexrad.Phidp,calibrated,flag)
        title('Scatter plot for \phi_{DP}');

    end

    
function sectorplot(xx,yy,moment,flag)

    aznum = size(moment,1);
    if aznum == 360
        pcolor([xx;xx(1,:)],[yy;yy(1,:)],[moment(1:aznum,:);moment(1,:)]); 
    else
        pcolor(xx,yy,moment);
    end
    shading flat
    axis tight
    colormap('jet');
    colorbar;
    xlabel('Zonal Distance (km)');
    ylabel('Meridional Distance (km)');

    if strcmp(flag,'Z') 
        cmax = 60;
        cmin = -10;
        c = colorbar;
        c.Label.String = 'dBZ';
    elseif strcmp(flag,'Vr') 
        cmax = 15;
        cmin = 0;   
        c = colorbar;
        c.Label.String = 'm/s';
    elseif strcmp(flag,'sigmav') 
        cmax = 3;
        cmin = 0; 
        c = colorbar;
        c.Label.String = 'm/s';
    elseif strcmp(flag,'ZDR') 
        cmax = 5;
        cmin = -3;
        c = colorbar;
        c.Label.String = 'dB';
    elseif strcmp(flag,'Rhohv') 
        cmax = 1.0;
        cmin = 0.9;
    elseif strcmp(flag,'Phidp') 
        cmax = 90;
        cmin = 0;
        c = colorbar;
        c.Label.String = 'degree';
    end
    caxis([cmin,cmax]);

    
function range_profile(range,nexrad,moment,flag)
    plot(range/1e3, nexrad(1,:),'bx-','LineWidth',2); 
    hold on;
    plot(range/1e3, moment(1,:),'ro-','LineWidth',2); 
    legend('NEXRAD','Estimate','Location','southeast');
    xlow = floor(range(1)/1e3);
    xhigh = ceil(range(end)/1e3);
    xlim([xlow xhigh]);
    xlabel('Range (km)');
    if strcmp(flag,'Z') 
        ylabel('Z_{H} (dBZ)');
    elseif strcmp(flag,'Vr') 
        ylabel('V_{r} (m/s)'); 
    elseif strcmp(flag,'sigmav') 
        ylabel('sigma_{v} (m/s)');
    elseif strcmp(flag,'ZDR') 
        ylabel('Z_{DR} (dB)');
    elseif strcmp(flag,'Rhohv') 
        ylabel('\rho_{hv}');
    elseif strcmp(flag,'Phidp') 
        ylabel('\phi_{DP} (deg)');
    end
    
function scatter_plot(nexrad,moment,flag)                      
    moment1 = reshape(moment,[],1);
    nexrad1 = reshape(nexrad,[],1);
    scatter(nexrad1,moment1);
    zmin = min(min(nexrad1),min(moment1),'omitnan');
    zmax = max(max(nexrad1),max(moment1),'omitnan');
    if strcmp(flag,'Rhohv')
        x = zmin:0.001:zmax;
    else
        x = floor(zmin):ceil(zmax);
    end
    y = x;
    hold on;
    plot(x,y);
    grid on;

    if strcmp(flag,'Z') 
        xlabel('NEXRAD Z_{H} (dBZ)');
        ylabel('Simulated Z_{H} (dBZ)');
    elseif strcmp(flag,'Vr') 
        xlabel('NEXRAD V_{r} (m/s)');
        ylabel('Simulated V_{r} (m/s)');
    elseif strcmp(flag,'sigmav') 
        xlabel('NEXRAD \sigma_{v} (m/s)');
        ylabel('Simulated \sigma_{v} (m/s)');
    elseif strcmp(flag,'ZDR') 
        xlabel('NEXRAD Z_{DR} (dB)');
        ylabel('Simulated Z_{DR} (dB)');
    elseif strcmp(flag,'Rhohv') 
        xlabel('NEXRAD \rho_{hv}');
        ylabel('Simulated \rho_{hv}');
    elseif strcmp(flag,'Phidp') 
        xlabel('NEXRAD \phi_{DP} (deg)');
        ylabel('Simulated \phi_{DP} (deg)');
    end
    grid on;