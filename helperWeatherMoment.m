function moment1 = helperWeatherMoment(rxh,rxv,nexrad,N,L,prt,lambda) 
% This function helperWeatherMoment is in support of WeatherRadarExample.
% helperWeatherMoment is to calculate radar moments based on pulse pair
% processing method. rxh and rxv are 2D complex data cubes, where the first
% dimension is range, and the second dimension is azimuth. nexrad is a
% struct object holding moment truth data, N is the number of processing
% pulses, L is the number of realizations. prt is the pulse repetition
% time, lambda is the radar wavelength. The returned moment1 is a 2D real
% data cube, containing estimated radar moments.
    s = size(nexrad.ZH);
    moment1 = struct('ZH',zeros(s),'Vr',zeros(s),'SW',zeros(s),'ZDR',zeros(s),'Rhohv',zeros(s),'Phidp',zeros(s));
    for kk = 1:nexrad.aznum
        % matrix size: rxh:M*1; Vh and Vhn:L*N; P and R:N*N; iqh:N*L; moment:M*1
        % N: pulse number; L: realization number
        M = size(rxh(:,kk),1);
        M;
        sigmah = sqrt((abs(rxh(:,kk))).^2/2);
        sigmav = sqrt((abs(rxv(:,kk))).^2/2);
        delta = 1e-4;
        % Calculate the maximum unambiguous velocity
        va = lambda/4/prt;
        for m = 1:M
            if isnan(nexrad.SW(kk,m))
                continue
            end
            Vh = sigmah(m)*(randn(L,N)+1i*randn(L,N));
            Vv = sigmav(m)*(randn(L,N)+1i*randn(L,N));  
            % Calculate covariance matrix   
            % P2P correlation:Dr.Doviak book eq.(6.4)
            pvec = exp(-8*(pi*nexrad.SW(kk,m)*(0:N-1)*prt/lambda).^2);  
            phivec = exp(1i*-4*pi*nexrad.Vr(kk,m)*(0:N-1)*prt/lambda);
            P = toeplitz(pvec.*phivec);
            P = P+delta*eye(size(P));
            % disp(P);
            % [R,p] = chol(P);
            % R = cholesky(P);
            [R, D] = eig(P);
            % [LL, D] = ldl(P);
            % if p~=0
            %     disp('feizhengding')
            % end
            % size(P)
            % R=LL;
            Vh = Vh*R;
            Vv = Vv*R;

            % Generate dualpol time series
            Vhn = Vh;
            % size(Vh)
            Vvn = (nexrad.Rhohv(kk,m)*Vh+sqrt(1-(nexrad.Rhohv(kk,m))^2)*Vv)*exp(1j*nexrad.Phidp(kk,m)*pi/180)/sqrt(10.^(nexrad.ZDR(kk,m)/10));
            iqh = Vhn.'; 
            iqv = Vvn.';
            size(iqh)
            % mom:(N-1)*1, N-1 is number of pulse pair
            mom = PulsePairProcessing(iqh,iqv,va);
            % average on pulse pair
            moment.ZH(m) = mean(mom.ZH,'omitnan');    
            moment.Vr(m) = mean(mom.Vr,'omitnan');
            moment.SW(m) = mean(mom.SW,'omitnan');
            moment.ZDR(m) = mean(mom.ZDR,'omitnan');
            moment.Rhohv(m) = mean(mom.Rhohv,'omitnan');
            moment.Phidp(m) = mean(mom.Phidp,'omitnan');
        end

        % Range correction for reflectivity
        range = (nexrad.b1-1:nexrad.b2-1)*250 + 2000;
        
        moment.ZH = moment.ZH + 20*log10(range/1e3);

        % For each NaN data in NEXRAD, radar estimate is also assigned as nan.
        [row,col] = find(isnan(nexrad.ZH(kk,:)));
        if not(isempty(row))
            rr = size(row,2);
            for nn = 1:rr
                moment.ZH(col(nn)) = nan;
                moment.Vr(col(nn)) = nan;
                moment.SW(col(nn)) = nan;
                moment.ZDR(col(nn)) = nan;
                moment.Rhohv(col(nn)) = nan;
                moment.Phidp(col(nn)) = nan;
            end
        end
        moment1.ZH(kk,:) = moment.ZH;
        moment1.Vr(kk,:) = moment.Vr;
        moment1.SW(kk,:) = moment.SW;
        moment1.ZDR(kk,:) = moment.ZDR;
        moment1.Rhohv(kk,:) = moment.Rhohv;
        moment1.Phidp(kk,:) = moment.Phidp;
    end

function moment = PulsePairProcessing(iqh,iqv,va)
    % matrix size: iqh:N*L; Rhh:N*N; R0_hh:(N-1)*1; moment:(N-1)*1
    np = size(iqh,1);
    ns = size(iqh,2);  % number of realization

    Rhh = iqh*iqh'/ns; % average on realization
    % size(Rhh)
    Rvv = iqv*iqv'/ns;
    Rhv = iqv*iqh'/ns;

    R0_hh = diag(Rhh(1:end-1,1:end-1));
    R0_vv = diag(Rvv(1:end-1,1:end-1));
    C0_hv = diag(Rhv(1:end-1,1:end-1));
    % Get the elements on the first subdiagonal (k=-1) of A. 
    % The result has one fewer element than the main diagonal.
    R1_hh = diag(Rhh,-1);
    R1_vv = diag(Rvv,-1);
    % size(R1_vv)
    % Pulse pair processing
    moment.Vr = -va./pi*angle(R1_hh); 
    % size(moment.Vr)
    moment.Vrv = -va./pi*angle(R1_vv);
    moment.SW = (2*va/pi/sqrt(2)).*(abs(log(R0_hh./abs(R1_hh)))).^0.5; 
    moment.SWv = (2*va/pi/sqrt(2)).*(abs(log(R0_vv./abs(R1_vv)))).^0.5; 
    moment.ZH = 10*log10(R0_hh);
    moment.ZV = 10*log10(R0_vv);
    moment.ZDR = moment.ZH-moment.ZV;
    moment.Rhohv = abs(C0_hv)./(R0_hh).^0.5./(R0_vv).^0.5;
    moment.Phidp = angle(C0_hv)*180/pi; 