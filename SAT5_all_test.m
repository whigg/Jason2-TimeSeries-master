% function [ TS, mat,TS_final,Cor] = Altprocess( vlat,vlon,SR,latmin, latmax, lonmin, lonmax, SAT,track,track2)

% ===========================================================
%
%  This function process satellite altimetry data for a given coordinate to
%  obtain water level time series
%

%  input:
%       vlat: Virtual station's latitude [deg]
%       vlon: virtual station's longitude [deg]
%       SR: search radius for virtual station [m]
%
%       ofndr: applying off nadir correction 1, not applying 0
%       latmin: minimum latitude [deg]
%       latmax: maximum latitude  [deg]
%       lonmin: minimum longitude [deg]
%       lonmax: maximum longitude [deg]
%       SAT: choice of satellite, 0: TOPEX 1: ENVISAT, 2: TOPEX new 3:
%       Jason 1, 4:SARAL, 5: Jason 2
%       track: track number []
%       track2: track number in ENVISAT ext []

%
%  output:
%       TS_final:
%               Column 1: year
%               Column 2: month
%               Column 3: surface water height from combined retracking
%               approach
%               Column 4: error of estimated water level
%       mat:

%       TS (ENVISAT):
%               Column 1: year
%               Column 2: month
%               Column 3: Time of measurement
%               Column 4: ocean retracker 1HZ
%               Column 5: median ocean retracker 20 HZ
%               Column 6: mean ocean retracker 20 HZ
%               Column 7: std (error) ocean retracker 20 HZ
%               Column 8: median Ice1 retracker 20 HZ
%               Column 9: mean Ice1 retracker 20 HZ
%               Column 10: std (error) Ice1 retracker 20 HZ
%               Column 11: median Ice2 retracker 20 HZ
%               Column 12: mean Ice2 retracker 20 HZ
%               Column 13: std (error) Ice2 retracker 20 HZ
%               Column 14: median SeaIce retracker 20 HZ
%               Column 15: mean SeaIce retracker 20 HZ
%               Column 16: std (error) SeaIce retracker 20 HZ
%        Time:
%               Time of measurements
%        mat: processed mat file within the defined box



% Mohammad J. Tourian
% tourian@gis.uni-stuttgart.de
% January 2014
% ===========================================================
load EGM2008_5 % loading EGM 2008 to provide height from Geoid

% if nargin>5

    %  track2=fix(track2)/2;


%     if SAT==5 % JASON 2
if lonmin<0
    lonmin=360+lonmin;
end
if lonmax<0
    lonmax=360+lonmax;
end
mat=JA2_PH_crt(latmin,latmax, lonmin, lonmax,track);
fclose all;
%     end

%% finding 20 Hz measurements inside the virtual station

[N1,E1]=ell2utm(vlat*pi/180,vlon*pi/180);
[m n]=size(mat);




if SAT==5
    if isempty(mat)==0
        pp=1;ic=0;

        TmSri=[];
        for i=1:m
            [s1 s2]=size(mat{i,2});
            if s1>1

                for j=1:s1
                    [N2,E2]=ell2utm(mat{i,2}(j,3)*pi/(180),mat{i,2}(j,2)*pi/(180));
                    Dist=sqrt((N1-N2)^2+(E1-E2)^2);
                    if Dist<=SR
                        if length(mat{i,2}(j,:))==179
                            TmSri(pp,:)=mat{i,2}(j,:);
                        end

                        ic=1;
                    end
                    if ic==1
                        pp=pp+1;
                    end
                    ic=0;
                end
            end
        end



        %% correcting the time series
        if isempty(TmSri)==0
            [m n]=size(TmSri);
            g=2;
            TmSric(1,:)=TmSri(1,:);
            for i=2:m
                f=find(TmSric(:,1)==TmSri(i,1) & TmSric(:,2)==TmSri(i,2) & TmSric(:,3)==TmSri(i,3));
                if isempty(f)==1
                    TmSric(g,:)=TmSri(i,:);
                    g=g+1;
                end
            end
        else
            TmSric=[];
        end
    else
        TmSric=[];
    end

end





%%Geoid height computation

if SAT==5
    if isempty(TmSric)==0
        lat=TmSric(:,3);
        lon=TmSric(:,2);
        f=find( lon>180);
        lon(f)=lon(f)-360;
        TmSric(:,180)= griddata(XX(:),YY(:),GH(:),lon,lat);
    end

end



%% geophysical correction



if SAT==5
    if isempty(TmSric)==0
        [m,n]=size(TmSric);

        for i=1:m
            %% Altitude
            % 1Hz
            Alt_1hz(i,1)=TmSric(i,26);
            % 20 Hz
            % Alt_20hz(i,1:20)=TmSric(i,154:173)/10^4;

            %% Ocean Range

            %%%%%%%%%%%%%%% Ku %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % 1Hz
            %% Corrections
            %Inverse barometric correction
            if TmSric(i,99)~=32767 & TmSric(i,99)~=-32768 & TmSric(i,99)~=2.147483647000000e+05
                InvBar_ku=TmSric(i,99);
            else
                InvBar_ku=0;
            end
            %Sea state bias
            if TmSric(i,51)~=32767 & TmSric(i,51)~=-32768 & TmSric(i,51)~=2.147483647000000e+05
                SeSbias_ku=TmSric(i,51);
            else
                SeSbias_ku=0;
            end
            % Ionospheric correction
            if TmSric(i,47)~=32767 & TmSric(i,47)~=-32768 & TmSric(i,47)~=2.147483647000000e+05
                IonCor_ku=TmSric(i,47);
            else
                IonCor_ku=0;
            end
            % Ocean Tide
            if TmSric(i,101)~=32767 & TmSric(i,101)~=-32768 & TmSric(i,101)~=2.147483647000000e+05
                OcTide_ku=TmSric(i,101);
            else
                OcTide_ku=0;
            end
            % Polar Tide
            if TmSric(i,110)~=32767 & TmSric(i,110)~=-32768 & TmSric(i,110)~=2.147483647000000e+05
                PoTide_ku=TmSric(i,110);
            else
                PoTide_ku=0;
            end
            % Earth Tide
            if TmSric(i,109)~=32767 & TmSric(i,109)~=-32768 & TmSric(i,109)~=2.147483647000000e+05
                ETide_ku=TmSric(i,109);
            else
                ETide_ku=0;
            end
            % Wet tropospheric correction
            if TmSric(i,42)~=32767 & TmSric(i,42)~=-32768 & TmSric(i,42)~=2.147483647000000e+05
                WTCor_ku=TmSric(i,42);
            else
                WTCor_ku=0;
            end
            % Dry tropospheric correction
            if TmSric(i,40)~=32767 & TmSric(i,40)~=-32768 & TmSric(i,40)~=2.147483647000000e+05
                DTCor_ku=TmSric(i,40);
            else
                DTCor_ku=0;
            end
            %%%% performing corrections
            Correction_ku=InvBar_ku+SeSbias_ku+IonCor_ku+OcTide_ku+PoTide_ku+ETide_ku+WTCor_ku+DTCor_ku+TmSric(i,180);
            Cor(i,:)=[InvBar_ku SeSbias_ku IonCor_ku OcTide_ku PoTide_ku ETide_ku WTCor_ku DTCor_ku];

            Range_ku(i,1)=TmSric(i,28)+Correction_ku;
            Range_c(i,1)=TmSric(i,29)+Correction_ku;
            Range_oce3_ku(i,1)=TmSric(i,30)+Correction_ku;
            Range_oce3_c(i,1)=TmSric(i,31)+Correction_ku;
            Range_red3_ku(i,1)=TmSric(i,32)+Correction_ku;
            Range_red3_c(i,1)=TmSric(i,33)+Correction_ku;
            Range_ice3_ku(i,1)=TmSric(i,34)+Correction_ku;
            Range_ice3_c(i,1)=TmSric(i,35)+Correction_ku;





        end

        %% SSH


        SSH_ku=Alt_1hz-Range_ku;
        SSH_c=Alt_1hz-Range_c;
        SSH_oce3_ku=Alt_1hz-Range_oce3_ku;
        SSH_oce3_c=Alt_1hz-Range_oce3_c;
        SSH_red3_ku=Alt_1hz-Range_red3_ku;
        SSH_red3_c=Alt_1hz-Range_red3_c;
        SSH_ice3_ku=Alt_1hz-Range_ice3_ku;
        SSH_ice3_c=Alt_1hz-Range_ice3_c;



        % Creating Time matrix
        s=1;
        ind=1;

        a=datestr(TmSric(1,1)/86400+datenum(2000,1,1),'yyyymmdd');
        y=str2num(a(1:4));
        mo=str2num(a(5:6));
        d=str2num(a(7:8));
        for i=2:m
            a=datestr(TmSric(i,1)/86400+datenum(2000,1,1),'yyyymmdd');
            y(i,1)=str2num(a(1:4));
            mo(i,1)=str2num(a(5:6));
            d(i,1)=str2num(a(7:8));
            if d(i,1)~=d(i-1,1) | i==m
                s=s+1;
                ind(s,1)=i;
                TS(s-1,1)=y(i-1,1);
                TS(s-1,2)=mo(i-1,1);
                TS(s-1,3)=d(i-1,1);

                days=TmSric(i-1,1)/86400+datenum(2000,1,1)-datenum(y(i-1,1),1,1);
                if y(i-1,1)==2008 || y(i-1,1)==1996 ||y(i-1,1)==2000 ||y(i-1,1)==2004 ||y(i-1,1)==2012 || y(i-1,1)==2016
                    le=366;
                else
                    le=365;
                end
                TS(s-1,4)=days/le+TS(s-1,1);
                TS(s-1,5)=TmSric(i-1,1)/86400+datenum(2000,1,1);



                TS(s-1,6)=nanmean(SSH_ku(ind(s-1):ind(s)-1));
                TS(s-1,7)=nanmedian(SSH_ku(ind(s-1):ind(s)-1));
                TS(s-1,8)=std(SSH_ku(ind(s-1):ind(s)-1));

                TS(s-1,9)=nanmean(SSH_c(ind(s-1):ind(s)-1));
                TS(s-1,10)=nanmedian(SSH_c(ind(s-1):ind(s)-1));
                TS(s-1,11)=std(SSH_c(ind(s-1):ind(s)-1));

                TS(s-1,12)=nanmean(SSH_oce3_ku(ind(s-1):ind(s)-1));
                TS(s-1,13)=nanmedian(SSH_oce3_ku(ind(s-1):ind(s)-1));
                TS(s-1,14)=std(SSH_oce3_ku(ind(s-1):ind(s)-1));

                TS(s-1,12)=nanmean(SSH_oce3_c(ind(s-1):ind(s)-1));
                TS(s-1,13)=nanmedian(SSH_oce3_c(ind(s-1):ind(s)-1));
                TS(s-1,14)=std(SSH_oce3_c(ind(s-1):ind(s)-1));

                TS(s-1,15)=nanmean(SSH_red3_ku(ind(s-1):ind(s)-1));
                TS(s-1,16)=nanmedian(SSH_red3_ku(ind(s-1):ind(s)-1));
                TS(s-1,17)=std(SSH_red3_ku(ind(s-1):ind(s)-1));

                TS(s-1,15)=nanmean(SSH_red3_c(ind(s-1):ind(s)-1));
                TS(s-1,16)=nanmedian(SSH_red3_c(ind(s-1):ind(s)-1));
                TS(s-1,17)=std(SSH_red3_c(ind(s-1):ind(s)-1));

                TS(s-1,18)=nanmean(SSH_ice3_ku(ind(s-1):ind(s)-1));
                TS(s-1,19)=nanmedian(SSH_ice3_ku(ind(s-1):ind(s)-1));
                TS(s-1,20)=std(SSH_ice3_ku(ind(s-1):ind(s)-1));

                TS(s-1,21)=nanmean(SSH_ice3_c(ind(s-1):ind(s)-1));
                TS(s-1,22)=nanmedian(SSH_ice3_c(ind(s-1):ind(s)-1));
                TS(s-1,23)=std(SSH_ice3_c(ind(s-1):ind(s)-1));

            end
        end
        TS=TStcor(TS);
        TS_t=Outcor(TS,19,2.9,1);
        s=find(TS_t(:,20)>0.7);
        s0=find(TS_t(:,20)<0.7);
        TS_t(s,20)=nanmean(TS_t(s0,20));
        TS_final=[TS_t(:,5) TS_t(:,19:20)];
    else
        TS_final=[];
        TS=[];
    end


end
