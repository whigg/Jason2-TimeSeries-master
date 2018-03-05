clc;
clear all;
close all;

addpath('D:\jason2\gdr\s_gdr\cycle000');
mat={};
tp_track = 225;
% latmin = 82;
% latmax = 83;
% lonmin = -67;
% lonmax = -66;
latmin = 65;
latmax = 66;
lonmin = -128;
lonmax = -127;
if tp_track<10
    t1='0';
else
    t1='';
end
if tp_track<100
    t2='0';
else
    t2='';
end
% id_tp_track=['_00' t1 t2 num2str(tp_track) '_'];

k=1;
for i=1:327
    B=[];
    s=1;
%     i

    if i<10
        i1='0';
    else
        i1='';
    end
    if i<100
        i2='0';
    else
        i2='';
    end
    id_tp_track=['JA2_GPS_2PdP000_225_20080710_211338_20080710_220951'];
    %progressbar(i/265, 0);
    if i<10
        nam1='0';
    else
        nam1='';
    end
    if i<100
        nn='0';
    else
        nn='';
    end
    nam2= [nn num2str(i)];
    nam=['cycle_' nam1 nam2];
    drctry=['D:\jason2\gdr\s_gdr\cycle000'];
    cd(drctry)
    d=dir;
    for j=3:length(d)
        fndr=regexp(d(j,1).name,id_tp_track,'match');
        if isempty(fndr)==0
            fnam=d(j,1).name;
            [header data]=read_jason2_PH_nc(fnam);

            [m n]=size(data);
            for j=1:m
                if data(j,3)<latmax & data(j,3)>latmin & data(j,2)<lonmax & data(j,2)>lonmin
                    B(s,:)=data(j,:);
                    s=s+1;
                end
            end
            if isempty(B)==0
                mat{k,1}= i;
                mat{k,2}=B;
                k=k+1;
            end
        end
    end
    fclose all;
end