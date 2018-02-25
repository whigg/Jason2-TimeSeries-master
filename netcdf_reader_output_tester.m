% function [ncData,ncid]=netcdf_reader(fnam)
clc;
close all;
clear all;

fnam = 'JA2_GPS_2PdP000_225_20080710_211338_20080710_220951.nc'
%cyc01  JA2_GPS_2PdP000_225_20080710_211338_20080710_220951.nc
%cyc02  JA2_GPS_2PdP001_225_20080720_191210_20080720_200823.nc
ncid = netcdf.open(fnam);

[numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);

for i=1:numvars
    sc=1;
    off=0;

    [varname, xtype, dimids, numatts]= netcdf.inqVar(ncid,i-1);

    varid = netcdf.inqVarID(ncid,varname);

    for j=1:numatts
        attname = netcdf.inqAttName(ncid,varid,j-1);

        attval = netcdf.getAtt(ncid,varid,attname);

        data = netcdf.getVar(ncid,varid);
        if strcmp(attname,'units')==1
            ncData{i,2}=attname;
            ncData{i,3}=attval;
        end
        if strcmp(attname,'scale_factor')==1
                        sc=attval;
        end
        if strcmp(attname,'add_offset')==1
                        off=attval;
        end
    end
    ncData{i,1}=varname;
%     ncData{i,4}=sc;
    ncData{i,4}=double(data)*sc+off;
end

s=1;
for i=1:179
    header{i,1}=ncData{i,1};
    [m,n]=size(ncData{i,4});

    if i==1
    c=m;
    end
    if m==c
    header{i,2}=n;
    data(:,s:s+n-1)=double(ncData{i,4});
    s=s+n;
    end
end
netcdf.close(ncid);
