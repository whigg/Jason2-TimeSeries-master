[m,n]=size(TS);
d=datenum(TS(:,1),TS(:,2),TS(:,3));
u=unique(d);

TSc=nan(length(u),n);

for i=1:length(u)
    f=find(d==u(i));
    TSc(i,1:5)=TS(f(1),1:5);
    if length(f)>1
    TSc(i,4:n)=nanmean(TS(f,4:end));
    else
        TSc(i,4:n)=TS(f(1),4:end);
    end
end
