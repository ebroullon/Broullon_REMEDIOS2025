close all
clear
load diffusivity_cruise_remedios.mat 
load LTd.mat

minS=min(min(S));
maxS=max(max(S));
minfluo=min(min(fluo));
maxfluo=max(max(fluo));

chl = 1.460*(fluo)-0.248; 
minchl = min(min(chl));
maxchl = max(max(chl));

i01 = find(date(:,3)==6,1,'last');
i02 = find(date(:,3)==9,1,'last');

dateV_I01 = date(1:i01,:);
dateV_I02 = date(i01+1:i02,:);
dateV_I03 = date(i02+1:end,:);

chl1 = chl(:,1:i01);
chl1 = [chl1,chl1(:,end)]; 
chl2 = chl(:,i01+1:i02);
chl2 = [chl2,chl2(:,end)];

chl3 = chl(:,i02+1:end);
chl3 = [chl3,chl3(:,end)];

T1 = T(:,1:i01);
T2 = T(:,i01+1:i02);
T3 = T(:,i02+1:end);

T1 = [T1,T1(:,end)];
T2 = [T2,T2(:,end)];
T3 = [T3,T3(:,end)];

S1 = S(:,1:i01);
S2 = S(:,i01+1:i02);
S3 = S(:,i02+1:end);

S1 = [S1,S1(:,end)]; 
S2 = [S2,S2(:,end)];
S3 = [S3,S3(:,end)];


LTd01=[NaN(1,712); LTd01];
LTd02=[NaN(1,259); LTd02];
LTd03=[NaN(1,687); LTd03];

LTd01=[LTd01; NaN(5,712)];
LTd02=[LTd02; NaN(5,259)];
LTd03=[LTd03; NaN(5,687)];

LT=[LTd01 LTd02 LTd03];

LTd01=[LTd01,LTd01(:,end)];
LTd02=[LTd02,LTd02(:,end)];
LTd03=[LTd03,LTd03(:,end)];

N21 = N2(:,1:i01);
N22 = N2(:,i01+1:i02);
N23 = N2(:,i02+1:end);

N21 = [N21,N21(:,end)]; 
N22 = [N22,N22(:,end)];
N23 = [N23,N23(:,end)];

dateN_I01=datenum(dateV_I01);
dateN_I02=datenum(dateV_I02);
dateN_I03=datenum(dateV_I03);
dateN_I01=[dateN_I01; dateN_I01(end,:)+0.005];
dateN_I02=[dateN_I02; dateN_I02(end,:)+0.005];
dateN_I03=[dateN_I03; dateN_I03(end,:)+0.005];

pres=pres(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure() %thorpe scale LT

dens1 = dens(:,1:i01); dens1 = [dens1 dens1(:,end)];
dens2 = dens(:,i01+1:i02); dens2 = [dens2 dens2(:,end)];
dens3 = dens(:,i02+1:end); dens3 = [dens3 dens3(:,end)];

vdens = [26.55 29];


set(gcf,'PaperUnits','centimeters');
xSize=21;
ySize=4.5;
xLeft=(21-xSize)/2;
yTop=(4.5-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[xLeft yTop xSize*50 ySize*50]);
set(gcf,'color',[1 1 1]);

xt=1.05/xSize;
yt=.5/ySize;
wt=(7.558-.5)/xSize; %18.4 total
ht=3.2/ySize;
axes('position',[xt yt wt ht]);

pcolor(dateN_I01,-pres,log10(LTd01));
shading flat;
datetick('x','dd');
xlim([dateN_I01(1,1)-0.005 dateN_I01(end)+0.005]);
ylim([-34.5 0]);

ylabel('pressure (dbar)');

map = cmocean('deep',10); map=flip(map);
colormap(map);
caxis([-2 0.5]);

hold on; 
[C,h]=contour(dateN_I01,-pres,dens1,vdens,'LineWidth',1.5,'LineColor','k');
clabel(C,'manual','color','k')
title('I01');

xt=(8.908-.5)/xSize; 
yt=.5/ySize;
wt=(2.749-.5)/xSize;
ht=3.2/ySize;
axes('position',[xt yt wt ht]);

pcolor(dateN_I02,-pres,log10(LTd02));
shading flat;
datetick('x','dd');
xlim([dateN_I02(1,1)-0.005 dateN_I02(end)+0.005]);
ylim([-34.5 0]);

set(gca,'Yticklabel',[]);

map = cmocean('deep',10); map=flip(map);
colormap(map);
caxis([-2 0.5]);

hold on; 
[C,h]=contour(dateN_I02,-pres,dens2,vdens,'LineWidth',1.5,'LineColor','k');
clabel(C,'manual','color','k')
title('I02');

xt=(11.957-1)/xSize; 
yt=.5/ySize;
wt=(7.293-.5)/xSize; 
ht=3.2/ySize;
axes('position',[xt yt wt ht]);

pcolor(dateN_I03,-pres,log10(LTd03));
shading flat;
datetick('x','dd');
xlim([dateN_I03(1,1)-0.005 dateN_I03(end)+0.005]);
ylim([-34.5 0]);

set(gca,'Yticklabel',[]);

hold on; 
[C,h]=contour(dateN_I03,-pres,dens3,vdens,'LineWidth',1.5,'LineColor','k');
clabel(C,'manual','color','k')
title('I03');

map = cmocean('deep',10); map=flip(map);
colormap(map);
caxis([-2 0.5]);

xc=(19.6-1.5)/xSize;
yc=.5/ySize;
wc=0.5/xSize;
hc=3.2/ySize;

c=colorbar;
label=ylabel(c,'log_{10}L_T (m)');
set(c,'Position',[xc yc wc hc]);

supersizeme(1.7);

%% N2 + isopycnals contour
minN2=min(min(N2));
maxN2=max(max(N2));

figure() 

set(gcf,'PaperUnits','centimeters');
xSize=21;
ySize=4.5;
xLeft=(21-xSize)/2;
yTop=(4.5-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[xLeft yTop xSize*50 ySize*50]);
set(gcf,'color',[1 1 1]);

xt=1.05/xSize;
yt=.5/ySize;
wt=(7.558-.5)/xSize; %18.4 total
ht=3.2/ySize;
axes('position',[xt yt wt ht]);

pcolor(dateN_I01,-pres,N21);
shading flat;
datetick('x','dd');
xlim([dateN_I01(1,1)-0.005 dateN_I01(end)+0.005]);
ylim([-34.5 0]);
ylabel('pressure (dbar)');

map = cmocean('thermal',20); 
colormap(map);
caxis([minN2 0.0040]);

hold on; 
[C,h]=contour(dateN_I01,-pres,dens1,vdens,'LineWidth',2,'LineColor',[.6 .6 .6]);
clabel(C,'manual','color',[.6 .6 .6])
title('I01');

xt=(8.908-.5)/xSize; 
yt=.5/ySize;
wt=(2.749-.5)/xSize;
ht=3.2/ySize;
axes('position',[xt yt wt ht]);

pcolor(dateN_I02,-pres,(N22));
shading flat;
datetick('x','dd');
xlim([dateN_I02(1,1)-0.005 dateN_I02(end)+0.005]);
ylim([-34.5 0]);

set(gca,'Yticklabel',[]);

map = cmocean('thermal',20); 
colormap(map);
caxis([minN2 0.0040]);

hold on; 
[C,h]=contour(dateN_I02,-pres,dens2,vdens,'LineWidth',2,'LineColor',[.6 .6 .6]);
clabel(C,'manual','color',[.6 .6 .6])
title('I02');

xt=(11.957-1)/xSize; 
yt=.5/ySize;
wt=(7.293-.5)/xSize; 
ht=3.2/ySize;
axes('position',[xt yt wt ht]);

pcolor(dateN_I03,-pres,(N23));
shading flat;
datetick('x','dd');
xlim([dateN_I03(1,1)-0.005 dateN_I03(end)+0.005]);
ylim([-34.5 0]);

set(gca,'Yticklabel',[]);

hold on; 
[C,h]=contour(dateN_I03,-pres,dens3,vdens,'LineWidth',2,'LineColor',[.6 .6 .6]);
clabel(C,'manual','color',[.6 .6 .6])
title('I03');

map = cmocean('thermal',20);
colormap(map);
caxis([minN2 0.0040]);

xc=(19.6-1.5)/xSize;
yc=.5/ySize;
wc=0.5/xSize;
hc=3.2/ySize;

c=colorbar;
label=ylabel(c,'N^2 (s^{-2})');
set(c,'Position',[xc yc wc hc]);

supersizeme(1.7);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure() %chl + isopycnals contour


set(gcf,'PaperUnits','centimeters');
xSize=21;
ySize=4.5;
xLeft=(21-xSize)/2;
yTop=(4.5-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[xLeft yTop xSize*50 ySize*50]);
set(gcf,'color',[1 1 1]);

xt=1.05/xSize;
yt=1/ySize;
wt=(7.558-.5)/xSize; 
ht=3.2/ySize;
axes('position',[xt yt wt ht]);

pcolor(dateN_I01,-pres,(chl1));
shading flat;
datetick('x','dd');
xlim([dateN_I01(1,1)-0.005 dateN_I01(end)+0.005]);
ylim([-34.5 0]);

ylabel('pressure (dbar)');
xlab=xlabel('July (days)');
set(xlab,'Position',[dateN_I01(end,1)+1 -39]);

map = brewermap(10,'Greens'); 
colormap(map);
caxis([0 3]);

hold on; 
[C,h]=contour(dateN_I01,-pres,dens1,vdens,'LineWidth',2,'LineColor','k');
clabel(C,'manual')

xt=(8.908-.5)/xSize; 
yt=1/ySize;
wt=(2.749-.5)/xSize; 
ht=3.2/ySize;
axes('position',[xt yt wt ht]);

pcolor(dateN_I02,-pres,(chl2));
shading flat;
datetick('x','dd');
xlim([dateN_I02(1,1)-0.005 dateN_I02(end)+0.005]);
ylim([-34.5 0]);

set(gca,'Yticklabel',[]);

caxis([0 3]);

hold on; 
[C,h]=contour(dateN_I02,-pres,dens2,vdens,'LineWidth',2,'LineColor','k');
clabel(C,'manual')

xt=(11.957-1)/xSize;
yt=1/ySize;
wt=(7.293-.5)/xSize; 
ht=3.2/ySize;
axes('position',[xt yt wt ht]);

pcolor(dateN_I03,-pres,(chl3));
shading flat;
datetick('x','dd');
xlim([dateN_I03(1,1)-0.005 dateN_I03(end)+0.005]);
ylim([-34.5 0]);

set(gca,'Yticklabel',[]);

caxis([0 3]);

hold on; 
[C,h]=contour(dateN_I03,-pres,dens3,vdens,'LineWidth',2,'LineColor','k');
clabel(C,'manual')

xc=(19.6-1.5)/xSize;
yc=1/ySize;
wc=0.5/xSize;
hc=3.2/ySize;

c=colorbar;
label=ylabel(c,'chlorophyll {\ita} (\mug L^-^1)');
set(c,'Position',[xc yc wc hc]);

supersizeme(1.7);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure() %T + isopycnals contour

set(gcf,'PaperUnits','centimeters');
xSize=21;
ySize=4.5;
xLeft=(21-xSize)/2;
yTop=(4.5-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[xLeft yTop xSize*50 ySize*50]);
set(gcf,'color',[1 1 1]);

xt=1.05/xSize;
yt=.5/ySize;
wt=(7.558-.5)/xSize; %18.4 total
ht=3.2/ySize;
axes('position',[xt yt wt ht]);

pcolor(dateN_I01,-pres,T1);
shading flat;
datetick('x','dd');
xlim([dateN_I01(1,1)-0.005 dateN_I01(end)+0.005]);
ylim([-34.5 0]);

ylabel('pressure (dbar)');

map = brewermap(20,'Spectral'); map=flip(map);
colormap(map);
caxis([13 22]);

hold on; 
[C,h]=contour(dateN_I01,-pres,dens1,vdens,'LineWidth',2,'LineColor','k');
clabel(C,'manual')
title('I01');

xt=(8.908-.5)/xSize; 
yt=.5/ySize;
wt=(2.749-.5)/xSize;
ht=3.2/ySize;
axes('position',[xt yt wt ht]);

pcolor(dateN_I02,-pres,(T2));
shading flat;
datetick('x','dd');
xlim([dateN_I02(1,1)-0.005 dateN_I02(end)+0.005]);
ylim([-34.5 0]);

set(gca,'Yticklabel',[]);

caxis([13 22]);

hold on; 
[C,h]=contour(dateN_I02,-pres,dens2,vdens,'LineWidth',2,'LineColor','k');
clabel(C,'manual')
title('I02');

xt=(11.957-1)/xSize; 
yt=.5/ySize;
wt=(7.293-.5)/xSize; 
ht=3.2/ySize;
axes('position',[xt yt wt ht]);

pcolor(dateN_I03,-pres,(T3));
shading flat;
datetick('x','dd');
xlim([dateN_I03(1,1)-0.005 dateN_I03(end)+0.005]);
ylim([-34.5 0]);

set(gca,'Yticklabel',[]);

hold on; 
[C,h]=contour(dateN_I03,-pres,dens3,vdens,'LineWidth',2,'LineColor','k');
clabel(C,'manual')
title('I03');

caxis([13 22]);

xc=(19.6-1.5)/xSize;
yc=.5/ySize;
wc=0.5/xSize;
hc=3.2/ySize;

c=colorbar;
label=ylabel(c,'temperature (ºC)');
set(c,'Position',[xc yc wc hc]);

supersizeme(1.7);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure() %S + isopycnals contour

set(gcf,'PaperUnits','centimeters');
xSize=21;
ySize=4.5;
xLeft=(21-xSize)/2;
yTop=(4.5-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[xLeft yTop xSize*50 ySize*50]);
set(gcf,'color',[1 1 1]);

xt=1.05/xSize;
yt=.5/ySize;
wt=(7.558-.5)/xSize; 
ht=3.2/ySize;
axes('position',[xt yt wt ht]);

pcolor(dateN_I01,-pres,S1);
shading flat;
datetick('x','dd');
xlim([dateN_I01(1,1)-0.005 dateN_I01(end)+0.005]);
ylim([-34.5 0]);

ylabel('pressure (dbar)');

map = cmocean('haline',120);
colormap(map);
caxis([34.3 35.7]);

hold on; 
[C,h]=contour(dateN_I01,-pres,dens1,vdens,'LineWidth',2,'LineColor','k');
clabel(C,'manual')
title('I01');

xt=(8.908-.5)/xSize; 
yt=.5/ySize;
wt=(2.749-.5)/xSize;
ht=3.2/ySize;
axes('position',[xt yt wt ht]);

pcolor(dateN_I02,-pres,(S2));
shading flat;
datetick('x','dd');
xlim([dateN_I02(1,1)-0.005 dateN_I02(end)+0.005]);
ylim([-34.5 0]);

set(gca,'Yticklabel',[]);

caxis([34.3 35.7]);

hold on; 
[C,h]=contour(dateN_I02,-pres,dens2,vdens,'LineWidth',2,'LineColor','k');
clabel(C,'manual')
title('I02');

xt=(11.957-1)/xSize; 
yt=.5/ySize;
wt=(7.293-.5)/xSize; 
ht=3.2/ySize;
axes('position',[xt yt wt ht]);

pcolor(dateN_I03,-pres,(S3));
shading flat;
datetick('x','dd');
xlim([dateN_I03(1,1)-0.005 dateN_I03(end)+0.005]);
ylim([-34.5 0]);

set(gca,'Yticklabel',[]);

hold on; 
[C,h]=contour(dateN_I03,-pres,dens3,vdens,'LineWidth',2,'LineColor','k');
clabel(C,'manual')
title('I03');

caxis([34.3 35.7]);

xc=(19.6-1.5)/xSize;
yc=.5/ySize;
wc=0.5/xSize;
hc=3.2/ySize;

c=colorbar;
label=ylabel(c,'salinity (PSU)');
set(c,'Position',[xc yc wc hc]);

supersizeme(1.7);

%%
figure()
%%%%%% nitrate %%%%%%
set(gcf,'PaperUnits','centimeters');
xSize=21;
ySize=4.5;
xLeft=(21-xSize)/2;
yTop=(4.5-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[xLeft yTop xSize*50 ySize*50]);
set(gcf,'color',[1 1 1]);

xt=1.05/xSize;
yt=1/ySize;
wt=(7.558-.5)/xSize; %1
ht=3.2/ySize;
axes('position',[xt yt wt ht]);

load nutrients_remedios_cruise.mat

i0F(1)=find(cruise==1,1,'first');
j0F(1)=find(cruise==1,1,'last');
i0F(2)=find(cruise==2,1,'first');
j0F(2)=find(cruise==2,1,'last');
i0F(3)=find(cruise==3,1,'first');
j0F(3)=find(cruise==3,1,'last');

time=nan(size(cruise));

for i=1:length(cruise)
    time(i)=(datenum(date(i,:))-datenum(date(i0F(cruise(i)),:)))*24;
end

%time
time1=time(i0F(1):j0F(1));
time2=time(i0F(2):j0F(2));
time3=time(i0F(3):j0F(3));

%date
date1=date(i0F(1):j0F(1),:);
date2=date(i0F(2):j0F(2),:);
date3=date(i0F(3):j0F(3),:);

%date matlab
date1_Vn=datenum(date1);
date2_Vn=datenum(date2);
date3_Vn=datenum(date3);

%Pres
pres1=pres(i0F(1):j0F(1));
pres2=pres(i0F(2):j0F(2));
pres3=pres(i0F(3):j0F(3));

%Nitrate
NO3_1=NO3(i0F(1):j0F(1));
NO3_2=NO3(i0F(2):j0F(2));
NO3_3=NO3(i0F(3):j0F(3));
    m_NO3=min(min(NO3));
    M_NO3=max(max(NO3));

kk_t1 = find(time1(2:end,1)-time1(1:end-1,1)~=0);
kk_t1=[kk_t1;length(time1)];
date_1= date(kk_t1,:);
date1_Vn_d=date1_Vn(kk_t1,:); 

kk_t2 = find(time2(2:end,1)-time2(1:end-1,1)~=0);
kk_t2=[kk_t2;length(time2)];
date_2= date2(kk_t2,:);
date2_Vn_d=date2_Vn(kk_t2,:); 

kk_t3 = find(time3(2:end,1)-time3(1:end-1,1)~=0);
kk_t3=[kk_t3;length(time3)];
date_3= date3(kk_t3,:);
date3_Vn_d=date3_Vn(kk_t3,:); 

date01 = datenum(date1(1,:));
dates01 = (datenum(date1(8,:)):1/2:datenum(date1(end,:)));
xtick_tick1 = (dates01-date01)*24;
xtick_dateL1 = datevec(dates01);
    
date02 = datenum(date2(1,:));
dates02 = (datenum(date2(1,:)):1/2:datenum(date2(end,:)));
xtick_tick2 = (dates02-date02)*24;
xtick_dateL2 = datevec(dates02);

date03 = datenum(date3(1,:));
dates03 = (datenum(date3(1,:)):1/2:datenum(date3(end,:)));
xtick_tick3 = (dates03-date03)*24;
xtick_dateL3 = datevec(dates03);

timeINT1=[0:4:max(time1),max(time1)];
timeINT2=[0:4:max(time2),max(time2)];
timeINT3=[0:4:max(time3),max(time3)];

dateINT1=[linspace(min(date1_Vn),max(date1_Vn),25)];
dateINT2=[linspace(min(date2_Vn),max(date2_Vn),4)];
dateINT3=[linspace(min(date3_Vn),max(date3_Vn),25)];

presINT=(1:2.5:35)';

[x1,y1]=meshgrid(timeINT1,presINT);
[x2,y2]=meshgrid(timeINT2,presINT);
[x3,y3]=meshgrid(timeINT3,presINT);
[x2,y2]=meshgrid(timeINT2,presINT);
[x3,y3]=meshgrid(timeINT3,presINT);

[x1d,y1d]=meshgrid(dateINT1,presINT);
[x2d,y2d]=meshgrid(dateINT2,presINT);
[x3d,y3d]=meshgrid(dateINT3,presINT);

NO3_INT_1=gridfit(time1(isfinite(NO3_1)),pres1(isfinite(NO3_1)),NO3_1(isfinite(NO3_1)),x1(1,:),y1(:,1));
NO3_INT_2=gridfit(time2(isfinite(NO3_2)),pres2(isfinite(NO3_2)),NO3_2(isfinite(NO3_2)),x2(1,:),y2(:,1));
NO3_INT_3=gridfit(time3(isfinite(NO3_3)),pres3(isfinite(NO3_3)),NO3_3(isfinite(NO3_3)),x3(1,:),y3(:,1));

NO3_INT_1=gridfit(date1_Vn(isfinite(NO3_1)),pres1(isfinite(NO3_1)),NO3_1(isfinite(NO3_1)),x1d(1,:),y1d(:,1));
NO3_INT_2=gridfit(date2_Vn(isfinite(NO3_2)),pres2(isfinite(NO3_2)),NO3_2(isfinite(NO3_2)),x2d(1,:),y2d(:,1));
NO3_INT_3=gridfit(date3_Vn(isfinite(NO3_3)),pres3(isfinite(NO3_3)),NO3_3(isfinite(NO3_3)),x3d(1,:),y3d(:,1));

NO3_INT_1(NO3_INT_1<0)=0;
NO3_INT_2(NO3_INT_2<0)=0;
NO3_INT_3(NO3_INT_3<0)=0;
NO3_INT=[NO3_INT_1,NO3_INT_2,NO3_INT_3];
m_NO3_INT=min(min(NO3_INT));
M_NO3_INT=max(max(NO3_INT));

dateINT1 = [dateINT1,dateINT1(1,end)];
NO3_INT_1 = [NO3_INT_1,NO3_INT_1(:,end)];


pcolor(dateINT1,-presINT,NO3_INT_1);
shading flat;
datetick('x','dd');
xlim([dateN_I01(1)-0.008 dateN_I01(end)+0.08]);
ylim([-34.5 0]);

ylabel('pressure (dbar)');
xlab=xlabel('July (days)');
set(xlab,'Position',[dateN_I01(end,1)+1 -39]);

map3=cmocean('matter');
colormap(map3);
caxis([m_NO3 M_NO3]);

hold on;
scatter(date1_Vn(isfinite(NO3_1)),-pres1(isfinite(NO3_1)),20,NO3_1(isfinite(NO3_1)),'o','filled','MarkerEdgeColor',[.6 .6 .6],'LineWidth',1);

pres=[1:1:45]';
hold on; 
[C,h]=contour(dateN_I01,-pres,dens1,vdens,'LineWidth',2,'LineColor','k');
clabel(C,'manual')


xt=(8.908-.5)/xSize;
yt=1/ySize;
wt=(2.749-.5)/xSize; 
ht=3.2/ySize;
axes('position',[xt yt wt ht]);

NO3_INT_2 = [NO3_INT_2,NO3_INT_2(:,end)];
dateINT2= [dateINT2,dateINT2(1,end)];

pcolor(dateINT2,-presINT,NO3_INT_2);
shading flat;
datetick('x','dd');
xlim([dateN_I02(1,1)-0.005 dateN_I02(end)+0.005]);
ylim([-34.5 0]);

set(gca,'Yticklabel',[]);
map3=cmocean('matter');
colormap(map3);
caxis([m_NO3 M_NO3]);

hold on;
scatter(date2_Vn(isfinite(NO3_2)),-pres2(isfinite(NO3_2)),20,NO3_2(isfinite(NO3_2)),'o','filled','MarkerEdgeColor',[.6 .6 .6],'LineWidth',1);

hold on; 
[C,h]=contour(dateN_I02,-pres,dens2,vdens,'LineWidth',2,'LineColor','k');
clabel(C,'manual')


xt=(11.957-1)/xSize; 
yt=1/ySize;
wt=(7.293-.5)/xSize; 
ht=3.2/ySize;
axes('position',[xt yt wt ht]);

NO3_INT_3 = [NO3_INT_3,NO3_INT_3(:,end)];
dateINT3 = [dateINT3,dateINT3(1,end)];

pcolor(dateINT3,-presINT,NO3_INT_3);
shading flat;
datetick('x','dd');
xlim([dateN_I03(1,1)-0.005 dateN_I03(end)+0.02]);
ylim([-34.5 0]);

set(gca,'Yticklabel',[]);
map3=cmocean('matter');
colormap(map3);
caxis([m_NO3 M_NO3]);

hold on;
scatter(date3_Vn(isfinite(NO3_3)),-pres3(isfinite(NO3_3)),20,NO3_3(isfinite(NO3_3)),'o','filled','MarkerEdgeColor',[.6 .6 .6],'LineWidth',1);

hold on; 
[C,h]=contour(dateN_I03,-pres,dens3,vdens,'LineWidth',2,'LineColor','k');
clabel(C,'manual')


xc=(19.6-1.5)/xSize;
yc=1/ySize;
wc=0.5/xSize;
hc=3.2/ySize;

c=colorbar;
label=ylabel(c,'nitrate (\muM)');
set(c,'Position',[xc yc wc hc]);

supersizeme(1.7);
