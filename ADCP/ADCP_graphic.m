clear
close all

load ADCP_remedios.mat
MAB_matrix=repmat(MAB,6543,1)';
depth=Depth'-MAB_matrix;
depth=depth';


%%
figure (1)

theta=(0:5:90)*(pi/180);

for i=1:length(theta)
u_prima{i}=uf.*cos(theta(i))+vf.*sin(theta(i));
v_prima{i}=-uf.*sin(theta(i))+vf.*cos(theta(i));

var_uprima(i)=nanvar(u_prima{1,i}(:));
var_vprima(i)=nanvar(v_prima{1,i}(:));

figure(1)
plot(theta(i),var_uprima(i),'.r');
hold on;
plot(theta(i),var_vprima(i),'.b');

end

theta_sel=theta(find(var_uprima==max(var_uprima))); %

X_along=uf.*cos(theta_sel)+vf.*sin(theta_sel);
Y_transv=-uf.*sin(theta_sel)+vf.*cos(theta_sel);




%%
X_along=X_along.*100;
depthf=depthf';
timeV=datevec(time);
i01i=1645;i01f=2805;
i02i=3255;i02f=3630;
i03i=4067;i03f=5239;
time=repmat(time,70,1);


%% DETIDED EASTWARD VELOCITY 
X_original=X_along;
X_along=smooth2a(X_along,1,1);

figure('color','w')

set(gcf,'PaperUnits','centimeters');
xSize=21;
ySize=4.5;
xLeft=(21-xSize)/2;
yTop=(4.5-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[xLeft yTop xSize*50 ySize*50]);

xt=1.05/xSize;
yt=.5/ySize;
wt=(7.558-.5)/xSize; 
ht=3.2/ySize;
axes('position',[xt yt wt ht]);

contourf(time(1:48,i01i:i01f),-depthf(1:48,i01i:i01f),X_along(1:48,i01i:i01f),-(-15:1.5:15),'LineWidth',.3,'LineStyle','none'); %negativa porque asi e cara o este (fora da ria)
hold on;
contour(time(1:48,i01i:i01f),-depthf(1:48,i01i:i01f),X_along(1:48,i01i:i01f),-(-15:1.5:0),'LineWidth',.3,'LineStyle','-','color',[.7 .7 .7]); %positiv
contour(time(1:48,i01i:i01f),-depthf(1:48,i01i:i01f),X_along(1:48,i01i:i01f),-(0:1.5:15),'LineWidth',.3,'LineStyle',':','color',[.7 .7 .7]); %negativ
contour(time(1:48,i01i:i01f),-depthf(1:48,i01i:i01f),X_along(1:48,i01i:i01f),-(-30:30:30),'LineWidth',1,'LineStyle','-','color',[.4 .4 .4]); %linea 0
ylabel('depth (m)');
xlim([time(1,i01i) time(1,i01f)]);
datetick('x','dd/mmm','keeplimits');
map=brewermap(30,'RdBu');map=flip(map);
colormap(map);
caxis([-11 11]);
ylim([-35 0]) 
set(gca,'Xticklabel',[]);

xt=(8.908-.5)/xSize; 
yt=.5/ySize;
wt=(2.749-.5)/xSize; 
ht=3.2/ySize;
axes('position',[xt yt wt ht]);

contourf(time(1:48,i02i:i02f),-depthf(1:48,i02i:i02f),X_along(1:48,i02i:i02f),-(-15:1.5:15),'LineWidth',.3,'LineStyle','none');
hold on;
contour(time(1:48,i02i:i02f),-depthf(1:48,i02i:i02f),X_along(1:48,i02i:i02f),-(-15:1.5:0),'LineWidth',.3,'LineStyle','-','color',[.7 .7 .7]); %positiv
contour(time(1:48,i02i:i02f),-depthf(1:48,i02i:i02f),X_along(1:48,i02i:i02f),-(0:1.5:15),'LineWidth',.3,'LineStyle',':','color',[.7 .7 .7]); %negativ
contour(time(1:48,i02i:i02f),-depthf(1:48,i02i:i02f),X_along(1:48,i02i:i02f),-(-30:30:30),'LineWidth',1,'LineStyle','-','color',[.4 .4 .4]); %linea 0
xlim([time(1,i02i) time(1,i02f)]);
datetick('x','dd/mmm','keeplimits');
map=brewermap(30,'RdBu');map=flip(map);
colormap(map);
caxis([-11 11]);
ylim([-35 0]) 
set(gca,'Yticklabel',[]);set(gca,'Xticklabel',[]);


xt=(11.957-1)/xSize; %
yt=.5/ySize;
wt=(7.293-.5)/xSize; 
ht=3.2/ySize;
axes('position',[xt yt wt ht]);

contourf(time(1:48,i03i:i03f),-depthf(1:48,i03i:i03f),X_along(1:48,i03i:i03f),-(-15:1:15),'LineWidth',.3,'LineStyle','none');
hold on;
contour(time(1:48,i03i:i03f),-depthf(1:48,i03i:i03f),X_along(1:48,i03i:i03f),-(-15:1:0),'LineWidth',.3,'LineStyle','-','color',[.7 .7 .7]); %positiv
contour(time(1:48,i03i:i03f),-depthf(1:48,i03i:i03f),X_along(1:48,i03i:i03f),-(0:1:15),'LineWidth',.3,'LineStyle',':','color',[.7 .7 .7]); %negativ
contour(time(1:48,i03i:i03f),-depthf(1:48,i03i:i03f),X_along(1:48,i03i:i03f),-(-30:30:30),'LineWidth',1,'LineStyle','-','color',[.4 .4 .4]); %linea 0
xlim([time(1,i03i) time(1,i03f)]);
datetick('x','dd/mmm','keeplimits');
map=brewermap(30,'RdBu');map=flip(map);
colormap(map);
caxis([-11 11]);
ylim([-35 0]) 
set(gca,'Yticklabel',[]);
set(gca,'Xticklabel',[]);

xc=(19.6-1.5)/xSize;
yc=.5/ySize;
wc=0.5/xSize;
hc=3.2/ySize;

c=colorbar;
label=ylabel(c,{'u (cm s^-^1)'}); %eastward velocity
set(c,'Position',[xc yc wc hc]);

supersizeme(1.6)