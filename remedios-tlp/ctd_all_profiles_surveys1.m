clear

load ctd_cruise_remedios_1.mat


idx1 = all(ismember(sample,'P01'),2); %survey S01
Pidx1 = find(idx1);
idx2 = all(ismember(sample,'P02'),2); %survey S02
Pidx2 = find(idx2);
idx3 = all(ismember(sample,'P03'),2); %survey S03
Pidx3 = find(idx3);
idx4 = all(ismember(sample,'P04'),2); %survey S04
Pidx4 = find(idx4);

pres=repmat(pres,1,366);


sigma(sigma<10)=NaN;

sigma1=sigma(:,Pidx1);
sigma2=sigma(:,Pidx2);
sigma3=sigma(:,Pidx3);
sigma4=sigma(:,Pidx4);

chl_S1=chla(:,Pidx1);
chl_S2=chla(:,Pidx2);
chl_S3=chla(:,Pidx3);
chl_S4=chla(:,Pidx4);

[chl_mean]=mean(chla,'omitnan');

chl_mean1=chl_mean(:,Pidx1);
chl_mean2=chl_mean(:,Pidx2);
chl_mean3=chl_mean(:,Pidx3);
chl_mean4=chl_mean(:,Pidx4);

[chl_max idx_max]=max(chla);

[chl_max1 idx_max1]=max(chl_S1);
[chl_max2 idx_max2]=max(chl_S2);
[chl_max3 idx_max3]=max(chl_S3);
[chl_max4 idx_max4]=max(chl_S4);

pres1=pres(:,Pidx1);
pres2=pres(:,Pidx2);
preS2=pres(:,Pidx3);
pres4=pres(:,Pidx4);

pres_max1=pres1(idx_max1);
pres_max2=pres2(idx_max2);
pres_max3=preS2(idx_max3);
pres_max4=pres4(idx_max4);

pres_max1(isnan(chl_max1))=NaN;
pres_max2(isnan(chl_max2))=NaN;
pres_max3(isnan(chl_max3))=NaN;
pres_max4(isnan(chl_max4))=NaN;

t1=t(:,Pidx1);
t2=t(:,Pidx2);
t3=t(:,Pidx3);
t4=t(:,Pidx4);

s1=s(:,Pidx1);
s2=s(:,Pidx2);
s3=s(:,Pidx3);
s4=s(:,Pidx4);

N21=N2_sort(:,Pidx1);
N22=N2_sort(:,Pidx2);
N23=N2_sort(:,Pidx3);
N24=N2_sort(:,Pidx4);

st1=st(Pidx1);
st2=st(Pidx2);
st3=st(Pidx3);
st4=st(Pidx4);

%% s4
%stations vigo st4(1:30,82) from #1 to #31 + #111
%stations shelf st4(31:63,83,84) from #32 to #64 + #333
%stations pontev st4(rest +222)

idx4V=[find(st4(1,:)>=1 & st4(1,:)<=31) find(st4(1,:)==111)];
st4V=st4(idx4V);
idx4S=[find(st4(1,:)>=32 & st4(1,:)<=65) find(st4(1,:)==333)];
st4S=st4(idx4S);
idx4P=setdiff(find(st4),[idx4V idx4S]);

idx4_sort=[idx4V idx4S idx4P];
st4_sort=st4(idx4_sort); 

chl_S4_sort=chl_S4(:,idx4_sort); chl_S4_sort(50,55)=NaN;
t_S4_sort=t4(:,idx4_sort); t_S4_sort(50,55)=NaN;
s_S4_sort=s4(:,idx4_sort); s_S4_sort(50,55)=NaN;
N2_S4_sort=N24(:,idx4_sort); N2_S4_sort(50,55)=NaN;
sigma_S4_sort=sigma4(:,idx4_sort); sigma_S4_sort(50,55)=NaN;

%%%%%%  SO4 %%%%%%
fig=figure('color','w','Position',[926 320 400 700]);

%%% chl s04 %%%
ax1=subplot(4,1,1); posi1=ax1.Position;
pcolor(chl_S4_sort);shading flat; hold on;
sup=suptitle('S04');
cmocean('algae');
set(gca, 'YDir','reverse')

patch([idx4S(1) idx4S(1)+length(idx4S) idx4S(1)+length(idx4S) idx4S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children')));

[C,h]=contour(sigma_S4_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color','k');
[C,h]=contour(sigma_S4_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color','k');

caxis([0.2 5.7]);
clabel(C,h,'color','k','fontsize',12);
c=colorbar;
c.Label.String='chlorophyll {\ita} (µg L^-^1)';
set(gca,'XTickLabel',[],'YTickLabel',[]);

xlim([.5 84.5]);
ylim([0 50]);

text((idx4V(1)+idx4V(end-1))/2,-6,'Vigo')
text((idx4S(1)+idx4S(end-2))/2-5,-6,'Shelf')
text((idx4P(1)+idx4P(end))/2-10,-6,'Pontevedra')

%%% t S04 %%%
ax2=subplot(4,1,2); posi2=ax2.Position;
pcolor(t_S4_sort);shading flat; hold on;
map=brewermap(25,'Spectral'); map=flipud(map);colormap(map)
set(gca, 'YDir','reverse')

patch([idx4S(1) idx4S(1)+length(idx4S) idx4S(1)+length(idx4S) idx4S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S4_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color',[0 0 0.6]);
[C,h]=contour(sigma_S4_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color',[0 0 0.6]);

caxis([12 20]); 
clabel(C,h,'color','k','fontsize',12);
c=colorbar;
c.Label.String='temperature (^oC)';
set(gca,'XTickLabel',[],'YTickLabel',[]);

xlim([.5 84.5]);
ylim([0 50]);

%%% salt s04 %%%
ax3=subplot(4,1,3); posi3=ax3.Position;
pcolor(s_S4_sort);shading flat; hold on;
cmocean('haline');
set(gca, 'YDir','reverse')

patch([idx4S(1) idx4S(1)+length(idx4S) idx4S(1)+length(idx4S) idx4S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S4_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color',[0 .3 0]);
[C,h]=contour(sigma_S4_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color',[0 .3 0]);

caxis([34.2 35.8]);
clabel(C,h,'color','k','fontsize',12);
c=colorbar;
c.Label.String='salinity (psu)';
set(gca,'XTickLabel',[],'YTickLabel',[]);

xlim([.5 84.5]);
ylim([0 50]);

%%% N2 s04 %%%
ax4=subplot(4,1,4); posi4=ax4.Position;
pcolor(N2_S4_sort);shading flat; hold on;
cmocean('thermal');
set(gca, 'YDir','reverse')

patch([idx4S(1) idx4S(1)+length(idx4S) idx4S(1)+length(idx4S) idx4S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S4_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color','y');
[C,h]=contour(sigma_S4_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color','y');

caxis([0 0.0039]);
clabel(C,h,'color','y','fontsize',12);
c=colorbar;
c.Label.String='N^2 (s^-^2)';
xlabel('CTD cast');
set(gca,'YTickLabel',[]);

xlim([.5 84.5]);
ylim([0 50]);


set(sup,'Position',[0.4500 -0.0200 0])
posf1=set(ax1,'Position',[.1300 .7573 .6500 .1800]);
posf2=set(ax2,'Position',[.1300    0.5381    .6500 .1800]);
posf3=set(ax3,'Position',[.1300    0.3205   .6500 .1800]);
posf4=set(ax4,'Position',[.1300    0.1014   .6500 .1800]);

supersizeme(1.6)

%% S3
idx3V=[find(st3(1,:)>=1 & st3(1,:)<=31) find(st3(1,:)==111)];
st3V=st3(idx3V);
idx3S=[find(st3(1,:)>=32 & st3(1,:)<=65) find(st3(1,:)==333)];
st3S=st3(idx3S);
idx3P=setdiff(find(st3),[idx3V idx3S]);
st3P=st3(idx3P);

idx3_sort=[idx3V idx3S idx3P];
st3_sort=st3(idx3_sort);

chl_S3_sort=chl_S3(:,idx3_sort);chl_S3_sort(63,41)=NaN;
t_S3_sort=t3(:,idx3_sort);t_S3_sort(63,41)=NaN;
s_S3_sort=s3(:,idx3_sort);s_S3_sort(63,41)=NaN;
N2_S3_sort=N23(:,idx3_sort);N2_S3_sort(63,41)=NaN;
sigma_S3_sort=sigma3(:,idx3_sort);sigma_S3_sort(63,41)=NaN;

%%%%%%  SO3 %%%%%%
fig3=figure('color','w','Position',[926 320 400 700]);
%%% chl s03 %%%
ax1=subplot(4,1,1); posi1=ax1.Position;
pcolor(chl_S3_sort);shading flat; hold on;
sup=suptitle('S03');
cmocean('algae');
set(gca, 'YDir','reverse')

patch([idx3S(1) idx3S(1)+length(idx3S) idx3S(1)+length(idx3S) idx3S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S3_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color','k');
[C,h]=contour(sigma_S3_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color','k');

caxis([0.2 5.7]);
clabel(C,h,'color','k','fontsize',12);
set(gca,'XTickLabel',[],'YTickLabel',[]);

xlim([.5 63.5]);
ylim([0 50]);

text((idx3V(1)+idx3V(end-1))/2,-6,'Vigo')
text((idx3S(1)+idx3S(end-2))/2-10,-6,'Shelf')
text((idx3P(1)+idx3P(end))/2-7,-6,'Pontevedra')

%%% t S03 %%%
ax2=subplot(4,1,2); posi2=ax2.Position;
pcolor(t_S3_sort);shading flat; hold on;
map=brewermap(25,'Spectral'); map=flipud(map);colormap(map)
set(gca, 'YDir','reverse')

patch([idx3S(1) idx3S(1)+length(idx3S) idx3S(1)+length(idx3S) idx3S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S3_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color',[0 0 0.6]);
[C,h]=contour(sigma_S3_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color',[0 0 0.6]);

caxis([12 20]); 
clabel(C,h,'color','k','fontsize',12);

set(gca,'XTickLabel',[],'YTickLabel',[]);

xlim([.5 63.5]);

ylim([0 50]);


%%% salt s03 %%%
ax3=subplot(4,1,3); posi3=ax3.Position;
pcolor(s_S3_sort);shading flat; hold on;
cmocean('haline');
set(gca, 'YDir','reverse')

patch([idx3S(1) idx3S(1)+length(idx3S) idx3S(1)+length(idx3S) idx3S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S3_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color',[0 .3 0]);
[C,h]=contour(sigma_S3_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color',[0 .3 0]);

caxis([34.2 35.8]);
clabel(C,h,'color','k','fontsize',12);
set(gca,'XTickLabel',[],'YTickLabel',[]);

xlim([.5 63.5]);
ylim([0 50]);


%%% N2 s03 %%%
ax4=subplot(4,1,4); posi4=ax4.Position;
pcolor(N2_S3_sort);shading flat; hold on;
cmocean('thermal');
set(gca, 'YDir','reverse')

patch([idx3S(1) idx3S(1)+length(idx3S) idx3S(1)+length(idx3S) idx3S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S3_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color','y');
[C,h]=contour(sigma_S3_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color','y');

caxis([0 0.0039]);
clabel(C,h,'color','y','fontsize',12);
xlabel('CTD cast');
set(gca,'YTickLabel',[]);

xlim([.5 63.5]);
ylim([0 50]);


set(sup,'Position',[0.45000 -0.0200 0]) 
posf1=set(ax1,'Position',[.1300 .7573 .6500 .1800]);
posf2=set(ax2,'Position',[.1300    0.5381    .6500 .1800]);
posf3=set(ax3,'Position',[.1300    0.3205   .6500 .1800]);
posf4=set(ax4,'Position',[.1300    0.1014   .6500 .1800]);

supersizeme(1.6)


%% s2
idx2V=[find(st2(1,:)>=1 & st2(1,:)<=31) find(st2(1,:)==111)];
st2V=st2(idx2V);
idx2S=[find(st2(1,:)>=32 & st2(1,:)<=65) find(st2(1,:)==333)];
st2S=st2(idx2S);
idx2P=setdiff(find(st2),[idx2V idx2S]);
st2P=st2(idx2P);

idx2_sort=[idx2V idx2S idx2P];
st2_sort=st2(idx2_sort);

chl_S2_sort=chl_S2(:,idx2_sort);
t_S2_sort=t2(:,idx2_sort);
s_S2_sort=s2(:,idx2_sort);
N2_S2_sort=N22(:,idx2_sort);
sigma_S2_sort=sigma2(:,idx2_sort);

%%%%%%  SO2 %%%%%%
fig2=figure('color','w','Position',[926 320 400 700]);

%%% chl s02 %%%
ax1=subplot(4,1,1); posi1=ax1.Position;
pcolor(chl_S2_sort);shading flat; hold on;
sup=suptitle('S02');
 
cmocean('algae');
set(gca, 'YDir','reverse')
set(gca,'XTickLabel',[],'YTickLabel',[]);


patch([idx2S(1) idx2S(1)+length(idx2S) idx2S(1)+length(idx2S) idx2S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S2_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color','k');
[C,h]=contour(sigma_S2_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color','k');

caxis([0.2 5.7]);
clabel(C,h,'color','k','fontsize',12);

xlim([.5 81.5]);
ylim([0 50]);


text((idx2V(1)+idx2V(end-1))/2,-6,'Vigo')
text((idx2S(1)+idx2S(end-2))/2-5,-6,'Shelf')
text((idx2P(1)+idx2P(end))/2-7,-6,'Pontevedra')

%%% t S02 %%%
ax2=subplot(4,1,2); posi2=ax2.Position;
pcolor(t_S2_sort);shading flat; hold on;
map=brewermap(25,'Spectral'); map=flipud(map);colormap(map)
set(gca, 'YDir','reverse')
set(gca,'XTickLabel',[],'YTickLabel',[]);

patch([idx2S(1) idx2S(1)+length(idx2S) idx2S(1)+length(idx2S) idx2S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S2_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color',[0 0 0.6]);
[C,h]=contour(sigma_S2_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color',[0 0 0.6]);

caxis([12 20]); 
clabel(C,h,'color','k','fontsize',12);

xlim([.5 81.5]);

ylim([0 50]);


%%% salt s02 %%%
ax3=subplot(4,1,3); posi3=ax3.Position;
pcolor(s_S2_sort);shading flat; hold on;
cmocean('haline');
set(gca, 'YDir','reverse')


patch([idx2S(1) idx2S(1)+length(idx2S) idx2S(1)+length(idx2S) idx2S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S2_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color',[0 .5 0]);
[C,h]=contour(sigma_S2_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-');

caxis([34.2 35.8]);
clabel(C,h,'color','k','fontsize',12);

set(gca,'XTickLabel',[],'YTickLabel',[]);

xlim([.5 81.5]);
ylim([0 50]);


%%% N2 s02 %%%
ax4=subplot(4,1,4); posi4=ax4.Position;
pcolor(N2_S2_sort);shading flat; hold on;
cmocean('thermal');
set(gca, 'YDir','reverse')

patch([idx2S(1) idx2S(1)+length(idx2S) idx2S(1)+length(idx2S) idx2S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S2_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color',[0 .5 0]);
[C,h]=contour(sigma_S2_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-');

caxis([0 0.0039]);
clabel(C,h,'color','y','fontsize',12);

c.Label.String='N^2 (s^-^2)';

xlabel('CTD cast');
set(gca,'YTickLabel',[]);

xlim([.5 81.5]);

ylim([0 50]);


set(sup,'Position',[0.45000 -0.0200 0]) 
posf1=set(ax1,'Position',[.1300 .7573 .6500 .1800]);
posf2=set(ax2,'Position',[.1300    0.5381    .6500 .1800]);
posf3=set(ax3,'Position',[.1300    0.3205   .6500 .1800]);
posf4=set(ax4,'Position',[.1300    0.1014   .6500 .1800]);

supersizeme(1.6)

%% s1
idx1V=[find(st1(1,:)>=1 & st1(1,:)<=31) find(st1(1,:)==111)];
st1V=st1(idx1V);
idx1S=[find(st1(1,:)>=32 & st1(1,:)<=65) find(st1(1,:)==333)];
st1S=st1(idx1S);
idx1P=setdiff(find(st1),[idx1V idx1S]);
st1P=st1(idx1P);

idx1_sort=[idx1V idx1S idx1P];
st1_sort=st1(idx1_sort);

chl_S1_sort=chl_S1(:,idx1_sort); chl_S1_sort(97:98,66)=NaN;chl_S1_sort(99,43)=NaN;
t_S1_sort=t1(:,idx1_sort);t_S1_sort(97:98,66)=NaN;t_S1_sort(99,43)=NaN;
s_S1_sort=s1(:,idx1_sort);s_S1_sort(97:98,66)=NaN;s_S1_sort(99,43)=NaN;
N2_S1_sort=N21(:,idx1_sort);N2_S1_sort(97:98,66)=NaN;N2_S1_sort(99,43)=NaN;
sigma_S1_sort=sigma1(:,idx1_sort); sigma_S1_sort(97:98,66)=NaN; sigma_S1_sort(99,43)=NaN;

%%%%%%  SO1 %%%%%%
fig1=figure('color','w','Position',[926 320 400 700]);
%%% chl s01 %%%
ax1=subplot(4,1,1); posi1=ax1.Position;
pcolor(chl_S1_sort);shading flat; hold on;
sup=suptitle('S01');

cmocean('algae');
set(gca, 'YDir','reverse')

patch([idx1S(1) idx1S(1)+length(idx1S) idx1S(1)+length(idx1S) idx1S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S1_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color','k');
[C,h]=contour(sigma_S1_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color','k');

caxis([0.2 5.7]);
clabel(C,h,'color','k','fontsize',12);


ylabel('Depth (m)');

set(gca,'XTickLabel',[])

xlim([.5 89.5]);

ylim([0 50]);

text((idx1V(1)+idx1V(end-1))/2,-6,'Vigo')
text((idx1S(1)+idx1S(end-2))/2-10,-6,'Shelf')
text((idx1P(1)+idx1P(end))/2-7,-6,'Pontevedra')

%%% t S01 %%%
ax2=subplot(4,1,2); posi2=ax2.Position;
pcolor(t_S1_sort);shading flat; hold on;
map=brewermap(25,'Spectral'); map=flipud(map);colormap(map)
set(gca, 'YDir','reverse')

patch([idx1S(1) idx1S(1)+length(idx1S) idx1S(1)+length(idx1S) idx1S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S1_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color',[0 0 0.6]);
[C,h]=contour(sigma_S1_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color',[0 0 0.6]);

caxis([12 20]); 
clabel(C,h,'color','k','fontsize',12);


ylabel('Depth (m)');

set(gca,'XTickLabel',[])

xlim([.5 89.5]);

ylim([0 50]);


%%% salt s01 %%%
ax3=subplot(4,1,3); posi3=ax3.Position;
pcolor(s_S1_sort);shading flat; hold on;
cmocean('haline');
set(gca, 'YDir','reverse')

patch([idx1S(1) idx1S(1)+length(idx1S) idx1S(1)+length(idx1S) idx1S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S1_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color',[0 .5 0]);
[C,h]=contour(sigma_S1_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-');

caxis([34.2 35.8]);
clabel(C,h,'color','k','fontsize',12);

ylabel('Depth (m)');

set(gca,'XTickLabel',[])

xlim([.5 89.5]);

ylim([0 50]);

%%% N2 s01 %%%
ax4=subplot(4,1,4); posi4=ax4.Position;
pcolor(N2_S1_sort);shading flat; hold on;
cmocean('thermal');
set(gca, 'YDir','reverse')

patch([idx1S(1) idx1S(1)+length(idx1S) idx1S(1)+length(idx1S) idx1S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S1_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color',[0 .5 0]);
[C,h]=contour(sigma_S1_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-');

caxis([0 0.0039]);
clabel(C,h,'color','y','fontsize',12);

ylabel('Depth (m)');
xlabel('CTD cast');

xlim([.5 89.5]);
ylim([0 50]);

set(sup,'Position',[0.45000 -0.0200 0]) 
posf1=set(ax1,'Position',[.1300 .7573 .6500 .1800]);
posf2=set(ax2,'Position',[.1300    0.5381    .6500 .1800]);
posf3=set(ax3,'Position',[.1300    0.3205   .6500 .1800]);
posf4=set(ax4,'Position',[.1300    0.1014   .6500 .1800]);

supersizeme(1.6)

% close all

%%
%%%%%%  SO1 %%%%%%
fig5=figure('color','w','Position',[926 320 400 700]);
%%% chl s01A %%%
ax1=subplot(8,1,1); posi1=ax1.Position;
pcolor(chl_S1_sort);shading flat; hold on;
sup=suptitle('S01');

cmocean('algae');
set(gca, 'YDir','reverse')

patch([idx1S(1) idx1S(1)+length(idx1S) idx1S(1)+length(idx1S) idx1S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S1_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color','k');
[C,h]=contour(sigma_S1_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color','k');

caxis([0.2 5.7]);
clabel(C,h,'color','k','fontsize',12);


ylabel('Depth (m)   ');

set(gca,'XTickLabel',[])
set(gca,'YTick',[0 10 20 30],'YTickLabel',[0 10 20 30]);

xlim([.5 89.5]);

ylim([-2 40]); 

text((idx1V(1)+idx1V(end-1))/2-1,-6,'Vigo')
text((idx1S(1)+idx1S(end-2))/2-15,-6,'Shelf')
text((idx1P(1)+idx1P(end))/2-7,-6,'Pontevedra')

%%% chl s01B %%%
ax12=subplot(8,1,2); posi12=ax12.Position;
pcolor(chl_S1_sort);shading flat; hold on;
sup=suptitle('S01');

cmocean('algae');
set(gca, 'YDir','reverse')

patch([idx1S(1) idx1S(1)+length(idx1S) idx1S(1)+length(idx1S) idx1S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S1_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color','k');
[C,h]=contour(sigma_S1_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color','k');

caxis([0.2 5.7]);
clabel(C,h,'color','k','fontsize',12);

set(gca,'XTickLabel',[])
set(gca,'YTick',[60 100],'YTickLabel',[60 100]);

xlim([.5 89.5]);
ylim([59 115]);

%%% t S01a %%%

ax2=subplot(8,1,3); posi2=ax2.Position;
pcolor(t_S1_sort);shading flat; hold on;
map=brewermap(25,'Spectral'); map=flipud(map);colormap(map)
set(gca, 'YDir','reverse')

patch([idx1S(1) idx1S(1)+length(idx1S) idx1S(1)+length(idx1S) idx1S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S1_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color',[0 0 0.6]);
[C,h]=contour(sigma_S1_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color',[0 0 0.6]);

caxis([12 20]); 
clabel(C,h,'color','k','fontsize',12);


ylabel('Depth (m)   ');

set(gca,'XTickLabel',[])
set(gca,'YTick',[0 10 20 30],'YTickLabel',[0 10 20 30]);

xlim([.5 89.5]);

ylim([-2 40]); 


%%% t S01b %%%

ax22=subplot(8,1,4); posi22=ax22.Position;
pcolor(t_S1_sort);shading flat; hold on;
map=brewermap(25,'Spectral'); map=flipud(map);colormap(map)
set(gca, 'YDir','reverse')

patch([idx1S(1) idx1S(1)+length(idx1S) idx1S(1)+length(idx1S) idx1S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S1_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color',[0 0 .6]);
[C,h]=contour(sigma_S1_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color',[0 0 .6]);

caxis([12 20]); 
clabel(C,h,'color','k','fontsize',12);
set(gca,'XTickLabel',[])
set(gca,'YTick',[60 100],'YTickLabel',[60 100]);


xlim([.5 89.5]);

ylim([59 115]);


%%% salt s01a %%%
ax3=subplot(8,1,5); posi3=ax3.Position;
pcolor(s_S1_sort);shading flat; hold on;
cmocean('haline');
set(gca, 'YDir','reverse')

patch([idx1S(1) idx1S(1)+length(idx1S) idx1S(1)+length(idx1S) idx1S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S1_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color',[0 .3 0]);
[C,h]=contour(sigma_S1_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color',[0 .3 0]);

caxis([34.2 35.8]);
clabel(C,h,'color','k','fontsize',12);

% c.Label.String='salinity (psu)';
ylabel('Depth (m)   ');

set(gca,'XTickLabel',[])
set(gca,'YTick',[0 10 20 30],'YTickLabel',[0 10 20 30]);

xlim([.5 89.5]);

ylim([-2 40]); 

%%% salt s01b %%%
ax32=subplot(8,1,6); posi32=ax32.Position;
pcolor(s_S1_sort);shading flat; hold on;
cmocean('haline');
set(gca, 'YDir','reverse')

patch([idx1S(1) idx1S(1)+length(idx1S) idx1S(1)+length(idx1S) idx1S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S1_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color',[0 .3 0]);
[C,h]=contour(sigma_S1_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color',[0 .3 0]);

caxis([34.2 35.8]);
clabel(C,h,'color','k','fontsize',12);
set(gca,'XTickLabel',[])
set(gca,'YTick',[60 100],'YTickLabel',[60 100]);

xlim([.5 89.5]);

ylim([59 115]);


%%% N2 s01A %%%
ax4=subplot(8,1,7); posi4=ax4.Position;
pcolor(N2_S1_sort);shading flat; hold on;
cmocean('thermal');
set(gca, 'YDir','reverse')

patch([idx1S(1) idx1S(1)+length(idx1S) idx1S(1)+length(idx1S) idx1S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S1_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color','y');
[C,h]=contour(sigma_S1_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color','y');

caxis([0 0.0039]);
clabel(C,h,'color','y','fontsize',12);

% c.Label.String='N^2 (s^-^2)';
ylabel('Depth (m)   ');

set(gca,'XTickLabel',[]);
set(gca,'YTick',[0 10 20 30],'YTickLabel',[0 10 20 30]);

xlim([.5 89.5]);
ylim([-2 40]); 

%%% N2 s01B %%%
ax42=subplot(8,1,8); posi42=ax42.Position;
pcolor(N2_S1_sort);shading flat; hold on;
cmocean('thermal');
set(gca, 'YDir','reverse')

patch([idx1S(1) idx1S(1)+length(idx1S) idx1S(1)+length(idx1S) idx1S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S1_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color','y');
[C,h]=contour(sigma_S1_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color','y');

caxis([0 0.0039]);
clabel(C,h,'color','y','fontsize',12);

% c.Label.String='N^2 (s^-^2)';

xlabel('CTD cast');
set(gca,'YTick',[60 100],'YTickLabel',[60 100]);

xlim([.5 89.5]);
ylim([59 115]);



set(sup,'Position',[0.45000 -0.0200 0]) 
posf1=set(ax1,'Position',[.1300 .7614 .6500 .1600]);
posf12=set(ax12,'Position',[.1300 .7254 .6500 .0300]);
posf2=set(ax2,'Position',[.1300    0.5414    .6500 .1600]);
posf22=set(ax22,'Position',[.1300    0.5054    .6500 .0300]);
posf3=set(ax3,'Position',[.1300    0.3214   .6500 .1600]);
posf32=set(ax32,'Position',[.1300    0.2854   .6500 .0300]);
posf4=set(ax4,'Position',[.1300    0.1014   .6500 .1600]);
posf42=set(ax42,'Position',[.1300    0.0654   .6500 .0300]);

supersizeme(1.6)

%%
%%%%%%  SO2 %%%%%%
fig6=figure('color','w','Position',[926 320 400 700]);
%%% chl s02A %%%
ax1=subplot(8,1,1); posi1=ax1.Position;
pcolor(chl_S2_sort);shading flat; hold on;
cmocean('algae');
set(gca, 'YDir','reverse')

patch([idx2S(1) idx2S(1)+length(idx2S) idx2S(1)+length(idx2S) idx2S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S2_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color','k');
[C,h]=contour(sigma_S2_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color','k');

caxis([0.2 5.7]);
clabel(C,h,'color','k','fontsize',12);
set(gca,'XTickLabel',[])
set(gca,'YTick',[0 10 20 30],'YTickLabel',[]);

xlim([.5 81.5]);
ylim([-2 40]); 

text((idx2V(1)+idx2V(end-1))/2-2,-6,'Vigo')
text((idx2S(1)+idx2S(end-2))/2-6,-6,'Shelf')
text((idx2P(1)+idx2P(end))/2-7,-6,'Pontevedra')

%%% chl s02B %%%
ax12=subplot(8,1,2); posi12=ax12.Position;
pcolor(chl_S2_sort);shading flat; hold on;
sup=suptitle('S02');
cmocean('algae');
set(gca, 'YDir','reverse')

patch([idx2S(1) idx2S(1)+length(idx2S) idx2S(1)+length(idx2S) idx2S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S2_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color','k');
[C,h]=contour(sigma_S2_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color','k');

caxis([0.2 5.7]);
clabel(C,h,'color','k','fontsize',12);

set(gca,'XTickLabel',[])
set(gca,'YTick',[60 100],'YTickLabel',[]);

xlim([.5 81.5]);
ylim([59 115]);

%%% t S02a %%%

ax2=subplot(8,1,3); posi2=ax2.Position;
pcolor(t_S2_sort);shading flat; hold on;
map=brewermap(25,'Spectral'); map=flipud(map);colormap(map)
set(gca, 'YDir','reverse')

patch([idx2S(1) idx2S(1)+length(idx2S) idx2S(1)+length(idx2S) idx2S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S2_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color',[0 0 .6]);
[C,h]=contour(sigma_S2_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color',[0 0 .6]);

caxis([12 20]); 
clabel(C,h,'color','k','fontsize',12);
set(gca,'XTickLabel',[])
set(gca,'YTick',[0 10 20 30],'YTickLabel',[]);

xlim([.5 81.5]);
ylim([-2 40]); 


%%% t S02b %%%
ax22=subplot(8,1,4); posi22=ax22.Position;
pcolor(t_S2_sort);shading flat; hold on;
map=brewermap(25,'Spectral'); map=flipud(map);colormap(map)
set(gca, 'YDir','reverse')

patch([idx2S(1) idx2S(1)+length(idx2S) idx2S(1)+length(idx2S) idx2S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S2_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color',[0 0 .6]);
[C,h]=contour(sigma_S2_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color',[0 0 .6]);

caxis([12 20]); 
clabel(C,h,'color','k','fontsize',12);
set(gca,'XTickLabel',[])
set(gca,'YTick',[60 100],'YTickLabel',[]);

xlim([.5 81.5]);
ylim([59 115]);


%%% salt s02a %%%
ax3=subplot(8,1,5); posi3=ax3.Position;
pcolor(s_S2_sort);shading flat; hold on;
cmocean('haline');
set(gca, 'YDir','reverse')

patch([idx2S(1) idx2S(1)+length(idx2S) idx2S(1)+length(idx2S) idx2S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S2_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color',[0 .3 0]);
[C,h]=contour(sigma_S2_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color',[0 .3 0]);

caxis([34.2 35.8]);
clabel(C,h,'color','k','fontsize',12);
set(gca,'XTickLabel',[])
set(gca,'YTick',[0 10 20 30],'YTickLabel',[]);

xlim([.5 81.5]);
ylim([-2 40]); 

%%% salt s02b %%%
ax32=subplot(8,1,6); posi32=ax32.Position;
pcolor(s_S2_sort);shading flat; hold on;
cmocean('haline');
set(gca, 'YDir','reverse')

patch([idx2S(1) idx2S(1)+length(idx2S) idx2S(1)+length(idx2S) idx2S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S2_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color',[0 .3 0]);
[C,h]=contour(sigma_S2_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color',[0 .3 0]);

caxis([34.2 35.8]);
clabel(C,h,'color','k','fontsize',12);
set(gca,'XTickLabel',[])
set(gca,'YTick',[60 100],'YTickLabel',[]);

xlim([.5 81.5]);
ylim([59 115]);


%%% N2 s01A %%%
ax4=subplot(8,1,7); posi4=ax4.Position;
pcolor(N2_S2_sort);shading flat; hold on;
cmocean('thermal');
set(gca, 'YDir','reverse')

patch([idx2S(1) idx2S(1)+length(idx2S) idx2S(1)+length(idx2S) idx2S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S2_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color','y');
[C,h]=contour(sigma_S2_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color','y');

caxis([0 0.0039]);
clabel(C,h,'color','y','fontsize',12);
set(gca,'XTickLabel',[]);
set(gca,'YTick',[0 10 20 30],'YTickLabel',[]);

xlim([.5 81.5]);
ylim([-2 40]); 

%%% N2 s02B %%%
ax42=subplot(8,1,8); posi42=ax42.Position;
pcolor(N2_S2_sort);shading flat; hold on;
cmocean('thermal');
set(gca, 'YDir','reverse')

patch([idx2S(1) idx2S(1)+length(idx2S) idx2S(1)+length(idx2S) idx2S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S2_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color','y');
[C,h]=contour(sigma_S2_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color','y');

caxis([0 0.0039]);
clabel(C,h,'color','y','fontsize',12);
xlabel('CTD cast');
set(gca,'YTick',[60 100],'YTickLabel',[]);

xlim([.5 81.5]);
ylim([59 115]);

set(sup,'Position',[0.45000 -0.0200 0]) 
posf1=set(ax1,'Position',[.1300 .7614 .6500 .1600]);
posf12=set(ax12,'Position',[.1300 .7254 .6500 .0300]);
posf2=set(ax2,'Position',[.1300    0.5414    .6500 .1600]);
posf22=set(ax22,'Position',[.1300    0.5054    .6500 .0300]);
posf3=set(ax3,'Position',[.1300    0.3214   .6500 .1600]);
posf32=set(ax32,'Position',[.1300    0.2854   .6500 .0300]);
posf4=set(ax4,'Position',[.1300    0.1014   .6500 .1600]);
posf42=set(ax42,'Position',[.1300    0.0654   .6500 .0300]);

supersizeme(1.6)

%%
%%%%%%  SO3 %%%%%%
fig7=figure('color','w','Position',[926 320 400 700]);
%%% chl s03A %%%
ax1=subplot(8,1,1); posi1=ax1.Position;
pcolor(chl_S3_sort);shading flat; hold on;
cmocean('algae');
set(gca, 'YDir','reverse')

patch([idx3S(1) idx3S(1)+length(idx3S) idx3S(1)+length(idx3S) idx3S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S3_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color','k');
[C,h]=contour(sigma_S3_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color','k');

caxis([0.2 5.7]);
clabel(C,h,'color','k','fontsize',12);
set(gca,'XTickLabel',[])
set(gca,'YTick',[0 10 20 30],'YTickLabel',[]);

xlim([.5 63.5]);
ylim([-2 40]); 

text((idx3V(1)+idx3V(end-1))/2,-6,'Vigo')
text((idx3S(1)+idx3S(end-2))/2-15,-6,'Shelf')
text((idx3P(1)+idx3P(end))/2-7,-6,'Pontevedra')

%%% chl s03B %%%
ax12=subplot(8,1,2); posi12=ax12.Position;
pcolor(chl_S3_sort);shading flat; hold on;
sup=suptitle('S03');
cmocean('algae');
set(gca, 'YDir','reverse')

patch([idx3S(1) idx3S(1)+length(idx3S) idx3S(1)+length(idx3S) idx3S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S3_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color','k');
[C,h]=contour(sigma_S3_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color','k');

caxis([0.2 5.7]);
clabel(C,h,'color','k','fontsize',12);

set(gca,'XTickLabel',[])
set(gca,'YTick',[60 100],'YTickLabel',[]);

xlim([.5 63.5]);
ylim([59 115]);

%%% t S02a %%%

ax2=subplot(8,1,3); posi2=ax2.Position;
pcolor(t_S3_sort);shading flat; hold on;
map=brewermap(25,'Spectral'); map=flipud(map);colormap(map)
set(gca, 'YDir','reverse')

patch([idx3S(1) idx3S(1)+length(idx3S) idx3S(1)+length(idx3S) idx3S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S3_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color',[0 0 0.6]);
[C,h]=contour(sigma_S3_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color',[0 0 0.6]);

caxis([12 20]); 
clabel(C,h,'color','k','fontsize',12);
set(gca,'XTickLabel',[])
set(gca,'YTick',[0 10 20 30],'YTickLabel',[]);

xlim([.5 63.5]);
ylim([-2 40]); 


%%% t S02b %%%
ax22=subplot(8,1,4); posi22=ax22.Position;
pcolor(t_S3_sort);shading flat; hold on;
map=brewermap(25,'Spectral'); map=flipud(map);colormap(map)
set(gca, 'YDir','reverse')

patch([idx3S(1) idx3S(1)+length(idx3S) idx3S(1)+length(idx3S) idx3S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S3_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color',[0 0 0.6]);
[C,h]=contour(sigma_S3_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color',[0 0 0.6]);

caxis([12 20]); 
clabel(C,h,'color','k','fontsize',12);
set(gca,'XTickLabel',[])
set(gca,'YTick',[60 100],'YTickLabel',[]);

xlim([.5 63.5]);
ylim([59 115]);


%%% salt s02a %%%
ax3=subplot(8,1,5); posi3=ax3.Position;
pcolor(s_S3_sort);shading flat; hold on;
cmocean('haline');
set(gca, 'YDir','reverse')

patch([idx3S(1) idx3S(1)+length(idx3S) idx3S(1)+length(idx3S) idx3S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S3_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color',[0 .3 0]);
[C,h]=contour(sigma_S3_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color',[0 .3 0]);

caxis([34.2 35.8]);
clabel(C,h,'color','k','fontsize',12);
set(gca,'XTickLabel',[])
set(gca,'YTick',[0 10 20 30],'YTickLabel',[]);

xlim([.5 63.5]);
ylim([-2 40]); 

%%% salt s03b %%%
ax32=subplot(8,1,6); posi32=ax32.Position;
pcolor(s_S3_sort);shading flat; hold on;
cmocean('haline');
set(gca, 'YDir','reverse')

patch([idx3S(1) idx3S(1)+length(idx3S) idx3S(1)+length(idx3S) idx3S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S3_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color',[0 .3 0]);
[C,h]=contour(sigma_S3_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color',[0 .3 0]);

caxis([34.2 35.8]);
clabel(C,h,'color','k','fontsize',12);
set(gca,'XTickLabel',[])
set(gca,'YTick',[60 100],'YTickLabel',[]);

xlim([.5 63.5]);
ylim([59 115]);


%%% N2 s03A %%%
ax4=subplot(8,1,7); posi4=ax4.Position;
pcolor(N2_S3_sort);shading flat; hold on;
cmocean('thermal');
set(gca, 'YDir','reverse')

patch([idx3S(1) idx3S(1)+length(idx3S) idx3S(1)+length(idx3S) idx3S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S3_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color','y');
[C,h]=contour(sigma_S3_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color','y');

caxis([0 0.0039]);
clabel(C,h,'color','y','fontsize',12);
set(gca,'XTickLabel',[]);
set(gca,'YTick',[0 10 20 30],'YTickLabel',[]);

xlim([.5 63.5]);
ylim([-2 40]); 

%%% N2 s03B %%%
ax42=subplot(8,1,8); posi42=ax42.Position;
pcolor(N2_S3_sort);shading flat; hold on;
cmocean('thermal');
set(gca, 'YDir','reverse')

patch([idx3S(1) idx3S(1)+length(idx3S) idx3S(1)+length(idx3S) idx3S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S3_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color','y');
[C,h]=contour(sigma_S3_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color','y');

caxis([0 0.0039]);
clabel(C,h,'color','y','fontsize',12);
xlabel('CTD cast');
set(gca,'YTick',[60 100],'YTickLabel',[]);

xlim([.5 63.5]);
ylim([59 115]);

set(sup,'Position',[0.45000 -0.0200 0]) 
posf1=set(ax1,'Position',[.1300 .7614 .6500 .1600]);
posf12=set(ax12,'Position',[.1300 .7254 .6500 .0300]);
posf2=set(ax2,'Position',[.1300    0.5414    .6500 .1600]);
posf22=set(ax22,'Position',[.1300    0.5054    .6500 .0300]);
posf3=set(ax3,'Position',[.1300    0.3214   .6500 .1600]);
posf32=set(ax32,'Position',[.1300    0.2854   .6500 .0300]);
posf4=set(ax4,'Position',[.1300    0.1014   .6500 .1600]);
posf42=set(ax42,'Position',[.1300    0.0654   .6500 .0300]);

supersizeme(1.6)

%%
%%%%%%  SO4 %%%%%%
fig8=figure('color','w','Position',[926 320 400 700]);
%%% chl s04A %%%
ax1=subplot(8,1,1); posi1=ax1.Position;
pcolor(chl_S4_sort);shading flat; hold on;
cmocean('algae');
set(gca, 'YDir','reverse')

patch([idx4S(1) idx4S(1)+length(idx4S) idx4S(1)+length(idx4S) idx4S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S4_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color','k');
[C,h]=contour(sigma_S4_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color','k');

caxis([0.2 5.7]);
clabel(C,h,'color','k','fontsize',12);
cc=colorbar;
cc.Label.String='chlorophyll {\ita} (µg L^-^1)';
set(gca,'XTickLabel',[])
set(gca,'YTick',[0 10 20 30],'YTickLabel',[]);

xlim([.5 84.5]);
ylim([-2 40]); 

text((idx4V(1)+idx4V(end-1))/2-5,-6,'Vigo')
text((idx4S(1)+idx4S(end-2))/2-5,-6,'Shelf')
text((idx4P(1)+idx4P(end))/2-7,-6,'Pontevedra')

%%% chl s04B %%%
ax12=subplot(8,1,2); posi12=ax12.Position;
pcolor(chl_S4_sort);shading flat; hold on;
sup=suptitle('S04');
cmocean('algae');
set(gca, 'YDir','reverse')

patch([idx4S(1) idx4S(1)+length(idx4S) idx4S(1)+length(idx4S) idx4S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S4_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color','k');
[C,h]=contour(sigma_S4_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color','k');

caxis([0.2 5.7]);
clabel(C,h,'color','k','fontsize',12);

set(gca,'XTickLabel',[])
set(gca,'YTick',[60 100],'YTickLabel',[]);

xlim([.5 84.5]);
ylim([59 115]);

%%% t S02a %%%

ax2=subplot(8,1,3); posi2=ax2.Position;
pcolor(t_S4_sort);shading flat; hold on;
map=brewermap(25,'Spectral'); map=flipud(map);colormap(map)
set(gca, 'YDir','reverse')

patch([idx4S(1) idx4S(1)+length(idx4S) idx4S(1)+length(idx4S) idx4S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S4_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color',[0 0 0.6]);
[C,h]=contour(sigma_S4_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color',[0 0 0.6]);

caxis([12 20]); 
clabel(C,h,'color','k','fontsize',12);
ct=colorbar;
ct.Label.String='temperature (^oC)';
set(gca,'XTickLabel',[])
set(gca,'YTick',[0 10 20 30],'YTickLabel',[]);

xlim([.5 84.5]);
ylim([-2 40]); 


%%% t S02b %%%
ax22=subplot(8,1,4); posi22=ax22.Position;
pcolor(t_S4_sort);shading flat; hold on;
map=brewermap(25,'Spectral'); map=flipud(map);colormap(map)
set(gca, 'YDir','reverse')

patch([idx4S(1) idx4S(1)+length(idx4S) idx4S(1)+length(idx4S) idx4S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S4_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color',[0 0 0.6]);
[C,h]=contour(sigma_S4_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color',[0 0 0.6]);

caxis([12 20]); 
clabel(C,h,'color','k','fontsize',12);
set(gca,'XTickLabel',[])
set(gca,'YTick',[60 100],'YTickLabel',[]);

xlim([.5 84.5]);
ylim([59 115]);


%%% salt s04a %%%
ax3=subplot(8,1,5); posi3=ax3.Position;
pcolor(s_S4_sort);shading flat; hold on;
cmocean('haline');
set(gca, 'YDir','reverse')

patch([idx4S(1) idx4S(1)+length(idx4S) idx4S(1)+length(idx4S) idx4S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S4_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color',[0 .3 0]);
[C,h]=contour(sigma_S4_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color',[0 .3 0]);

caxis([34.2 35.8]);
clabel(C,h,'color','k','fontsize',12);
cs=colorbar;
cs.Label.String='salinity (psu)';
set(gca,'XTickLabel',[])
set(gca,'YTick',[0 10 20 30],'YTickLabel',[]);

xlim([.5 84.5]);
ylim([-2 40]); 

%%% salt s04b %%%
ax32=subplot(8,1,6); posi32=ax32.Position;
pcolor(s_S4_sort);shading flat; hold on;
cmocean('haline');
set(gca, 'YDir','reverse')

patch([idx4S(1) idx4S(1)+length(idx4S) idx4S(1)+length(idx4S) idx4S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S4_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color',[0 .3 0]);
[C,h]=contour(sigma_S4_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color',[0 .3 0]);

caxis([34.2 35.8]);
clabel(C,h,'color','k','fontsize',12);
set(gca,'XTickLabel',[])
set(gca,'YTick',[60 100],'YTickLabel',[]);

xlim([.5 84.5]);
ylim([59 115]);


%%% N2 s04A %%%
ax4=subplot(8,1,7); posi4=ax4.Position;
pcolor(N2_S4_sort);shading flat; hold on;
cmocean('thermal');
set(gca, 'YDir','reverse')

patch([idx4S(1) idx4S(1)+length(idx4S) idx4S(1)+length(idx4S) idx4S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S4_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color','y');
[C,h]=contour(sigma_S4_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color','y');

caxis([0 0.0039]);
clabel(C,h,'color','y','fontsize',12);
cn2=colorbar;
cn2.Label.String='N^2 (s^-^2)';
set(gca,'XTickLabel',[]);
set(gca,'YTick',[0 10 20 30],'YTickLabel',[]);

xlim([.5 84.5]);
ylim([-2 40]); 

%%% N2 s04B %%%
ax42=subplot(8,1,8); posi42=ax42.Position;
pcolor(N2_S4_sort);shading flat; hold on;
cmocean('thermal');
set(gca, 'YDir','reverse')

patch([idx4S(1) idx4S(1)+length(idx4S) idx4S(1)+length(idx4S) idx4S(1)],[115 115 -2 -2],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'children',flipud(get(gca,'children'))); 

[C,h]=contour(sigma_S4_sort,[25 26 27 28],'LineWidth',.1,'LineStyle','-.','color','y');
[C,h]=contour(sigma_S4_sort,[26.4 26.8],'LineWidth',1,'LineStyle','-','color','y');

caxis([0 0.0039]);
clabel(C,h,'color','y','fontsize',12);
xlabel('CTD cast');
set(gca,'YTick',[60 100],'YTickLabel',[]);

xlim([.5 84.5]);
ylim([59 115]);

set(sup,'Position',[0.45000 -0.0200 0]) 
posf1=set(ax1,'Position',[.1300 .7614 .6500 .1600]);
posf12=set(ax12,'Position',[.1300 .7254 .6500 .0300]);
posf2=set(ax2,'Position',[.1300    0.5414    .6500 .1600]);
posf22=set(ax22,'Position',[.1300    0.5054    .6500 .0300]);
posf3=set(ax3,'Position',[.1300    0.3214   .6500 .1600]);
posf32=set(ax32,'Position',[.1300    0.2854   .6500 .0300]);
posf4=set(ax4,'Position',[.1300    0.1014   .6500 .1600]);
posf42=set(ax42,'Position',[.1300    0.0654   .6500 .0300]);
set(cc,'Position',[0.8050    0.7254    0.0500    0.1950]);
set(ct,'Position',[0.8050    0.5054    0.0500    0.1950]);
set(cs,'Position',[0.8050    0.2854    0.0500    0.1950]);
set(cn2,'Position',[0.8050    0.0654    0.0500    0.1950]);

supersizeme(1.6)