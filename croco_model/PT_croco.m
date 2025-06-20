clear

load croco_output.mat


%% PT: power input to the turbulent motions
u=u(1:115,1:190,:);
v=v(1:115,1:190,:);
kb=.0025; %simpson & sharples pag 179
rho_0=1025; %reference density kg m-3

U=zeros(115,190,1105); U_cubic=zeros(115,190,1105); PT_time=zeros(115,190,1105);
for i=1:1:1105
    U(:,:,i)=(u(:,:,i).^2+v(:,:,i).^2).^.5; %m s-1
    U_cubic(:,:,i)=U(:,:,i).^3; %m3 s-3
    PT_time(:,:,i)=kb*rho_0*U_cubic(:,:,i); % W m-2 – stirring power of tidal flow 
end

PT_mean=mean(PT_time,3);
mask=mask(1:115,1:190);
PT_mean(mask==0)=NaN;

figure('color','w','Position',[1553 917 400 450])
pcolor(PT_mean');shading flat;
clim([0 0.01]);cmocean('thermal')
c=colorbar;
c.Label.String='mean P_T bottom (W m^-^2)';


%% PT/m

h=h(1:115,1:190);
pcolor(h');
shading flat

PT_m3=PT_mean./h;

%% 
figure('color','w','Position',[1553 917 700 500])

m_proj('Miller','lat',[42.05 42.64],'lon',[-9.115 -8.581]);
m_grid('tickdir','in','LineWidth',0.5,'color',[0.1 0.1 0.1]); 
m_usercoast('costa_REMEDIOS','patch',[.7 .7 .7],'edgecolor','none');
hold on;

[~,h2]=m_contourf(lon1(1:115,1:190),lat1(1:115,1:190),PT_m3,(0:0.00002:0.0012)); 

set(h2,'color','none')
cmocean('thermal')
c=colorbar;
clim([0 .0002]);
c.Label.String='mean P_T bottom (W m^-^3)'; 

m_scatter(-8.640021,42.433526,40,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'LineWidth',1);
m_scatter(-8.719693,42.239253,40,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'LineWidth',1);
m_scatter(-8.766947,42.597714,40,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'LineWidth',1);
m_scatter(-9.057725,42.774447,40,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'LineWidth',1);
m_scatter(-8.886364,42.785046,40,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'LineWidth',1);

m_text(-9.100,42.789,'Muros');
m_text(-8.873,42.785,'Noia');
m_text(-8.7450,42.450,'Pontevedra');
m_text(-8.700,42.239,'Vigo');
m_text(-8.749,42.597,'Vilagarcía de');
m_text(-8.749,42.570,'Arousa');

supersizeme(1.6)
