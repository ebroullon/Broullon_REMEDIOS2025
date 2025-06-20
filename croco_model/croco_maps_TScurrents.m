clear

load croco_output.mat

 %% averages current data every 24h since 24/06/2018 - 00:00 h

 usurfm=nan(115,190,45);
     for i=13:24:1069
        usurfm(:,:,i)=mean(usurf(:,:,i:(i+23)),3); %usurf mean 
     end
usurfm_23=usurfm(:,:,13:24:1069);

 vsurfm=nan(115,190,45);
     for i=13:24:1069
        vsurfm(:,:,i)=mean(vsurf(:,:,i:(i+23)),3); %vsurf 24 h mean
     end
vsurfm_23=vsurfm(:,:,13:24:1069);

tempm=nan(116,191,45);
     for i=13:24:1069
        tempm(:,:,i)=mean(temp(:,:,i:(i+23)),3); %usurf mean 
     end
tempm_23=tempm(:,:,13:24:1069);

 saltm=nan(116,191,45);
     for i=13:24:1069
        saltm(:,:,i)=mean(salt(:,:,i:(i+23)),3); %vsurf mean 
     end
saltm_23=saltm(:,:,13:24:1069);
 
temp_botm=nan(116,191,45);
     for i=13:24:1069
        temp_botm(:,:,i)=mean(temp_bot(:,:,i:(i+23)),3); %usurf mean 
     end
temp_botm_23=temp_botm(:,:,13:24:1069);

 salt_botm=nan(116,191,45);
     for i=13:24:1069
        salt_botm(:,:,i)=mean(salt_bot(:,:,i:(i+23)),3); %vsurf mean 
     end
salt_botm_23=salt_botm(:,:,13:24:1069);

ubotm=nan(115,190,45);
     for i=13:24:1069
        ubotm(:,:,i)=mean(ubot(:,:,i:(i+23)),3); %ubot mean 
     end
ubotm_23=ubotm(:,:,13:24:1069); %solo 23 capas

vbotm=nan(115,190,45);
     for i=13:24:1069
        vbotm(:,:,i)=mean(vbot(:,:,i:(i+23)),3); %vbot 24 h mean
     end
vbotm_23=vbotm(:,:,13:24:1069);

timeV_23=timevec(13:24:1069,:);

%% daily mean temp salt surf_bot

for i=23:1:45
    f=figure(i); set(gcf,'Color','w');
    f.Position=[1553 917 800 600]; 
    supt=suptitle(sprintf('%s',num2str(timeV_23(i,1)),'/',num2str(timeV_23(i,2)),'/',num2str(timeV_23(i,3)),' ',num2str(timeV_23(i,4)),':30'));
    set(supt,'Position',[0.5 -.03 0]);
    
    ax(1)=subplot(2,2,1);
    m_proj('Miller','lat',[42.05 42.44],'lon',[-9.065 -8.581]); 
    m_grid('tickdir','in','LineWidth',0.5,'color',[0.1 0.1 0.1],'xticklabel',[]);
    m_usercoast('costa_REMEDIOS','patch',[.5 .5 .5],'edgecolor','none');
    hold on;

    [~,h]=m_contourf(lon1,lat1,tempm_23(:,:,i),12:.01:22);
    set(h,'color','none')
    title('surface');
    ylabel('latitude');
    map1 = brewermap(25,'Spectral'); map1=flipud(map1);  
    colormap(ax(1),map1)
    posits1f=set(ax(1),'Position',[0.1800,0.5838,0.3347,0.3412]); %SUBPLOT FINAL POSITION
    caxis([15 21]);
    
    hold on;
    m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),usurfm_23(1:4:end,1:4:end,i),vsurfm_23(1:4:end,1:4:end,i),3,'LineWidth',1,'Color','k');

    ax(2)=subplot(2,2,2);
    m_proj('Miller','lat',[42.05 42.44],'lon',[-9.065 -8.581]); 
    m_grid('tickdir','in','LineWidth',0.5,'color',[0.1 0.1 0.1],'yticklabel',[],'xticklabel',[]);
    m_usercoast('costa_REMEDIOS','patch',[.5 .5 .5],'edgecolor','none');
    hold on;

    [~,h]=m_contourf(lon1,lat1,temp_botm_23(:,:,i),12:.01:22);
    set(h,'color','none')
    title('bottom');
    map2 = brewermap(25,'Spectral'); map2=flipud(map2); 
    colormap(ax(2),map2)
    cbar2=colorbar;
    posit2i=cbar2.Position; set(cbar2,'Position',[0.8339,0.5838,0.0357,0.3412]);
    posits2i=ax(2).Position; set(ax(2),'Position',[0.5003,0.5885,0.3347,0.3364]);
    posits2f=set(ax(2),'Position',[0.5003,0.5838,0.3347,0.3412]); %SUBPLOT FINAL POSITION
    cbar2.Label.String = 'temperature (^oC)';
    caxis([15 21]);
    
    hold on;
    m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),ubotm_23(1:4:end,1:4:end,i),vbotm_23(1:4:end,1:4:end,i),4,'LineWidth',1,'Color','k');    
    
    
    ax(3)=subplot(2,2,3);
    m_proj('Miller','lat',[42.05 42.44],'lon',[-9.065 -8.581]); 
    m_grid('tickdir','in','LineWidth',0.5,'color',[0.1 0.1 0.1]);
    m_usercoast('costa_REMEDIOS','patch',[.5 .5 .5],'edgecolor','none');
    hold on;

    [~,h]=m_contourf(lon1,lat1,saltm_23(:,:,i),25:.03:35.7);
    set(h,'color','none')
    title('surface');
    xlabel('longitude');
    ylabel('latitude');
    map3=cmocean('haline',25);
    colormap(ax(3),map3)
    posits3f=set(ax(3),'Position',[0.1800,0.1147,0.3347,0.3412]); 
    caxis([33.5 35.7]);
    
    hold on;
    m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),usurfm_23(1:4:end,1:4:end,i),vsurfm_23(1:4:end,1:4:end,i),4,'LineWidth',1,'Color','k');
    
    
    ax(4)=subplot(2,2,4);
    m_proj('Miller','lat',[42.05 42.44],'lon',[-9.065 -8.581]); 
    m_grid('tickdir','in','LineWidth',0.5,'color',[0.1 0.1 0.1],'yticklabel',[]);
    m_usercoast('costa_REMEDIOS','patch',[.5 .5 .5],'edgecolor','none');
    hold on;

    [~,h]=m_contourf(lon1,lat1,salt_botm_23(:,:,i),25:.03:35.7);
    set(h,'color','none')
    title('bottom');
    xlabel('longitude');
    map4=cmocean('haline',25);
    colormap(ax(4),map4)
    cbar4=colorbar;
    posit4i=cbar4.Position; set(cbar4,'Position',[0.8339,0.1147,0.0357,0.3412]);    
    posits4i=ax(4).Position; set(ax(4),'Position',[0.5003,0.1147,0.3347,0.3364]);
    posits4f=set(ax(4),'Position',[0.5003,0.1147,0.3347,0.3412]); %SUBPLOT FINAL POSITION
    cbar4.Label.String = 'salinity (psu)';
    caxis([33.5 35.7]);

    hold on;
    m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),ubotm_23(1:4:end,1:4:end,i),vbotm_23(1:4:end,1:4:end,i),4,'LineWidth',1,'Color','k');
    
    supersizeme(1.6); supersizeme(cbar2,1.2); supersizeme(cbar4,1.2)
    
    
    fpath = fullfile(pwd, 'Figures_4h');
        
    if ~exist(fpath, 'dir')
        mkdir(fpath);
    end

    saveas(gcf,fullfile(fpath,['TScurrents_',num2str(timeV_23(i,1)),'_',num2str(timeV_23(i,2)),'_',num2str(timeV_23(i,3)),'_',num2str(timeV_23(i,4))]),'png');
    close all
end




%% daily mean temp salt surface currents

for i=1:1:45 
    f=figure(i); set(gcf,'Color','w');
    f.Position=[1553 917 400 600]; 
    supt=suptitle(sprintf('%s',num2str(timeV_23(i,1)),'/',num2str(timeV_23(i,2)),'/',num2str(timeV_23(i,3)),' ',num2str(timeV_23(i,4)),':30'));
    set(supt,'Position',[0.42 -.03 0]);
    
    ax(1)=subplot(2,1,1);
    m_proj('Miller','lat',[42.05 42.44],'lon',[-9.065 -8.581]); 
    m_grid('tickdir','in','LineWidth',0.5,'color',[0.1 0.1 0.1],'xticklabel',[]);
    m_usercoast('costa_REMEDIOS','patch',[.5 .5 .5],'edgecolor','none');
    hold on;

    [~,h]=m_contourf(lon1,lat1,tempm_23(:,:,i),12:.01:22);
    set(h,'color','none')
    title('surface');
    ylabel('latitude');
    map1 = brewermap(25,'Spectral'); map1=flipud(map1);  
    colormap(ax(1),map1)
    posits1f=set(ax(1),'Position',[0.1800,0.5238,0.3347*2,0.3412]); %SUBPLOT FINAL POSITION
    caxis([15 21]);
    
    hold on;
    m_quiver(lon(1:6:end,1:6:end),lat(1:6:end,1:6:end),usurfm_23(1:6:end,1:6:end,i),vsurfm_23(1:6:end,1:6:end,i),3,'LineWidth',1,'Color','k');
    cbar2=colorbar;
    cbar2.Label.String = 'temperature (^oC)';
    caxis([15 21]);

    
    
    ax(3)=subplot(2,1,2);
    m_proj('Miller','lat',[42.05 42.44],'lon',[-9.065 -8.581]); 
    m_grid('tickdir','in','LineWidth',0.5,'color',[0.1 0.1 0.1]);
    m_usercoast('costa_REMEDIOS','patch',[.5 .5 .5],'edgecolor','none');
    hold on;

    [~,h]=m_contourf(lon1,lat1,saltm_23(:,:,i),25:.03:35.7);
    set(h,'color','none')
    xlabel('longitude');
    ylabel('latitude');

    colormap(ax(3),map3)
    posits3f=set(ax(3),'Position',[0.1800,0.1147,0.3347*2,0.3412]); %SUBPLOT FINAL POSITION
    caxis([33.5 35.7]);
    
    hold on;
    m_quiver(lon(1:6:end,1:6:end),lat(1:6:end,1:6:end),usurfm_23(1:6:end,1:6:end,i),vsurfm_23(1:6:end,1:6:end,i),3,'LineWidth',1,'Color','k');
   
    cbar4=colorbar;
    cbar4.Label.String = 'salinity (psu)';
%     set(cbar4,'YDir','reverse');
    caxis([33.5 35.7]);

    supersizeme(1.6); supersizeme(cbar2,1.2); supersizeme(cbar4,1.2)
    
    fpath = fullfile(pwd, 'Figures_6h');
        
    if ~exist(fpath, 'dir')
        mkdir(fpath);
    end

    saveas(gcf,fullfile(fpath,['TScurrents_',num2str(timeV_23(i,1)),'_',num2str(timeV_23(i,2)),'_',num2str(timeV_23(i,3)),'_',num2str(timeV_23(i,4))]),'png');
    close all
end
