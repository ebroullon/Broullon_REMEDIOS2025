clear
load ctd_cruise_remedios_1.mat

idx1 = all(ismember(sample,'P01'),2);
Pidx1 = find(idx1);
idx2 = all(ismember(sample,'P02'),2);
Pidx2 = find(idx2);
idx3 = all(ismember(sample,'P03'),2);
Pidx3 = find(idx3);
idx4 = all(ismember(sample,'P04'),2);
Pidx4 = find(idx4);

pres=repmat(pres,1,366);
pres1=pres(:,Pidx1);
pres2=pres(:,Pidx2);
pres3=pres(:,Pidx3);
pres4=pres(:,Pidx4);


chl=chla;

chl_scm=NaN(142,366);
chl_scm_top=NaN(1,366);
chl_scm_bot=NaN(1,366);

for i=1:1:size(chla,2)
    
   if max(chla(:,i))>2 %subsurface chlorophyll maximum criteria
       chl_scm_top(1,i)=find(chla(:,i)>=2,1,'first');
       chl_scm_bot(1,i)=find(chla(:,i)>=2,1,'last');

   else
   end
end


for i=1:1:length(chl_scm_top)
    
    if ~isnan(chl_scm_top(i))>0
    chla(1:chl_scm_top(i)-1,i)=NaN;
    chla(chl_scm_bot(i)+1:end,i)=NaN;
    
    else
        chla(:,i)=NaN;
    end
end

idx_st_vigo=sort([find(st(1,:)>=1 & st(1,:)<=31), find(st(1,:)==111)]);
idx_st_shelf=sort([find(st(1,:)>=32 & st(1,:)<=71), find(st(1,:)==333)]);
idx_st_ponte=sort([find(st(1,:)>=72 & st(1,:)<=94), find(st(1,:)==222)]);

idxV1=intersect(Pidx1',idx_st_vigo); stV1=st(idxV1);
idxV2=intersect(Pidx2',idx_st_vigo); stV2=st(idxV2);
idxV3=intersect(Pidx3',idx_st_vigo); stV3=st(idxV3);
idxV4=intersect(Pidx4',idx_st_vigo); stV4=st(idxV4);

idxP1=intersect(Pidx1',idx_st_ponte);
idxP2=intersect(Pidx2',idx_st_ponte);
idxP3=intersect(Pidx3',idx_st_ponte);
idxP4=intersect(Pidx4',idx_st_ponte);

idxSH1=intersect(Pidx1',idx_st_shelf);
idxSH2=intersect(Pidx2',idx_st_shelf);
idxSH3=intersect(Pidx3',idx_st_shelf);
idxSH4=intersect(Pidx4',idx_st_shelf);

thickness=chl_scm_bot-chl_scm_top;
thickness1=thickness(Pidx1);
thickness2=thickness(Pidx2);
thickness3=thickness(Pidx3);
thickness4=thickness(Pidx4);

int_mean=mean(chla,1,'omitnan');
int_mean1=int_mean(Pidx1);
int_mean2=int_mean(Pidx2);
int_mean3=int_mean(Pidx3);
int_mean4=int_mean(Pidx4);

scm_mean_depth=(chl_scm_top+chl_scm_bot)/2;
scm_mean_depth1=scm_mean_depth(Pidx1);
scm_mean_depth2=scm_mean_depth(Pidx2);
scm_mean_depth3=scm_mean_depth(Pidx3);
scm_mean_depth4=scm_mean_depth(Pidx4);

int_meanV=int_mean(idx_st_vigo);
int_meanS=int_mean(idx_st_shelf);
int_meanP=int_mean(idx_st_ponte);

scm_mean_depthV=scm_mean_depth(idx_st_vigo);
scm_mean_detphS=scm_mean_depth(idx_st_shelf);
scm_mean_depthP=scm_mean_depth(idx_st_ponte);

%s01
int_mean1V=int_mean(idxV1);
int_mean1P=int_mean(idxP1);
int_mean1S=int_mean(idxSH1);
%S02
int_mean2V=int_mean(idxV2);
int_mean2P=int_mean(idxP2);
int_mean2S=int_mean(idxSH2);
%S03
int_mean3S=int_mean(idxSH3);
int_mean3V=int_mean(idxV3);
int_mean3P=int_mean(idxP3);
%S04
int_mean4V=int_mean(idxV4);
int_mean4P=int_mean(idxP4);
int_mean4S=int_mean(idxSH4);

%s01
scm_mean_depth1V=scm_mean_depth(idxV1);
scm_mean_depth1P=scm_mean_depth(idxP1);
scm_mean_depth1S=scm_mean_depth(idxSH1);
%s02
scm_mean_depth2V=scm_mean_depth(idxV2);
scm_mean_depth2P=scm_mean_depth(idxP2);
scm_mean_depth2S=scm_mean_depth(idxSH2);
%s03
scm_mean_depth3V=scm_mean_depth(idxV3);
scm_mean_depth3P=scm_mean_depth(idxP3);
scm_mean_depth3S=scm_mean_depth(idxSH3);
%S04
scm_mean_depth4V=scm_mean_depth(idxV4);
scm_mean_depth4P=scm_mean_depth(idxP4);
scm_mean_depth4S=scm_mean_depth(idxSH4);

%s01
thickness1V=thickness(idxV1);
thickness1P=thickness(idxP1);
thickness1S=thickness(idxSH1);
%S02
thickness2V=thickness(idxV2);
thickness2P=thickness(idxP2);
thickness2S=thickness(idxSH2);
%S03
thickness3S=thickness(idxSH3);
thickness3V=thickness(idxV3);
thickness3P=thickness(idxP3);
%S04
thickness4V=thickness(idxV4);
thickness4P=thickness(idxP4);
thickness4S=thickness(idxSH4);


%% S04
%intensity e depth 1matrix (1:2->V, 2:3->Shelf,4:5->P)
metricsV4=padconcatenation(scm_mean_depth4V,padconcatenation(int_mean4V,thickness4V,1),1);
metricsS4=padconcatenation(scm_mean_depth4S,padconcatenation(int_mean4S,thickness4S,1),1);
metricsP4=padconcatenation(scm_mean_depth4P,padconcatenation(int_mean4P,thickness4P,1),1);

%intensity & depth 1matrix SORTED by depth at each survey
metricsV4=sortrows(metricsV4',1)';
metricsS4=sortrows(metricsS4',1)';
metricsP4=sortrows(metricsP4',1)';
metrics4=padconcatenation(metricsV4,padconcatenation(metricsS4,metricsP4,1),1);

%%(1,4->V, 2,5->Shelf,3,6->P)
metrics4_separated=[];
metrics4_separated(1,:)=metrics4(1,:);%depth
metrics4_separated(2,:)=metrics4(4,:);%depth
metrics4_separated(3,:)=metrics4(7,:);%depth
metrics4_separated(4,:)=metrics4(2,:);%chl
metrics4_separated(5,:)=metrics4(5,:);%chl
metrics4_separated(6,:)=metrics4(8,:);%chl
metrics4_separated(7,:)=metrics4(3,:);%thick
metrics4_separated(8,:)=metrics4(6,:);%thick
metrics4_separated(9,:)=metrics4(9,:);%thick

%%BOXPLOTS S04
figure('color','w','Position',[1553 917 1100 300])  
ax1=subplot(1,3,1);
posi1=ax1.Position;
sup=suptitle('S04');
set(sup,'Position',[0.5000 -0.0300 0])
box1=violinplot(metrics4_separated(1:3,:)');
set(gca, 'XTickLabel',{'Vigo' 'Shelf' '  Pontevedra'});
ylabel('mean pressure (dbar)');
ylim([0 50]);

ax2=subplot(1,3,2);
posi2=ax2.Position;
box2=violinplot(metrics4_separated(4:6,:)');
set(gca, 'XTickLabel',{'Vigo' 'Shelf' '  Pontevedra'});
ylabel('mean chl {\ita} (µg L^-^1)');
ylim([0 5]);

ax3=subplot(1,3,3);
posi3=ax3.Position;
box3=violinplot(metrics4_separated(7:9,:)');
set(gca, 'XTickLabel',{'Vigo' 'Shelf' '  Pontevedra'});
ylabel('thickness (m)');
ylim([0 35]); 

supersizeme(1.6)

posf1=set(ax1,'Position',posi1);
posf2=set(ax2,'Position',posi2);
posf3=set(ax3,'Position',posi3);

%% S01
metricsV1=padconcatenation(scm_mean_depth1V,padconcatenation(int_mean1V,thickness1V,1),1);
metricsS1=padconcatenation(scm_mean_depth1S,padconcatenation(int_mean1S,thickness1S,1),1);
metricsP1=padconcatenation(scm_mean_depth1P,padconcatenation(int_mean1P,thickness1P,1),1);

metricsV1=sortrows(metricsV1',1)';
metricsS1=sortrows(metricsS1',1)';
metricsP1=sortrows(metricsP1',1)';
metrics1=padconcatenation(metricsV1,padconcatenation(metricsS1,metricsP1,1),1);

%%(1,4->V, 2,5->Shelf,3,6->P)
metrics1_separated=[];
metrics1_separated(1,:)=metrics1(1,:);%depth
metrics1_separated(2,:)=metrics1(4,:);%depth
metrics1_separated(3,:)=metrics1(7,:);%depth
metrics1_separated(4,:)=metrics1(2,:);%chl
metrics1_separated(5,:)=metrics1(5,:);%chl
metrics1_separated(6,:)=metrics1(8,:);%chl
metrics1_separated(7,:)=metrics1(3,:);%thick
metrics1_separated(8,:)=metrics1(6,:);%thick
metrics1_separated(9,:)=metrics1(9,:);%thick

%%BOXPLOTS S01
figure('color','w','Position',[1553 917 1100 300])  
ax4=subplot(1,3,1);
posi4=ax4.Position;
sup=suptitle('S01');
set(sup,'Position',[0.5000 -0.0300 0])
box4=violinplot(metrics1_separated(1:3,:)');
% xlabel('Location')
set(gca, 'XTickLabel',{'Vigo' 'Shelf' '  Pontevedra'});
ylabel('mean pressure (dbar)');
ylim([0 50]);

ax5=subplot(1,3,2);
posi5=ax5.Position;
box5=violinplot(metrics1_separated(4:6,:)'); 
set(gca, 'XTickLabel',{'Vigo' 'Shelf' '  Pontevedra'});
ylabel('mean chl {\ita} (µg L^-^1)');
ylim([0 5]);

ax6=subplot(1,3,3);
posi6=ax6.Position;
box6=violinplot(metrics1_separated(7:9,:)'); 
set(gca, 'XTickLabel',{'Vigo' 'Shelf' '  Pontevedra'});
ylabel('thickness (m)');
ylim([0 35]); 

supersizeme(1.6)

posf4=set(ax4,'Position',posi4);
posf5=set(ax5,'Position',posi5);
posf6=set(ax6,'Position',posi6);


%% S02
metricsV2=padconcatenation(scm_mean_depth2V,padconcatenation(int_mean2V,thickness2V,1),1);
metricsS2=padconcatenation(scm_mean_depth2S,padconcatenation(int_mean2S,thickness2S,1),1);
metricsP2=padconcatenation(scm_mean_depth2P,padconcatenation(int_mean2P,thickness2P,1),1);

metricsV2=sortrows(metricsV2',1)';
metricsS2=sortrows(metricsS2',1)';
metricsP2=sortrows(metricsP2',1)';
metrics2=padconcatenation(metricsV2,padconcatenation(metricsS2,metricsP2,1),1);

%%(1,4->V, 2,5->Shelf,3,6->P)
metrics2_separated=[];
metrics2_separated(1,:)=metrics2(1,:);%depth
metrics2_separated(2,:)=metrics2(4,:);%depth
metrics2_separated(3,:)=metrics2(7,:);%depth
metrics2_separated(4,:)=metrics2(2,:);%chl
metrics2_separated(5,:)=metrics2(5,:);%chl
metrics2_separated(6,:)=metrics2(8,:);%chl
metrics2_separated(7,:)=metrics2(3,:);%thick
metrics2_separated(8,:)=metrics2(6,:);%thick
metrics2_separated(9,:)=metrics2(9,:);%thick

%%BOXPLOTS S02
figure('color','w','Position',[1553 917 1100 300])  
ax7=subplot(1,3,1);
posi7=ax7.Position;
sup=suptitle('S02');
set(sup,'Position',[0.5000 -0.0300 0])
box7=violinplot(metrics2_separated([1:2 2],:)'); 


% xlabel('Location')
set(gca, 'XTickLabel',{'Vigo' 'Shelf' '  Pontevedra'});
ylabel('mean pressure (dbar)');
ylim([0 50]);

ax8=subplot(1,3,2);
posi8=ax8.Position;
box8=violinplot(metrics2_separated([4:5 5],:)'); 
set(gca, 'XTickLabel',{'Vigo' 'Shelf' '  Pontevedra'});
ylabel('mean chl {\ita} (µg L^-^1)');
ylim([0 5]);

ax9=subplot(1,3,3);
posi9=ax9.Position;
box9=violinplot(metrics2_separated([7:8 8],:)'); 
set(gca, 'XTickLabel',{'Vigo' 'Shelf' '  Pontevedra'});
ylabel('thickness (m)');
ylim([0 35]); 

supersizeme(1.6)

posf7=set(ax7,'Position',posi7);
posf8=set(ax8,'Position',posi8);
posf9=set(ax9,'Position',posi9);
%% S03
%intensidade e depth 1matrix (1:2->V, 2:3->Shelf,4:5->P)
metricsV3=padconcatenation(scm_mean_depth3V,padconcatenation(int_mean3V,thickness3V,1),1);
metricsS3=padconcatenation(scm_mean_depth3S,padconcatenation(int_mean3S,thickness3S,1),1);
metricsP3=padconcatenation(scm_mean_depth3P,padconcatenation(int_mean3P,thickness3P,1),1);

metricsV3=sortrows(metricsV3',1)';
metricsS3=sortrows(metricsS3',1)';
metricsP3=sortrows(metricsP3',1)';
metrics3=padconcatenation(metricsV3,padconcatenation(metricsS3,metricsP3,1),1);

%%(1,4->V, 2,5->Shelf,3,6->P)
metrics3_separated=[];
metrics3_separated(1,:)=metrics3(1,:);%depth
metrics3_separated(2,:)=metrics3(4,:);%depth
metrics3_separated(3,:)=metrics3(7,:);%depth
metrics3_separated(4,:)=metrics3(2,:);%chl
metrics3_separated(5,:)=metrics3(5,:);%chl
metrics3_separated(6,:)=metrics3(8,:);%chl
metrics3_separated(7,:)=metrics3(3,:);%thick
metrics3_separated(8,:)=metrics3(6,:);%thick
metrics3_separated(9,:)=metrics3(9,:);%thick

%%BOXPLOTS S03
figure('color','w','Position',[1553 917 1100 300])  
ax10=subplot(1,3,1);
posi10=ax10.Position;
sup=suptitle('S03');
set(sup,'Position',[0.5000 -0.0300 0])
box10=violinplot(metrics3_separated(1:3,:)'); 
set(gca, 'XTickLabel',{'Vigo' 'Shelf' '  Pontevedra'});
ylabel('mean pressure (dbar)');
ylim([0 50]);

ax11=subplot(1,3,2);
posi11=ax11.Position;
box11=violinplot(metrics3_separated(4:6,:)');
set(gca, 'XTickLabel',{'Vigo' 'Shelf' '  Pontevedra'});
ylabel('mean chl {\ita} (µg L^-^1)');
ylim([0 5]);

ax12=subplot(1,3,3);
posi12=ax12.Position;
box12=violinplot(metrics3_separated(7:9,:)');
set(gca, 'XTickLabel',{'Vigo' 'Shelf' '  Pontevedra'});
ylabel('thickness (m)');
ylim([0 35]); 


supersizeme(1.6)

posf10=set(ax10,'Position',posi10);
posf11=set(ax11,'Position',posi11);
posf12=set(ax12,'Position',posi12);

%% all toghether
metrics_all=[metrics1_separated metrics2_separated metrics3_separated metrics4_separated];

figure('color','w','Position',[1553 917 1100 300])  
ax10=subplot(1,3,1);
posi10=ax10.Position;
sup=suptitle('SCM');
set(sup,'Position',[0.5000 -0.0300 0])
box10=violinplot(metrics_all(1:3,:)'); 
set(gca, 'XTickLabel',{'Vigo' 'Shelf ' '  Pontevedra'});
ylabel('mean pressure (dbar)');
ylim([0 50]);

ax11=subplot(1,3,2);
posi11=ax11.Position;
box11=violinplot(metrics_all(4:6,:)'); 
set(gca, 'XTickLabel',{'Vigo' 'Shelf ' '  Pontevedra'});
ylabel('mean chl {\ita} (µg L^-^1)');
ylim([0 5]);

ax12=subplot(1,3,3);
posi12=ax12.Position;
box12=violinplot(metrics_all(7:9,:)');
set(gca, 'XTickLabel',{'Vigo' 'Shelf ' '  Pontevedra'});
ylabel('thickness (m)');
ylim([0 35]); 


supersizeme(1.6)

posf10=set(ax10,'Position',posi10);
posf11=set(ax11,'Position',posi11);
posf12=set(ax12,'Position',posi12);



