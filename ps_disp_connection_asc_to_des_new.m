clear all

%******************************************************************
% connection of PSI estimates (displacements): 

ph_asc=load('res_ph_corr_asc');
ph_des=load('res_ph_corr_des');

date_asc=textread('Asc\list_asc.txt');
date_des=textread('Des\list_des.txt');


common_date_interval=ph_asc.common_date_interval;

date_asc=ph_asc.date_asc;
date_des=ph_des.date_des;

%*********************************************************************
%***  METHOD 2:

fprintf('\n connection of PSI estimates (asc w.r.t. des)....');

ps2=load('Des\ps_aoi_lonlat_refined_des.mat');
master_ix=ps2.master_ix;
des_master=date_des(master_ix);
[master_date,BI,CI]=intersect(des_master,common_date_interval);

for i=1:length(common_date_interval)
ph_asc_to_des(:,i)=ph_asc.res_ph_asc(:,i)-ph_asc.res_ph_asc(:,CI);
end

% % test:
% figure;
% plot(common_date_interval,ph_asc_to_des(70000,:),'r');
% hold on;plot(common_date_interval,ph_asc.res_ph_asc(70000,:),'b');
% hold on;plot(common_date_interval,ph_des.res_ph_des(70000,:),'k');
% plot(common_date_interval(CI),0,'*')
ph_transformed_asc=ph_asc_to_des;
save('ph_transformed_asc','ph_transformed_asc','date_asc');




