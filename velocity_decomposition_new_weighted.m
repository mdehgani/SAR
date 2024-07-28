clear all

%********************************************************************
% Decomposition of the displacement field:

% Vertical component

% projection of ALD(des) horz. disp.  into ALD(asc) horz.
% ALD horz. disp. : projection of horizontal displacement in azimuth look direction

%***********************************************************************
a_asc = -12.6457890+360;
a_des = -167.3693825+360;

% da=-166.380249+13.487083;

ph_asc0=load('ph_transformed_asc');
ph_des0=load('res_ph_corr_des');

ph_asc=ph_asc0.ph_transformed_asc;
ph_des=ph_des0.res_ph_des;

date_asc=ph_asc0.date_asc;
date_des=ph_des0.date_des;

la_asc=load('Asc_final\la2.mat');
la_des=load('Des_final\la2.mat');
% ph_asc2=load('ph_grid_asc_coh');
% ph_des2=load('ph_grid_des_coh');

la_asc=la_asc.la*180/pi;
la_des=la_des.la*180/pi;

% pm2_asc=load('Asc_final\pm2.mat');
% pm2_des=load('des_final\pm2.mat');
%coh_asc=pm2_asc.coh_ps;
%clear pm2_asc;
%coh_des=pm2_des.coh_ps;
%clear pm2_des;

coh_asc=load('Asc_final\coh_asc.mat');
coh_des=load('Des_final\coh_des.mat');

coh_asc=coh_asc.coh_ps;
coh_des=coh_des.coh_ps;

% vario_asc=load('variogram_parm_asc');
% vario_des=load('variogram_parm_des');

common_date_interval=ph_des0.common_date_interval;
% [m,n,p]=size(ph_asc.res_ph_grid_asc);

ps2_asc = load('Asc_final\ps2.mat');
ps2_des = load('Des_final\ps2.mat');
%***********************************************************************
% PS from des (ps2_des) corresponding to asc points within aoi are
% selected:
aoi_lat = [29.8 30.05];
aoi_lon = [52.65 53];

[i_des,j_des]=find(ps2_des.lonlat(:,1)>= aoi_lon(1) & ps2_des.lonlat (:,1) <=aoi_lon(2) & ps2_des.lonlat(:,2)>= aoi_lat(1) & ps2_des.lonlat (:,2) <=aoi_lat(2));
[i_asc,j_asc]=find(ps2_asc.lonlat(:,1)>= aoi_lon(1) & ps2_asc.lonlat (:,1) <=aoi_lon(2) & ps2_asc.lonlat(:,2)>= aoi_lat(1) & ps2_asc.lonlat (:,2) <=aoi_lat(2));

ps_aoi_des_lonlat=ps2_des.lonlat(i_des,:);
ps_aoi_asc_lonlat=ps2_asc.lonlat(i_asc,:);

ph_aoi_des=ph_des(i_des,:);
ph_aoi_asc=ph_asc(i_asc,:);

la_aoi_des=la_des(i_des,:);
la_aoi_asc=la_asc(i_asc,:);

coh_aoi_des=coh_des(i_des,:);
coh_aoi_asc=coh_asc(i_asc,:);


n_asc=size(ph_aoi_asc,1);
n_des=size(ph_aoi_des,1);

%scatter(ps_aoi_asc_lonlat(:,1),ps_aoi_asc_lonlat(:,2),[],ph_aoi_asc(:,10),'filled')
% colormap jet
%******************************************************************
% estimate the des velocity:
ps2=load('Des_final\ps2.mat');
master_ix=ps2.master_ix;
des_master=date_des(master_ix);
[master_date,BI,CI]=intersect(des_master,common_date_interval);

G=[ones(length(common_date_interval),1),(common_date_interval-common_date_interval(CI))'];

for i=1:n_des
    ph=ph_aoi_des(i,:)';
    m = G\ph;
    mean_v_des(i)=m(2,1)*(-0.028/2/pi);
end
figure(10)
title('mean_v_des');
scatter(ps_aoi_des_lonlat(:,1),ps_aoi_des_lonlat(:,2),5,mean_v_des','filled')
% hold on;plot(ps_aoi_asc_lonlat(i,1),ps_aoi_asc_lonlat(i,2),'k*')
colormap jet
colorbar


% estimate the asc velocity:
ps2=load('Asc_final\ps2.mat');
master_ix=ps2.master_ix;
asc_master=date_asc(master_ix);
[master_date,BI,CI]=intersect(asc_master,common_date_interval);

G=[ones(length(common_date_interval),1),(common_date_interval-common_date_interval(CI))'];

for i=1:n_asc
    ph=ph_aoi_asc(i,:)';
    m = G\ph;
    mean_v_asc(i)=m(2,1)*(-0.028/2/pi);
end
figure(20)
title('mean_v_asc');
scatter(ps_aoi_asc_lonlat(:,1),ps_aoi_asc_lonlat(:,2),5,mean_v_asc','filled')
% hold on;plot(ps_aoi_asc_lonlat(i,1),ps_aoi_asc_lonlat(i,2),'k*')
colormap jet
colorbar
%************************************************************************
for i=1:n_asc
    index_des=dsearchn(ps_aoi_des_lonlat,ps_aoi_asc_lonlat(i,:));
       l=[mean_v_des(index_des);mean_v_asc(i)];
       A=[cosd(la_aoi_des(index_des)) -sind(la_aoi_des(index_des))*sind(a_des-270);...
      cosd(la_aoi_asc(i)) -sind(la_aoi_asc(i)*sind(a_asc-270))];
  p=[coh_aoi_des(index_des) 0;0 coh_aoi_asc(i)];
      d=inv(A'*p*A)*(A'*p*l);   
      du_vel(i)=d(1,1);
      de_vel(i)=d(2,1);
    end
% 
figure(1)
title('vertical component');
scatter(ps_aoi_asc_lonlat(:,1),ps_aoi_asc_lonlat(:,2),5,du_vel','filled')
colormap jet
colorbar

figure(2)
title('east-west component');
scatter(ps_aoi_asc_lonlat(:,1),ps_aoi_asc_lonlat(:,2),5,de_vel','filled')
colormap jet
colorbar
% 
% figure(3)
% title('des los defo');
% scatter(ps_aoi_des_lonlat(:,1),ps_aoi_des_lonlat(:,2),5,ph_aoi_des(:,end),'filled')
% % hold on;plot(ps_aoi_des_lonlat(index_des,1),ps_aoi_des_lonlat(index_des,2),'k*')
% colormap jet
% colorbar
% 
% figure(4)
% title('asc los defo');
% scatter(ps_aoi_asc_lonlat(:,1),ps_aoi_asc_lonlat(:,2),5,ph_aoi_asc(:,end),'filled')
% % hold on;plot(ps_aoi_asc_lonlat(i,1),ps_aoi_asc_lonlat(i,2),'k*')
% colormap jet
% colorbar
% 
% 
save('horz_vel_aoi','de_vel','ps_aoi_asc_lonlat');
save('vert_vel_aoi','du_vel','ps_aoi_asc_lonlat');
save('mean_v_des_aoi','mean_v_des','ps_aoi_des_lonlat');
save('mean_v_asc_aoi','mean_v_asc','ps_aoi_asc_lonlat');

