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


la_asc=load('Asc\la_aoi_lonlat_refined_asc.mat');
la_des=load('Des\la_aoi_lonlat_refined_des.mat');
la_asc=la_asc.la_aoi_lonlat_refined*180/pi;
la_des=la_des.la_aoi_lonlat_refined*180/pi;


coh_asc=load('Asc\coh_aoi_lonlat_refined_asc.mat');
coh_des=load('Des\coh_aoi_lonlat_refined_des.mat');

coh_asc=coh_asc.coh_aoi_lonlat_refined;
coh_des=coh_des.coh_aoi_lonlat_refined;

common_date_interval=ph_des0.common_date_interval;

ps2_asc = load('Asc\ps_aoi_lonlat_refined_asc.mat');
ps2_des = load('Des\ps_aoi_lonlat_refined_des.mat');

sigma0_2_cap_asc=load('Asc\sigma0_2_cap_aoi_lonlat_refined_asc.mat');
sigma0_2_cap_asc=sigma0_2_cap_asc.sigma0_2_cap_aoi_lonlat_refined;

sigma0_2_cap_des=load('Des\sigma0_2_cap_aoi_lonlat_refined_des.mat');
sigma0_2_cap_des=sigma0_2_cap_des.sigma0_2_cap_aoi_lonlat_refined;


%***********************************************************************
% PS from des (ps2_des) corresponding to asc points within aoi are
% selected:
% aoi_lat = [29.8 30.05];
% aoi_lon = [52.65 53];
% 
% [i_des,j_des]=find(ps2_des.lonlat(:,1)>= aoi_lon(1) & ps2_des.lonlat (:,1) <=aoi_lon(2) & ps2_des.lonlat(:,2)>= aoi_lat(1) & ps2_des.lonlat (:,2) <=aoi_lat(2));
% [i_asc,j_asc]=find(ps2_asc.lonlat(:,1)>= aoi_lon(1) & ps2_asc.lonlat (:,1) <=aoi_lon(2) & ps2_asc.lonlat(:,2)>= aoi_lat(1) & ps2_asc.lonlat (:,2) <=aoi_lat(2));

ps_aoi_des_lonlat=ps2_des.ps_aoi_lonlat_refined;
ps_aoi_asc_lonlat=ps2_asc.ps_aoi_lonlat_refined;

ph_aoi_des=ph_des;
ph_aoi_asc=ph_asc;

la_aoi_des=la_des;
la_aoi_asc=la_asc;

coh_aoi_des=coh_des;
coh_aoi_asc=coh_asc;

n_asc=size(ph_aoi_asc,1);
n_des=size(ph_aoi_des,1);

%scatter(ps_aoi_asc_lonlat(:,1),ps_aoi_asc_lonlat(:,2),[],ph_aoi_asc(:,10),'filled')
% colormap jet



for i=1:n_asc
    index_des=dsearchn(ps_aoi_des_lonlat,ps_aoi_asc_lonlat(i,:));
    for k=1:length(common_date_interval)
       l=[ph_aoi_des(index_des,k);ph_aoi_asc(i,k)];
       A=[cosd(la_aoi_des(index_des)) -sind(la_aoi_des(index_des))*sind(a_des-270);...
      cosd(la_aoi_asc(i)) -sind(la_aoi_asc(i)*sind(a_asc-270))];
      Qy=[sigma0_2_cap_des(index_des) 0;0 sigma0_2_cap_asc(i)];
      d=inv(A'*inv(Qy)*A)*(A'*inv(Qy)*l);   
      du(i,k)=d(1,1)*(-0.028/2/pi);
      de(i,k)=d(2,1)*(-0.028/2/pi);
    end
    Q_cap_x=inv(A'*inv(Qy)*A);
    sigma_2_cap_x(i,:)=diag(Q_cap_x);
    
    % Estimation of sigma_2_cap_x using error propagation method:
    a1=cosd(la_aoi_asc(i))/(sind(la_aoi_asc(i))*sind(a_asc-270)*cosd(la_aoi_des(index_des))-sind(la_aoi_des(index_des))*sind(a_des-270)*cosd(la_aoi_asc(i)));
    b1=-cosd(la_aoi_des(index_des))/(sind(la_aoi_asc(i))*sind(a_asc-270)*cosd(la_aoi_des(index_des))-sind(la_aoi_des(index_des))*sind(a_des-270)*cosd(la_aoi_asc(i)));
    delta_de(i)=(a1^2*sigma0_2_cap_des(index_des)+b1^2*sigma0_2_cap_asc(i));
    
    a2=sind(la_aoi_asc(i))*sind(a_asc-270)/(cosd(la_aoi_des(index_des))*sind(la_aoi_asc(i))*sind(a_asc-270)-cosd(la_aoi_asc(i))*sind(la_aoi_des(index_des))*sind(a_des-270));
    b2=-sind(la_aoi_des(index_des))*sind(a_des-270)/(cosd(la_aoi_des(index_des))*sind(la_aoi_asc(i))*sind(a_asc-270)-cosd(la_aoi_asc(i))*sind(la_aoi_des(index_des))*sind(a_des-270));
    delta_du(i)=(a2^2*sigma0_2_cap_des(index_des)+b2^2*sigma0_2_cap_asc(i));
end
sigma_2_cap_x_calculation=[delta_du' delta_de'];

du_compare_precision=zeros(n_asc,length(common_date_interval));
de_compare_precision=zeros(n_asc,length(common_date_interval));

for i=1:n_asc
    for k=1:length(common_date_interval)
        if abs(du(i,k))>=sqrt(sigma_2_cap_x_calculation(i,1))*1.96
            du_compare_precision(i,k)=1;
        end
        if abs(de(i,k))>=sqrt(sigma_2_cap_x_calculation(i,2))*1.96
            de_compare_precision(i,k)=1;
        end
    end
end

figure(1)
title('vertical component');
scatter(ps_aoi_asc_lonlat(:,1),ps_aoi_asc_lonlat(:,2),5,du(:,end),'filled')
colormap jet
colorbar

figure(2)
title('east-west component');
scatter(ps_aoi_asc_lonlat(:,1),ps_aoi_asc_lonlat(:,2),5,de(:,end),'filled')
colormap jet
colorbar

figure(3);
scatter(ps_aoi_asc_lonlat(:,1),ps_aoi_asc_lonlat(:,2),5,sum(de_compare_precision,2),'filled')
colormap jet
colorbar

figure(4);
scatter(ps_aoi_asc_lonlat(:,1),ps_aoi_asc_lonlat(:,2),5,sum(du_compare_precision,2),'filled')
colormap jet
colorbar
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
save('horz_disp_aoi','de','ps_aoi_asc_lonlat');
save('vert_disp_aoi','du','ps_aoi_asc_lonlat');
save('de_compare_precision.mat','de_compare_precision');
save('du_compare_precision.mat','du_compare_precision');
save('sigma_2_cap_x.mat','sigma_2_cap_x');
save('sigma_2_cap_x_calculation.mat','sigma_2_cap_x_calculation');
save('common_date_interval.mat','common_date_interval');


