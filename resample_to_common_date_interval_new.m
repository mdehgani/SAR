clear all


% %la_asc=load('la_red_asc');
% %la_des=load('la_red_des');
% %ph_asc=load('ph_sclam_grid_asc');
% %ph_des=load('ph_sclam_grid_des');
% 
% ph_asc=load('ph_grid_asc_coh');
% ph_des=load('ph_grid_des_coh');
% 
% la_asc=ph_asc.la_grid_asc;
% la_des=ph_des.la_grid_des;

% [m,n,p]=size(ph_asc.ph_grid_asc);

date_asc=textread('Asc\list_asc.txt');
date_des=textread('Des\list_des.txt');

month_day_asc=-(fix(date_asc/10000)*10000-date_asc);
year_asc=fix(date_asc/10000);
month_asc=fix(month_day_asc/100);
day_asc=-(fix(month_day_asc/100)*100-month_day_asc);
date_asc=year_asc+(month_asc-1)/12+day_asc/365;

month_day_des=-(fix(date_des/10000)*10000-date_des);
year_des=fix(date_des/10000);
month_des=fix(month_day_des/100);
day_des=-(fix(month_day_des/100)*100-month_day_des);
date_des=year_des+(month_des-1)/12+day_des/365;


all_date=[date_des' date_asc'];
[sorted_date,I]=sort(all_date);

min_date=max(date_asc(1),date_des(1));
max_date=min(date_asc(end),date_des(end));

date=find(sorted_date>=min_date & sorted_date<=max_date);

common_date_interval=all_date(I(date));

%*********************************************************************
% resample the displacement to common_date_interval:
fprintf('resampling to the common dates (ascending) ....');

ph_asc=load('Asc\ph_aoi_lonlat_refined_asc.mat');
ph_asc=ph_asc.ph_aoi_lonlat_refined;
n_ps_asc = size(ph_asc,1);

for i=1:n_ps_asc
ts=timeseries(ph_asc(i,:),date_asc');
res_ts=resample(ts,common_date_interval);
res_ph_asc(i,:)=res_ts.Data;
%plot(date_asc,shiftdim(ph_asc.ph_grid_asc(i,j,:),1));hold on;
%plot(common_date_interval,shiftdim(res_ph_grid_asc(i,j,:),1),'*k');
end

save('res_ph_corr_asc','res_ph_asc','common_date_interval','date_asc');

fprintf('resampling to the common dates (descending) ....');

ph_des=load('Des\ph_aoi_lonlat_refined_des.mat');
ph_des=ph_des.ph_aoi_lonlat_refined;
n_ps_des = size(ph_des,1);

for i=1:n_ps_des
ts=timeseries(ph_des(i,:),date_des');
res_ts=resample(ts,common_date_interval);
res_ph_des(i,:)=res_ts.Data;
%plot(date_asc,shiftdim(ph_asc.ph_grid_asc(i,j,:),1));hold on;
%plot(common_date_interval,shiftdim(res_ph_grid_asc(i,j,:),1),'*k');
end

save('res_ph_corr_des','res_ph_des','common_date_interval','date_des');



