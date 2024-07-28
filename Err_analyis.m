clear all
close all

%*************************************************************
% Load Data:
ps=load('ps2.mat');
ph=load('phuw2_corr_asc_flt_ref.mat');
la=load('la2.mat');
coh=load('coh_asc.mat');

%************************************************************
% To subset:
aoi_lat = [29.75 30.06];
aoi_lon = [52.7 53];

[i_des,j_des]=find(ps.lonlat(:,1)>= aoi_lon(1) & ps.lonlat (:,1) <=aoi_lon(2) & ps.lonlat(:,2)>= aoi_lat(1) & ps.lonlat (:,2) <=aoi_lat(2));
ps_aoi_lonlat=ps.lonlat(i_des,:);
ph_aoi_lonlat=ph.ph_all(i_des,:);
la_aoi_lonlat=la.la(i_des,:);
coh_aoi_lonlat=coh.coh_ps(i_des,:);

n=size(ph_aoi_lonlat,1); % Number of PS

%************************************************************
add=pwd;
add_file=strcat(add,'\list_asc.txt');
fid=fopen(add_file);
image_list=fscanf(fid,'%g',[1 inf]);
image_list=image_list';

day=image_list;
year=floor(day/10000);
month=floor((day-year*10000)/100);
monthday=day-year*10000-month*100;
imagedate=year+month/12+monthday/365;
imageday=datenum(year,month,monthday);

time_ts=ps.day-ps.day(ps.master_ix); % in terms of day
t=imagedate-imagedate(ps.master_ix); % in terms of year
%************************************************************
%*************************************************************
%No of outliers during elimation:
No_of_outlier=load('No_of_outlier3_final.mat');
No_of_outlier=No_of_outlier.No_of_outlier3_final;

rmse_test=load('rmse_test3_final.mat');
rmse_test_50=rmse_test.rmse_test;

rmse_test=load('rmse_test.mat');
rmse_test=rmse_test.rmse_test;

test_outlier=load('test_outlier2.mat');
test_outlier2=test_outlier.test_outlier2;

test_outlier_50=load('test_outlier3_final.mat');
test_outlier_50=test_outlier_50.test_outlier3_final;

index_of_NoOutlier=load('index_No_final.mat');
index_of_NoOutlier=index_of_NoOutlier.index_No_final;

e_cap=load('e_cap3_final.mat');
e_cap=e_cap.e_cap3_final;

sigma0_2_cap=load('sigma0_2_cap3_final.mat');
sigma0_2_cap=sigma0_2_cap.sigma0_2_cap3_final;


%***************************************************************
%No of outlier before elimination:
for i=1:n
   No_of_outlier_beforeEliminate(i)=size(test_outlier2,2)-nnz(test_outlier2(i,:));
end

%No of outliers remained after elimination:
for i=1:n
    vec_remained=test_outlier_50(i,1:nnz(index_of_NoOutlier(i,:)));
    No_of_outlier_remained(i)=length(vec_remained)-nnz(vec_remained);
end 

total_outlier_after_elimination=sum(No_of_outlier_remained);
total_outlier_before_elimination=sum(No_of_outlier_beforeEliminate);


% figure(100);
% scatter(ps_aoi_lonlat(:,1),ps_aoi_lonlat(:,2),2,No_of_outlier,'filled');colormap jet
% 
% [indx1,indx2]=find(No_of_outlier<10);
% ps_aoi_lonlat_good=ps_aoi_lonlat(indx2,:);
% No_of_outlier_good=No_of_outlier(indx2);
% figure(200)
% scatter(ps_aoi_lonlat_good(:,1),ps_aoi_lonlat_good(:,2),2,No_of_outlier_good,'filled');colormap jet

figure(300);
scatter(ps_aoi_lonlat(:,1),ps_aoi_lonlat(:,2),2,rmse_test_50,'filled');colormap jet

% figure(400);
% scatter(ps_aoi_lonlat(:,1),ps_aoi_lonlat(:,2),2,rmse_test,'filled');colormap jet

%******************************************************************************
% PS Elimination:
for i=1:n
    res(i)=e_cap(i,:)*e_cap(i,:)';
end
mean_res=mean(res);
std_res=std(res);

v=(res-mean_res)./std_res;
k=1;
for i=1:n
    if v(i)>-1.96 && v(i)<1.96
        ph_aoi_lonlat_refined(k,:)=ph_aoi_lonlat(i,:);
        ps_aoi_lonlat_refined(k,:)=ps_aoi_lonlat(i,:);
        la_aoi_lonlat_refined(k,:)=la_aoi_lonlat(i,:);
        coh_aoi_lonlat_refined(k,:)=coh_aoi_lonlat(i,:);
        sigma0_2_cap_aoi_lonlat_refined(k,:)=sigma0_2_cap(i);
        k=k+1;
    end
end
figure(500);
scatter(ps_aoi_lonlat_refined(:,1),ps_aoi_lonlat_refined(:,2),2,ph_aoi_lonlat_refined(:,end),'filled');colormap jet

figure(600);
scatter(ps_aoi_lonlat(:,1),ps_aoi_lonlat(:,2),2,ph_aoi_lonlat(:,end),'filled');colormap jet
bperp=ps.bperp;
day=ps.day;
ifgday=ps.ifgday;
ifgday_ix=ps.ifgday_ix;
ll0=ps.ll0;
master_day=ps.master_day;
master_ix=ps.master_ix;
mean_incidence=ps.mean_incidence;
mean_range=ps.mean_range;
n_ifg=ps.n_ifg;
n_image=ps.n_image;
n_ps=ps.n_ps;
xy=ps.xy;
ij=ps.ij;


save('ps_aoi_lonlat_refined_asc.mat','bperp','day','ifgday','ifgday_ix','ll0','master_day','master_ix','mean_incidence','mean_range','n_ifg','n_image','n_ps','xy','ij','ps_aoi_lonlat_refined');
save('ph_aoi_lonlat_refined_asc.mat','ph_aoi_lonlat_refined');
save('la_aoi_lonlat_refined_asc.mat','la_aoi_lonlat_refined');
save('coh_aoi_lonlat_refined_asc.mat','coh_aoi_lonlat_refined');
save('sigma0_2_cap_aoi_lonlat_refined_asc.mat','sigma0_2_cap_aoi_lonlat_refined');






