clear all
close all

%*************************************************************
% Load Data:
ps=load('ps2.mat');
ph=load('phuw2_corr_asc_flt_ref.mat');
coh=load('coh_asc.mat');

%************************************************************
% To subset:
aoi_lat = [29.75 30.06];
aoi_lon = [52.7 53];

[i_des,j_des]=find(ps.lonlat(:,1)>= aoi_lon(1) & ps.lonlat (:,1) <=aoi_lon(2) & ps.lonlat(:,2)>= aoi_lat(1) & ps.lonlat (:,2) <=aoi_lat(2));
ps_aoi_lonlat=ps.lonlat(i_des,:);
ph_aoi_lonlat=ph.ph_all(i_des,:);
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
alpha=0.05;

indx=load('indx_all_new.mat');
indx=indx.indx_all;
out_LS_index=[2 12 15 25 37 48 55 66 70];  
for i=1:n%90:100%length(indx)
    t=imagedate-imagedate(ps.master_ix);
    sigma0_2=coh_aoi_lonlat(i,:);
    
%     j=indx(i);
%     j=i;

    y0=ph_aoi_lonlat(i,:)*(-0.028/(2*pi));
    
    y0(out_LS_index)=[];
    t(out_LS_index)=[];
    y=y0';
    
    nu=length(t)-4;

    %l=ph_aoi_lonlat(i,:)*(-0.028/(2*pi));
    pl=eye(length(y));
    
    %Model solution with sigma0_2:
    [coeff_all(i,:),e_cap(i,:),sigma0_2_cap(i),sigma2_cap_e_cap(i,:),sigma2_x_cap(i,:)]=adjustment_LS(t,y,pl,sigma0_2);
    % sigma0 test:
    [y1,test_result1]=test_sigma2(sigma0_2,sigma0_2_cap(i),nu,alpha);
    test_result1_vec(i)=test_result1;
    y1_vec(i)=y1;
    % standardized residuals test:
    [r_tilta(i,:),test_outlier(i,:)]=test_residual(e_cap(i,:),sigma2_cap_e_cap(i,:),nu,alpha);
    
    % Model solution with sigma0_2_cap:
    pl2=inv(sigma0_2_cap(i)*pl);
    sigma0_2_2=1;
    [coeff_all2(i,:),e_cap2(i,:),sigma0_2_cap2(i),sigma2_cap_e_cap2(i,:),sigma2_x_cap2(i,:)]=adjustment_LS(t,y,pl2,sigma0_2_2);
    % sigma0 test:
    [y2,test_result2]=test_sigma2(sigma0_2_2,sigma0_2_cap2(i),nu,alpha);
    test_result2_vec(i)=test_result2;
    y2_vec(i)=y2;
    % standardized residuals test:
    [r_tilta2(i,:),test_outlier2(i,:)]=test_residual(e_cap2(i,:),sigma2_cap_e_cap2(i,:),nu,alpha);
% 
%     
%     figure(i+100)
%     y_model=coeff_all(i,1)*sin(2*pi*t+coeff_all(i,2))+coeff_all(i,3)*t+coeff_all(i,4);
%     plot(t,y_model,'r*');
%     hold on
%     plot(t,y,'ko');

      % Accuracy test:
      l0=ph_aoi_lonlat(i,:)*(-0.028/(2*pi));
      x=imagedate-imagedate(ps.master_ix);
      x_test=x(out_LS_index);
      l_test=ph_aoi_lonlat(i,out_LS_index)*(-0.028/(2*pi));
      res_test(i,:)=(coeff_all2(i,1)*sin(2*pi*x_test+coeff_all2(i,2))+coeff_all2(i,3)*x_test+coeff_all2(i,4))'-l_test;
      rmse_test(i)=sqrt(res_test(i,:)*res_test(i,:)'/length(out_LS_index));
%       plot(x_test,res_test(i,:));hold on
end

save('coeff_all.mat','coeff_all');
save('e_cap.mat','e_cap');
save('sigma0_2_cap.mat','sigma0_2_cap');
save('sigma0_2_cap2.mat','sigma0_2_cap2');
save('sigma2_cap_e_cap.mat','sigma2_cap_e_cap');
save('sigma2_x_cap.mat','sigma2_x_cap');
save('sigma2_x_cap2.mat','sigma2_x_cap2');
save('test_result1_vec.mat','test_result1_vec');
save('test_result2_vec.mat','test_result2_vec');
save('test_outlier.mat','test_outlier');
save('test_outlier2.mat','test_outlier2');
save('res_test.mat','res_test');
save('rmse_test.mat','rmse_test');








