clear all
% close all

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
coeff_all=load('coeff_all.mat');
e_cap=load('e_cap.mat');
sigma0_2_cap=load('sigma0_2_cap.mat');
% sigma0_2_cap2=load('sigma0_2_cap2.mat');
sigma2_cap_e_cap=load('sigma2_cap_e_cap.mat');
sigma2_x_cap=load('sigma2_x_cap.mat');
sigma2_x_cap2=load('sigma2_x_cap2.mat');
test_result1_vec=load('test_result1_vec.mat');
test_result2_vec=load('test_result2_vec.mat');
test_outlier=load('test_outlier.mat');
test_outlier2=load('test_outlier2.mat');
res_test=load('res_test.mat');
rmse_test=load('rmse_test.mat');

coeff_all=coeff_all.coeff_all;
e_cap=e_cap.e_cap;
sigma0_2_cap=sigma0_2_cap.sigma0_2_cap;
% sigma0_2_cap2=sigma0_2_cap2.sigma0_2_cap2;
sigma2_cap_e_cap=sigma2_cap_e_cap.sigma2_cap_e_cap;
sigma2_x_cap2=sigma2_x_cap.sigma2_x_cap;
test_result2_vec=test_result2_vec.test_result2_vec;
test_outlier2=test_outlier2.test_outlier2;
res_test=res_test.res_test;
rmse_test=rmse_test.rmse_test;
r_tilta2=e_cap./sqrt(sigma2_cap_e_cap);
%************************************************************
alpha=0.05;

indx=load('indx_all_new.mat');
indx=indx.indx_all;
% out_LS_index=[3 6 10 20 35 44 50 60 68];  
out_LS_index=[2 12 15 25 37 48 55 66 70];  


% r_tilta3_final=zeros(size(r_tilta2));
% test_outlier3_final=zeros(size(test_outlier2,1),size(test_outlier2,2)+1);
% test_outlier3_final=zeros(size(test_outlier2));
r_tilta3_final=r_tilta2;
test_outlier3_final=test_outlier2;
e_cap3_final=e_cap;
coeff_all3_final=coeff_all;
sigma0_2_cap3_final=sigma0_2_cap;
sigma2_x_cap3_final=sigma2_x_cap2;
sigma2_cap_e_cap3_final=sigma2_cap_e_cap;
index_No_final=zeros(size(e_cap));

for i=1:n
    max_r_tillta2(i)=max(abs(r_tilta2(i,:)));
end
thresh_r_tilta=mean(max_r_tillta2);


for i=1:n%78000:78020%1:length(indx)
    t=imagedate-imagedate(ps.master_ix);
    sigma0_2=1;
%     j=indx(i);
%     j=i;
    y0=ph_aoi_lonlat(i,:)*(-0.028/(2*pi));
    y0(out_LS_index)=[];
    t(out_LS_index)=[];

    y=y0';
    
    nu=length(t)-4;

    %l=ph_aoi_lonlat(i,:)*(-0.028/(2*pi));
    pl=eye(length(y));
    test_outlier3=test_outlier2(i,:);
    r_tilta3=r_tilta2(i,:);
    y3=y;
    t3=t;
    t3_final=t;
    k=0;
    index_No=(1:length(t));
    while (nnz(test_outlier3))~=length(y3) && (length(y3)>(length(y)-3))
        [val_max,ind_max]=max(abs(r_tilta3));
        if val_max > 2.2
            y3(ind_max)=[];
            t3(ind_max)=[];
            r_tilta3(ind_max)=[];
            index_No(ind_max)=[];
            pl=eye(length(y3));
            nu3=length(t3)-4;
            [coeff_all3,e_cap3,sigma0_2_cap3,sigma2_cap_e_cap3,sigma2_x_cap3]=adjustment_LS(t3,y3,pl,sigma0_2);
            e_cap3=e_cap3';
            sigma2_cap_e_cap3=sigma2_cap_e_cap3';
            sigma2_x_cap3=sigma2_x_cap3';
            [r_tilta3,test_outlier3]=test_residual(e_cap3,sigma2_cap_e_cap3,nu3,alpha);
            test_outlier3=test_outlier3';
            t3_final(ind_max)=0;
            k=k+1;
        else
            break;
        end
    end
    if k~=0
    No_of_outlier3_final(i)=k;
    coeff_all3_final(i,:)=coeff_all3;
    
    e_cap3_final(i,:)=0;
    e_cap3_final(i,1:length(e_cap3))=e_cap3;
    
    sigma0_2_cap3_final(i)=sigma0_2_cap3;
    
    sigma2_cap_e_cap3_final(i,:)=0;
    sigma2_cap_e_cap3_final(i,1:length(sigma2_cap_e_cap3))=sigma2_cap_e_cap3;
    
    sigma2_x_cap3_final(i,:)=0;
    sigma2_x_cap3_final(i,:)=sigma2_x_cap3;
    
    test_outlier3_final(i,:)=0;
    test_outlier3_final(i,1:length(test_outlier3))=test_outlier3;
    
    r_tilta3_final(i,:)=0;
    r_tilta3_final(i,1:length(r_tilta3))=r_tilta3;
    
    index_No_final(i,1:length(r_tilta3))=index_No;
    t3_final_allPSs(i,:)=t3_final;
    
%     figure(i+1000)
%     y_model=coeff_all3_final(i,1)*sin(2*pi*t3+coeff_all3_final(i,2))+coeff_all3_final(i,3)*t3+coeff_all3_final(i,4);
%     plot(t3,y_model,'r*');
%     hold on
%     plot(t3,y3,'ko');

    % Accuracy test:
%       l0=ph_aoi_lonlat(i,:)*(-0.028/(2*pi));
      x=imagedate-imagedate(ps.master_ix);
      x_test=x(out_LS_index);
      l_test=ph_aoi_lonlat(i,out_LS_index)*(-0.028/(2*pi));
      res_test(i,:)=(coeff_all3_final(i,1)*sin(2*pi*x_test+coeff_all3_final(i,2))+coeff_all3_final(i,3)*x_test+coeff_all3_final(i,4))'-l_test;
      rmse_test(i)=sqrt(res_test(i,:)*res_test(i,:)'/length(out_LS_index));
      plot(x_test,res_test(i,:),'ob--');hold on
    end
    if k==0
        t3_final_allPSs(i,:)=t3_final;
    end
end

% save('coeff_all3_final.mat','coeff_all3_final');
% 
% save('e_cap3_final.mat','e_cap3_final');
% 
% save('sigma0_2_cap3_final.mat','sigma0_2_cap3_final');
% 
% save('sigma2_cap_e_cap3_final.mat','sigma2_cap_e_cap3_final');
% 
% save('sigma2_x_cap3_final.mat','sigma2_x_cap3_final');
% 
% save('test_outlier3_final.mat','test_outlier3_final');
% 
% save('r_tilta3_final.mat','r_tilta3_final');
% 
% save('res_test3_final.mat','res_test');
% 
% save('rmse_test3_final.mat','rmse_test');
% 
% save('index_No_final.mat','index_No_final');
% 
% save('No_of_outlier3_final.mat','No_of_outlier3_final');

save('t3_final_allPSs.mat','t3_final_allPSs');
% 
% 
% 
% 
% 
% 
% 
% 
