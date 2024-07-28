clear all

vel_h=load('horz_vel_aoi');
vel_v=load('vert_vel_aoi');

disp_h=load('horz_disp_aoi');
disp_v=load('vert_disp_aoi');

figure(1)
scatter(vel_v.ps_aoi_asc_lonlat(:,1),vel_v.ps_aoi_asc_lonlat(:,2),5,vel_v.du_vel','filled')
title('vertical velcoity component');
colormap jet
colorbar

figure(2)
scatter(vel_h.ps_aoi_asc_lonlat(:,1),vel_h.ps_aoi_asc_lonlat(:,2),5,vel_h.de_vel','filled')
title('horizontal velocity component');
colormap jet
colorbar

figure(3)
scatter(disp_v.ps_aoi_asc_lonlat(:,1),disp_v.ps_aoi_asc_lonlat(:,2),5,disp_v.du(:,end)','filled')
title('vertical accumlated component');
colormap jet
colorbar

figure(4)
scatter(disp_h.ps_aoi_asc_lonlat(:,1),disp_h.ps_aoi_asc_lonlat(:,2),5,disp_h.de(:,end)','filled')
title('horizontal accumlated component');
colormap jet
colorbar

vel_v_forexcel=[vel_v.ps_aoi_asc_lonlat(:,1) vel_v.ps_aoi_asc_lonlat(:,2) vel_v.du_vel'];
vel_h_forexcel=[vel_h.ps_aoi_asc_lonlat(:,1) vel_h.ps_aoi_asc_lonlat(:,2) vel_h.de_vel'];
disp_v_forexcel=[disp_v.ps_aoi_asc_lonlat(:,1) disp_v.ps_aoi_asc_lonlat(:,2) disp_v.du(:,end)];
disp_h_forexcel=[disp_h.ps_aoi_asc_lonlat(:,1) disp_h.ps_aoi_asc_lonlat(:,2) disp_h.de(:,end)];

% xlswrite('vel_v.xlsx',vel_v_forexcel)
% xlswrite('vel_h.xlsx',vel_h_forexcel)
% xlswrite('disp_v.xlsx',disp_v_forexcel)
% xlswrite('disp_h.xlsx',disp_h_forexcel)


