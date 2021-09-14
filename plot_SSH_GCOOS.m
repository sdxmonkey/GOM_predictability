% plot SSH from GCOOS ssh data
clear
FS = 35;
file_in = 'ssh_gom_2016.nc';
irec = 117;%117;%282;

gmap=...
 [...
 10 50 120
 15 75 165
 30 110 200
 60 160 240
 80 180 250
130 210 255
160 240 255
200 250 255
255 255 255

255 255 255
255 232 120
255 192  60
255 160   0
255  96   0
255  50   0
225  20   0
192   0   0
165   0   0 ...
]./255;

sshlim = [-50 100];

tt = ncread(file_in,'time');
tnum = tt/3600/24+datenum('1970-1-1');
lon = ncread(file_in,'longitude');
lat = ncread(file_in,'latitude');

disp(['Reading ssh on ',datestr(tnum(irec))])
SSH = ncread(file_in,'ssh',[1 1 irec],[Inf Inf 1]);
SSH(SSH>1000) = NaN;

figure(1);clf
m_proj('miller','lon',[-98 -82],'lat',[18 31]);
m_pcolor(lon,lat,SSH');shading interp;colorbar vert;
hold on
m_contour(lon,lat,SSH',[17 17],'k','LineWidth',2);
m_gshhs_h('patch',[.8 .8 .8]);
caxis(sshlim)
colormap(gca,'jet')
m_grid('box','fancy','tickdir','out','fontsize',FS);
set(gca,'FontSize',FS)
title(datestr(tnum(irec)))
