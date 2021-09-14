% plot mean SSH (LC location) and SSS to see if LC location is related with
% SSS gradient
clear
FS = 24;

dir_surf = 'D:\work\processed_surf\';
cases = {'1km_bulk_priv','1km_bulk','1km_bulk_nor','3p5km_bulk_priv'};
Nc = length(cases);
casenames = {'SP active river','SP passive river','SP no river','MR active river'};
Nms = [3,3,2];
Ipart = 3;
tnum0 = datenum('2016-7-6');
tnum1 = datenum('2016-7-6');

SSHlim = [-.9 .9];
Slim = [24 38];
Sglim = [0 0.5e-3];

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

tmp = imread('SSS_gradient.PNG');
cmap_Sg = double(squeeze(tmp(end:-1:1,12,:)))/255.0;

load grd_GOM1km.mat
fp1 = rnt_2grid(grd.f,'r','p');
ff1 = mean(fp1(:));
grdr1.lat = grd.lat_rho;
grdr1.lon = grd.lon_rho;
grdr1.mask = grd.mask_rho;
grd1 = grd;

load grd_GOM3p5km.mat
fp2 = rnt_2grid(grd.f,'r','p');
ff2 = mean(fp2(:));
grdr2.lat = grd.lat_rho;
grdr2.lon = grd.lon_rho;
grdr2.mask = grd.mask_rho;
grd2 = grd;

SSH_crit = 0.17;
min_lat = 28;
max_lat = 26;
min_lon = -85;
max_lon = -95;
len_min = 300;
% vertices of the square that any Loop current contour should pass
xv = [max_lon min_lon min_lon max_lon max_lon];
yv = [max_lat max_lat min_lat min_lat max_lat];

figure(1);clf % compare SSH
ha = tight_subplot(2,2,[.1 .03],[.07 .05],[.05 .01]);
figure(2);clf % compare SSS
hb = tight_subplot(2,2,[.1 .03],[.07 .05],[.05 .01]);
figure(3);clf % compare SSS gradient
hc = tight_subplot(2,2,[.1 .03],[.07 .05],[.05 .01]);
for ic = 1:Nc
    disp(['Working on ',casenames{ic},' ...'])
    if ic <= 3
        grd = grd1;
        grdr = grdr1;
    else
        grd = grd2;
        grdr = grdr2;
    end
    tic
    surf_file = [dir_surf,'surf_',cases{ic},'_part',num2str(Ipart),'.mat'];
    load(surf_file,'tnum_rec')
    IIt = find(tnum_rec >= tnum0 & tnum_rec <= tnum1);
    m_case = matfile(surf_file);
    disp(['  Reading records:'])
    datestr(tnum_rec(IIt))
    SSH = m_case.SSH_surf(:,:,:,IIt);
    SSS = m_case.S_surf(:,:,:,IIt);
    SSH_mean = squeeze(mean(mean(SSH,3),4));
    SSS_mean = squeeze(mean(mean(SSS,3),4));
    [S_gradmean, grdSg] = sdx_gradient(SSS_mean,grdr);
    tmp = SSH_mean;
    tmp = tmp - nanmean(tmp(:));
    [clev, xx, yy, len,ifclosed] = sdx_contourc(grdr.lon(:,1),grdr.lat(1,:),tmp',[SSH_crit SSH_crit]);
    in = inpolygon(xx,yy,xv,yv);
    in = max(in,[],2);
    II = find(in & len>len_min);
%     II = find(ifclosed == 0 & in);
    clev1 = clev(II,:);
    xx1 = xx(II,:);
    yy1 = yy(II,:);
    len1 = len(II);
    nlev = size(clev1,1);
    yy_max = max(yy1,[],2);
    [M,II] = max(yy_max);
    if ~isempty(II)
        yy_choose = yy1(II,:);
        xx_choose = xx1(II,:);
    end
    
    axes(ha(ic))
    m_proj('miller','lon',[-98 -82],'lat',[24 31]);
    m_pcolor(grdr.lon,grdr.lat,SSH_mean);shading flat;%colorbar vert;
    for ii = 1:size(xx,1)
        xx_tmp = xx(ii,:);
        yy_tmp = yy(ii,:);
        xx_tmp(isnan(xx_tmp)) = [];
        yy_tmp(isnan(yy_tmp)) = [];
        m_line(xx_tmp,yy_tmp,'linewi',1,'color','k');
    end
    m_line(xx_choose,yy_choose,'linewi',3,'color','k');
    m_gshhs_h('patch',[.8 .8 .8]);
    caxis(SSHlim)
    colormap(gca,gmap)
    m_grid('box','fancy','tickdir','out','fontsize',FS,'ytick',[24:2:30]);
    set(gca,'FontSize',FS)
    title([casenames{ic}])
    
    axes(hb(ic))
    m_proj('miller','lon',[-98 -82],'lat',[24 31]);
    m_pcolor(grdr.lon,grdr.lat,SSS_mean);shading flat;%colorbar vert;
    for ii = 1:size(xx,1)
        xx_tmp = xx(ii,:);
        yy_tmp = yy(ii,:);
        xx_tmp(isnan(xx_tmp)) = [];
        yy_tmp(isnan(yy_tmp)) = [];
        m_line(xx_tmp,yy_tmp,'linewi',1,'color','k');
    end
    m_line(xx_choose,yy_choose,'linewi',3,'color','k');
    m_gshhs_h('patch',[.8 .8 .8]);
    caxis(Slim)
    colormap(gca,cmap_Sg)
    m_grid('box','fancy','tickdir','out','fontsize',FS,'ytick',[24:2:30]);
    set(gca,'FontSize',FS)
    title([casenames{ic}])
    
    axes(hc(ic))
    m_proj('miller','lon',[-98 -82],'lat',[24 31]);
    m_pcolor(grdSg.lon,grdSg.lat,S_gradmean);shading flat;%colorbar vert;
    for ii = 1:size(xx,1)
        xx_tmp = xx(ii,:);
        yy_tmp = yy(ii,:);
        xx_tmp(isnan(xx_tmp)) = [];
        yy_tmp(isnan(yy_tmp)) = [];
        m_line(xx_tmp,yy_tmp,'linewi',1,'color','k');
    end
    m_line(xx_choose,yy_choose,'linewi',3,'color','k');
    m_gshhs_h('patch',[.8 .8 .8]);
    caxis(Sglim)
    colormap(gca,cmap_Sg)
    m_grid('box','fancy','tickdir','out','fontsize',FS,'ytick',[24:2:30]);
    set(gca,'FontSize',FS)
    title([casenames{ic}])
    toc
end