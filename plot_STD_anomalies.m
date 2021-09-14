% calculate STD of anomalies (climatological monthly average removed)
clear
Re = 6370000;
FS = 24;

file_hycom1 = '_HYCOM_20131201_20140401.nc';
file_hycom2 = '_HYCOM_20140402_20160903.nc';

dir_surf = 'D:\work\processed_surf\';
case_check = '3p5km_bulk';
case_name = 'MR passive river';

mems = {'','m1','p1'};
Nm = 1;
Npar = 3;
Nrec_par = 72;

tnum0 = datenum('2013-12-02');
tnum1 = datenum('2016-08-10');

tnum_curl = datenum('2015-6-11');

load grd_GOM3p5km.mat
fp = rnt_2grid(grd.f,'r','p');
ff = mean(fp(:));
grdr.lat = grd.lat_rho;
grdr.lon = grd.lon_rho;
grdr.mask = grd.mask_rho;
grdp.lat = grd.lat_psi;
grdp.lon = grd.lon_psi;
grdp.mask = grd.mask_psi;

tnum_orig = datenum(2013,1,1);

% load hycom fields
tnum_orig2 = datenum(1900,12,31);
tnum_hyc1 = ncread(['STUV',file_hycom1],'MT')+tnum_orig2;
tnum_hyc2 = ncread(['STUV',file_hycom2],'MT')+tnum_orig2;
latr = ncread(['STUV',file_hycom1],'Latitude');
lonr = ncread(['STUV',file_hycom1],'Longitude');
[latr,lonr]=meshgrid(latr,lonr);
grdr_hyc.lat = double(latr);
grdr_hyc.lon = double(lonr);
grdr_hyc.mask = ones(size(latr));

load 150_sunshine_9lev.mat
% cmap=...
%  [...
% 255 255 255
% 255 232 120
% 255 192  60
% 255 160   0
% 255  96   0
% 255  50   0
% 225  20   0
% 192   0   0
% 165   0   0 ...
% ]./255;

Tlim = [0 2];%[22 32];
Slim = [0 2];
%% calculate climatology of hycom data
disp('Working on Hycom data ...')
T_hyc1 = squeeze(ncread(['STUV',file_hycom1],'temperature',[1 1 1 1],[Inf Inf 1 Inf]));
S_hyc1 = squeeze(ncread(['STUV',file_hycom1],'salinity',[1 1 1 1],[Inf Inf 1 Inf]));
T_hyc2 = squeeze(ncread(['STUV',file_hycom2],'temperature',[1 1 1 1],[Inf Inf 1 Inf]));
S_hyc2 = squeeze(ncread(['STUV',file_hycom2],'salinity',[1 1 1 1],[Inf Inf 1 Inf]));
tnum_hyc = [tnum_hyc1;tnum_hyc2];
T_hyc = cat(3,T_hyc1,T_hyc2);
S_hyc = cat(3,S_hyc1,S_hyc2);
IIt = find(tnum_hyc >= tnum0 & tnum_hyc<=tnum1);
tnum_hyc = tnum_hyc(IIt);
T_hyc = T_hyc(:,:,IIt);
S_hyc = S_hyc(:,:,IIt);
[Y,M,D] = datevec(tnum_hyc);
T2 = zeros(size(grdr_hyc.mask));
S2 = zeros(size(grdr_hyc.mask));
for imon = 1:12
    disp(['  Calculating anomalies on month ',num2str(imon),' ...'])
    II = find(M==imon);
    Tmean = mean(T_hyc(:,:,II),3);
    Tano = T_hyc(:,:,II) - Tmean;
    T2 = T2+sum(Tano.^2,3);
    Smean = mean(S_hyc(:,:,II),3);
    Sano = S_hyc(:,:,II) - Smean;
    S2 = S2+sum(Sano.^2,3);
end
Tstd_hyc = sqrt(sum(T2,3)/(length(IIt)-1));
Sstd_hyc = sqrt(sum(S2,3)/(length(IIt)-1));

% go through surface files
disp('Working on model data ...')
disp('  Reading date number ...')
tnum_all = [];
for ip = 1:Npar
    surf_file = [dir_surf,'surf_',case_check,'_part',num2str(ip),'.mat'];
    load(surf_file,'tnum_rec');
    tnum_all = [tnum_all;tnum_rec];
end
IIt = find(tnum_all >= tnum0 & tnum_all<=tnum1);
[Y,M,D] = datevec(tnum_all(IIt));
T2 = zeros(size(grdr.mask));
S2 = zeros(size(grdr.mask));
count = 0;
for imon = 1:12
    disp(['  Calculating anomalies on month ',num2str(imon),' ...'])
    II = find(M==imon);
    IIrec = IIt(II);
    Ttmp = [];
    Stmp = [];
    for ip = 1:Npar
        IItmp = IIrec - (ip-1)*Nrec_par;
        IItmp2 = find(IItmp>0 & IItmp<=Nrec_par);
        IIchoose = IItmp(IItmp2);
        if ~isempty(IIchoose)
            surf_file = [dir_surf,'surf_',case_check,'_part',num2str(ip),'.mat'];
            disp(['  Reading data from ',surf_file,' ...'])
            disp(['    ',num2str(length(IIchoose)),' record will be read'])
            m_check = matfile(surf_file);
            Ttmp = cat(3,Ttmp,squeeze(m_check.T_surf(:,:,1,IIchoose)));
            Stmp = cat(3,Stmp,squeeze(m_check.S_surf(:,:,1,IIchoose)));
            count = count+length(IIchoose);
        end
    end
    Tmean = mean(Ttmp,3);
    Tano = Ttmp - Tmean;
    T2 = T2+sum(Tano.^2,3);
    Smean = mean(Stmp,3);
    Sano = Stmp - Smean;
    S2 = S2+sum(Sano.^2,3);
end
if count ~= length(IIt)
    disp('ERROR: record number does not match!!!')
    return
end
Tstd = sqrt(sum(T2,3)/(length(IIt)-1));
Sstd = sqrt(sum(S2,3)/(length(IIt)-1));

figure(1);clf
ha = tight_subplot(2,2,[.1 .03],[.07 .05],[.05 .01]);
axes(ha(1))
m_proj('miller','lon',[-98 -82],'lat',[24 31]);
m_pcolor(grdr_hyc.lon,grdr_hyc.lat,Tstd_hyc);shading flat;%colorbar vert;
m_gshhs_h('patch',[.8 .8 .8]);
caxis(Tlim)
colormap(gca,cmap)
m_grid('box','fancy','tickdir','out','fontsize',FS,'ytick',[24:2:30]);
set(gca,'FontSize',FS)
title(['Hycom STD of SST anomalies'])
axes(ha(2))
m_proj('miller','lon',[-98 -82],'lat',[24 31]);
m_pcolor(grdr_hyc.lon,grdr_hyc.lat,Sstd_hyc);shading flat;%colorbar vert;
m_gshhs_h('patch',[.8 .8 .8]);
caxis(Slim)
colormap(gca,cmap)
m_grid('box','fancy','tickdir','out','fontsize',FS,'ytick',[24:2:30]);
set(gca,'FontSize',FS)
title(['Hycom STD of SSS anomalies'])

axes(ha(3))
m_proj('miller','lon',[-98 -82],'lat',[24 31]);
m_pcolor(grdr.lon,grdr.lat,Tstd);shading flat;%colorbar vert;
m_gshhs_h('patch',[.8 .8 .8]);
caxis(Tlim)
colormap(gca,cmap)
m_grid('box','fancy','tickdir','out','fontsize',FS,'ytick',[24:2:30]);
set(gca,'FontSize',FS)
title([case_name,' STD of SST anomalies'])
axes(ha(4))
m_proj('miller','lon',[-98 -82],'lat',[24 31]);
m_pcolor(grdr.lon,grdr.lat,Sstd);shading flat;
% colorbar('southoutside','Ticks',[0:0.2:2])
m_gshhs_h('patch',[.8 .8 .8]);
caxis(Slim)
colormap(gca,cmap)
m_grid('box','fancy','tickdir','out','fontsize',FS,'ytick',[24:2:30]);
set(gca,'FontSize',FS)
title([case_name,' STD of SSS anomalies'])
