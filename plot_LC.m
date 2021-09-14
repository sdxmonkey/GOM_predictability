% plot LC extension timeseries
clear
LW = 3;
FS = 20;
file_in = 'STD_LC_ensembles.mat';
hycom_in = 'tmp_LC_north.mat';
load(file_in)
cases = {'1km_bulk_nor','1km_bulk','1km_bulk_priv','3p5km_bulk','3p5km_bulk_priv'};
Nc = length(cases);
casenames = {'SP no river','SP passive river','SP active river','MR passive river','MR active river'};
lineC = {'b','r','k','r--','k--'};
Cgrey = [.9,.9,.9];
Chyc = [.6,.6,.6];
Np = 8; % 8 periods in total

LCrange = [25 30]
ymax_shading = 24;
ymin_shading = 32;

%build time ticks
tk_yr = 2013:2016;
tk_yr = repmat(tk_yr,[12,1]);
tk_yr = tk_yr(:);
tk_mon = 1:12;
tk_mon = repmat(tk_mon,[1,4]);
tk_num = datenum(tk_yr,tk_mon',1);

figure(1);clf % 1km runs
hold on
figure(2);clf % 3.5km runs
hold on

% find the overlap regions and plot them as grey shading
eval(['tnum_tmp = tnum_',cases{2},';'])
for ip = 2:Np
    xmin_shading = min(tnum_tmp(ip,:));
    xmax_shading = max(tnum_tmp(ip-1,:));
    for iplot = 1:2
        figure(iplot)
        patch([xmin_shading xmax_shading xmax_shading xmin_shading],[ymin_shading ymin_shading ymax_shading ymax_shading],Cgrey,'EdgeColor',Cgrey)
    end
end

mean_LC_STD = nan(Nc,1);
for ic = 1:Nc
    eval(['tnum_E = tnum_',cases{ic},';'])
    eval(['LC_tmp = LC_',cases{ic},';'])
    eval(['LC_STD_tmp = LC_STD_',cases{ic},';'])
    tnum_E = tnum_E';
    LC_tmp = LC_tmp';
    LC_STD_tmp = LC_STD_tmp';
    mean_LC_STD(ic) = nanmean(LC_STD_tmp(:));
    if ic <= 3
        figure(1)
    else
        figure(2)
    end
    ha(ic) = errorbar(tnum_E(:),LC_tmp(:),LC_STD_tmp(:),lineC{ic},'LineWidth',LW);
%     ha(ic) = plot(tnum_E(:),LC_tmp(:),lineC{ic},'LineWidth',LW);
end

% plot hycom data
load(hycom_in,'LC_north_hyc','tnum_E_hyc')
figure(1)
hb(1) = plot(tnum_E_hyc,LC_north_hyc,'LineWidth',LW,'Color',Chyc);
hold off
ylim(LCrange)
ylabel('max lat of LC')
xlim([datenum('2013-12-01') datenum('2016-08-31')])
%     xticks([datenum('2013-12-1'):90:datenum('2016-12-31')])
xticks(tk_num(4:3:end))
datetick('x','mm/yy','keeplimits','keepticks')
% title('mean of each case')
legend([ha(3:-1:1),hb(1)],{casenames{3:-1:1},'hycom'})
set(gca,'FontSize',FS)

figure(2)
hb(2) = plot(tnum_E_hyc,LC_north_hyc,'LineWidth',LW,'Color',Chyc);
hold off
ylim(LCrange)
ylabel('max lat of LC')
xlim([datenum('2013-12-01') datenum('2016-08-31')])
%     xticks([datenum('2013-12-1'):90:datenum('2016-12-31')])
xticks(tk_num(4:3:end))
datetick('x','mm/yy','keeplimits','keepticks')
% title('mean of each case')
legend([ha(5:-1:4),hb(2)],{casenames{5:-1:4},'hycom'})
set(gca,'FontSize',FS)