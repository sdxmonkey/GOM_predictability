% plot STD timeseries using data
clear
LW = 3;
FS = 20;
file_in = 'STD_LC_ensembles.mat';
load(file_in)
cases = {'1km_bulk_priv','1km_bulk','1km_bulk_nor','3p5km_bulk_priv','3p5km_bulk'};
Nc = length(cases);
casenames = {'SP active river','SP passive river','SP no river','MR active river','MR passive river'};
lineC = {'k','r','b','k--','r--','g'};
Cgrey = [.9,.9,.9];
Np = 8; % 8 periods in total

SSTrange1 = [0 1];
SSSrange1 = [0 2.5];
SSHrange1 = [0 0.1];
curlrange1 = [0 1];
ymax_shading = 3;
ymin_shading = 0;

%build time ticks
tk_yr = 2013:2016;
tk_yr = repmat(tk_yr,[12,1]);
tk_yr = tk_yr(:);
tk_mon = 1:12;
tk_mon = repmat(tk_mon,[1,4]);
tk_num = datenum(tk_yr,tk_mon',1);

figure(1);clf % check the overall change of std
hold on
figure(2);clf % check the overall change of std
hold on
figure(3);clf % check the overall change of std
hold on
figure(4);clf % check the overall change of std
hold on

% find the overlap regions and plot them as grey shading
eval(['tnum_tmp = tnum_',cases{2},';'])
for ip = 2:Np
    xmin_shading = min(tnum_tmp(ip,:));
    xmax_shading = max(tnum_tmp(ip-1,:));
    for iplot = 1:4
        figure(iplot)
        patch([xmin_shading xmax_shading xmax_shading xmin_shading],[ymin_shading ymin_shading ymax_shading ymax_shading],Cgrey,'EdgeColor',Cgrey)
    end
end

figure(1)
for ic = 1:Nc
    eval(['tnum_E = tnum_',cases{ic},';'])
    eval(['STD_tmp = STD_T_',cases{ic},';'])
    tnum_E = tnum_E';
    STD_tmp = STD_tmp';
    ha(ic) = plot(tnum_E(:),STD_tmp(:),lineC{ic},'LineWidth',LW);
end

figure(2)
for ic = 1:Nc
    eval(['tnum_E = tnum_',cases{ic},';'])
    eval(['STD_tmp = STD_S_',cases{ic},';'])
    tnum_E = tnum_E';
    STD_tmp = STD_tmp';
    plot(tnum_E(:),STD_tmp(:),lineC{ic},'LineWidth',LW);
end

figure(3)
for ic = 1:Nc
    eval(['tnum_E = tnum_',cases{ic},';'])
    eval(['STD_tmp = STD_SSH_',cases{ic},';'])
    tnum_E = tnum_E';
    STD_tmp = STD_tmp';
    plot(tnum_E(:),STD_tmp(:),lineC{ic},'LineWidth',LW);
end

figure(4)
for ic = 1:Nc
    eval(['tnum_E = tnum_',cases{ic},';'])
    eval(['STD_tmp = STD_curl_',cases{ic},';'])
    tnum_E = tnum_E';
    STD_tmp = STD_tmp';
    hb(ic) = plot(tnum_E(:),STD_tmp(:),lineC{ic},'LineWidth',LW);
end
figure(1)
hold off
ylabel('STD of SST')
ylim(SSTrange1)
xlim([datenum('2013-12-01') datenum('2016-08-31')])
% xticks([datenum('2013-12-1'):90:datenum('2016-12-31')])
xticks(tk_num(4:3:end))
datetick('x','mm/yy','keeplimits','keepticks')
legend(ha,[casenames])
% legend([casenames,'1km MEM','3.5 km MEM'])
set(gca,'FontSize',FS)

figure(2)
hold off
ylabel('STD of SSS')
ylim(SSSrange1)
xlim([datenum('2013-12-01') datenum('2016-08-31')])
% xticks([datenum('2013-12-1'):90:datenum('2016-12-31')])
xticks(tk_num(4:3:end))
datetick('x','mm/yy','keeplimits','keepticks')
% legend([casenames])
set(gca,'FontSize',FS)

figure(3)
hold off
ylabel('STD of SSH')
ylim(SSHrange1)
xlim([datenum('2013-12-01') datenum('2016-08-31')])
% xticks([datenum('2013-12-1'):90:datenum('2016-12-31')])
xticks(tk_num(4:3:end))
datetick('x','mm/yy','keeplimits','keepticks')
% legend([casenames])
set(gca,'FontSize',FS)

figure(4)
hold off
ylabel('STD of \zeta/f')
ylim(curlrange1)
xlim([datenum('2013-12-01') datenum('2016-08-31')])
% xticks([datenum('2013-12-1'):90:datenum('2016-12-31')])
xticks(tk_num(4:3:end))
datetick('x','mm/yy','keeplimits','keepticks')
legend(hb,[casenames])
set(gca,'FontSize',FS)