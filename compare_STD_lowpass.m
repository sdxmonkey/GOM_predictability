% compare the low-resolution data and low-passed high resolution data
clear
FS = 20;
LW = 3;
dir_surf = 'D:\work\processed_surf\';
case_checks = {'3p5km_bulk_priv','1km_bulk_priv'};
Ncase = length(case_checks);
Nm = 2;
ip_max = 3;

xrange = [-100 -80];
yrange = [27 35];

f_cut = 1/3.5/2;
N_butt = 4;

tnum0 = datenum('2013-12-02');
tnum1 = datenum('2016-08-10');
curlrange1 = [0 1];
Cgrey = [.9,.9,.9];

tnum_break = [datenum('2014-04-02'),datenum('2014-08-02'),datenum('2014-12-02'),...
    datenum('2015-04-02'),datenum('2015-08-02'),datenum('2015-12-02'),...
    datenum('2016-04-02')];

load grd_GOM1km.mat
fp = rnt_2grid(grd.f,'r','p');
ff = mean(fp(:));
grdp.lat = grd.lat_psi;
grdp.lon = grd.lon_psi;
grdp.mask = grd.mask_psi;
grdr.lat = grd.lat_rho;
grdr.lon = grd.lon_rho;
grdr.mask = grd.mask_rho;
IIr_high = find(grdr.lat >= yrange(1) & grdr.lat <= yrange(2) & grdr.lon >= xrange(1) & grdr.lon <= xrange(2) & grdr.mask == 1);
IIp_high = find(grdp.lat >= yrange(1) & grdp.lat <= yrange(2) & grdp.lon >= xrange(1) & grdp.lon <= xrange(2) & grdp.mask == 1);
grd_high = grd;
grdp_high = grdp;
grdr_high = grdr;
dx_high = 1;

load grd_GOM3p5km.mat
grdp.lat = grd.lat_psi;
grdp.lon = grd.lon_psi;
grdp.mask = grd.mask_psi;
grdr.lat = grd.lat_rho;
grdr.lon = grd.lon_rho;
grdr.mask = grd.mask_rho;
IIr_low = find(grdr.lat >= yrange(1) & grdr.lat <= yrange(2) & grdr.lon >= xrange(1) & grdr.lon <= xrange(2) & grdr.mask == 1);
IIp_low = find(grdp.lat >= yrange(1) & grdp.lat <= yrange(2) & grdp.lon >= xrange(1) & grdp.lon <= xrange(2) & grdp.mask == 1);
grd_low = grd;
grdp_low = grdp;
grdr_low = grdr;
dx_low = 3.5;

%build time ticks
tk_yr = 2013:2016;
tk_yr = repmat(tk_yr,[12,1]);
tk_yr = tk_yr(:);
tk_mon = 1:12;
tk_mon = repmat(tk_mon,[1,4]);
tk_num = datenum(tk_yr,tk_mon',1);

figure(1);clf % check the overall change of curl std

for ic = 1:Ncase
    if ic == 1
        case_end = '_low';
    else
        case_end = '_high';
    end
    eval(['grd = grd',case_end,';'])
    eval(['grdp = grdp',case_end,';'])
    eval(['grdr = grdr',case_end,';'])
    eval(['dx = dx',case_end,';'])
    eval(['IIp = IIp',case_end,';'])
    eval(['IIr = IIr',case_end,';'])
    ip = 1;
    it = 0;
    tnum = [];
    std_curl = [];
    while ip <= ip_max
        surf_file = [dir_surf,'surf_',case_checks{ic},'_part',num2str(ip),'.mat'];
        disp(['Reading from ',surf_file,' ...'])
        m_case = matfile(surf_file);
        if it == 0
            Irec = find(m_case.tnum_rec >= tnum0 & m_case.tnum_rec<= tnum1);
        else
            Irec = find(m_case.tnum_rec > tnum(it) & m_case.tnum_rec<= tnum1);
        end
        if ~isempty(Irec)
            tnum_rec = m_case.tnum_rec(Irec,1);
            U = m_case.U_surf(:,:,1:Nm,Irec);
            V = m_case.V_surf(:,:,1:Nm,Irec);
            for irec = 1:length(Irec)
                disp(['  Calulating on ',datestr(tnum_rec(irec)),' ...'])
                if it ~= 0 && ismember(tnum_rec(irec),tnum_break)
                    disp('    Breaking point!!!')
                    it = it + 1;
                    tnum(it) = NaN;
                    std_curl(it) = NaN;
                    if ic == 2
                        std_curl_intp(it) = NaN;
                        std_curl_fil(it) = NaN;
                        std_curl_filintp(it) = NaN;
                    end
                end
                it = it + 1;
                tnum(it) = tnum_rec(irec);
                U_rec = U(:,:,:,irec);
                V_rec = V(:,:,:,irec);
                curlf = rnt_curl(U_rec,V_rec,grd)./ff;
                curlstd = std(curlf,0,3);
                std_curl(it) = sqrt(sum(curlstd(IIp).^2)/length(IIp));
                if ic == 2
                    disp(['    Interpolate and filter the data ...'])
                    tic
                    IIUh = find(grd_high.mask_u == 1);
                    IIUl = find(grd_low.mask_u == 1);
                    IIVh = find(grd_high.mask_v == 1);
                    IIVl = find(grd_low.mask_v == 1);
                    curlf_intp = nan([size(grdp_low.mask),Nm]);
                    curlf_fil = nan([size(grdp_high.mask),Nm]);
                    curlf_filintp = nan([size(grdp_low.mask),Nm]);
                    f_butt_u = sdx_lowpassfilter(size(grd_high.mask_u), dx_high, f_cut, N_butt);
                    f_butt_v = sdx_lowpassfilter(size(grd_high.mask_v), dx_high, f_cut, N_butt);
                    for im = 1:Nm
                        U_tmp = U_rec(:,:,im);
                        V_tmp = V_rec(:,:,im);
                        
                        F_U = scatteredInterpolant(grd_high.lon_u(IIUh),grd_high.lat_u(IIUh),U_tmp(IIUh));
                        U_intp = nan(size(grd_low.mask_u));
                        U_intp(IIUl) = F_U(grd_low.lon_u(IIUl),grd_low.lat_u(IIUl));
                        F_V = scatteredInterpolant(grd_high.lon_v(IIVh),grd_high.lat_v(IIVh),V_tmp(IIVh));
                        V_intp = nan(size(grd_low.mask_v));
                        V_intp(IIVl) = F_V(grd_low.lon_u(IIVl),grd_low.lat_u(IIVl));
                        curlf_intp(:,:,im) = rnt_curl(U_intp,V_intp,grd_low)./ff;
                        
                        tmp = U_tmp;
                        tmp(grd_high.mask_u==0) = 0;
                        Y = fftshift(fft2(tmp));
                        Yhat = f_butt_u.*Y;
                        U_filtered = abs(ifft2(ifftshift(Yhat)));
                        U_filtered(grd_high.mask_u==0) = NaN;
                        tmp = V_tmp;
                        tmp(grd_high.mask_v==0) = 0;
                        Y = fftshift(fft2(tmp));
                        Yhat = f_butt_v.*Y;
                        V_filtered = abs(ifft2(ifftshift(Yhat)));
                        V_filtered(grd_high.mask_v==0) = NaN;
                        curlf_fil(:,:,im) = rnt_curl(U_filtered,V_filtered,grd_high)./ff;
                        
                        F_U = scatteredInterpolant(grd_high.lon_u(IIUh),grd_high.lat_u(IIUh),U_filtered(IIUh));
                        U_filintp = nan(size(grd_low.mask_u));
                        U_filintp(IIUl) = F_U(grd_low.lon_u(IIUl),grd_low.lat_u(IIUl));
                        F_V = scatteredInterpolant(grd_high.lon_v(IIVh),grd_high.lat_v(IIVh),V_filtered(IIVh));
                        V_filintp = nan(size(grd_low.mask_v));
                        V_filintp(IIVl) = F_V(grd_low.lon_u(IIVl),grd_low.lat_u(IIVl));
                        curlf_filintp(:,:,im) = rnt_curl(U_filintp,V_filintp,grd_low)./ff;
                    end
                    toc
                    curlstd = std(curlf_intp,0,3);
                    std_curl_intp(it) = sqrt(sum(curlstd(IIp_low).^2)/length(IIp_low));
                    curlstd = std(curlf_fil,0,3);
                    std_curl_fil(it) = sqrt(sum(curlstd(IIp_high).^2)/length(IIp_high));
                    curlstd = std(curlf_filintp,0,3);
                    std_curl_filintp(it) = sqrt(sum(curlstd(IIp_low).^2)/length(IIp_low));
                end
            end
        end
        ip = ip+1;
    end
    eval(['tnum',case_end,' = tnum;'])
    eval(['std_curl',case_end,' = std_curl;'])
end

figure(1)
plot(tnum_low,std_curl_low,'k--','LineWidth',LW)
hold on
plot(tnum_high,std_curl_high,'k','LineWidth',LW)
% plot(tnum_high,std_curl_intp,'g','LineWidth',2)
% plot(tnum_high,std_curl_fil,'r','LineWidth',LW)
plot(tnum_high,std_curl_filintp,'LineWidth',LW,'Color',[.7 .7 .7])
for ip = 1:length(tnum_break)
    plot([tnum_break(ip),tnum_break(ip)],[curlrange1(1) curlrange1(2)],'LineWidth',LW,'Color',Cgrey)
end
hold off
ylabel('STD of curl/f')
% xlim([tnum0 tnum1])
xlim([datenum('2013-12-01') datenum('2016-08-31')])
ylim(curlrange1)
xticks(tk_num(4:3:end))
datetick('x','mm/yy','keeplimits','keepticks')
% legend('3.5 km','1 km','1 l km interpolated','1 km low-pass','1 km low-pass & interpolated','Location','best')
legend('MR active river','SP active river','SP low-pass & interpolated','Location','best')
set(gca,'FontSize',FS)

II = find(~isnan(std_curl_low) & ~isnan(std_curl_high));
[rr,pp] = corrcoef(std_curl_low(II),std_curl_high(II))