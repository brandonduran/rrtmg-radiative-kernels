% calculate SW radiative kernel for ice clouds as a function of ice water
% path and ice-crystal effective radius
clearvars; clc;
%addpath('/home/cawall/RRTMG');
%addpath('/home/cawall/scripts/scripts_for_Brandon/functions/');

IWP_bins_MODIS=[0, 20, 50, 100, 200, 400, 1000, 20000]; % ice water path (g m-2)
Reff_bins_MODIS=[5, 10, 20, 30, 40, 50, 60, 90]; % ice-crystal effective radius (um)
IWP_mdpts_MODIS=(IWP_bins_MODIS(1:(end-1))+IWP_bins_MODIS(2:end))./2;
Reff_mdpts_MODIS=(Reff_bins_MODIS(1:(end-1))+Reff_bins_MODIS(2:end))./2;

% load seasonal climatology of temperature and humidity from GCM simulations for kernel calcualtion
d=load('/home/cawall/scripts/scripts_for_Brandon/data/CTL_seasonal_climatology_zonal_mean.mat',...
    'lat','months','p','q','t','ts');
lat=d.lat;
months=d.months;
p=d.p;
q=mean(d.q,4); % multi-model mean
t=mean(d.t,4); % multi-model mean
ts=mean(d.ts,3); % multi-model mean
clearvars d

% surface albedo values from Zelinka et al. 2012
albcs=[0 0.5 1];

% save data for loops
t_save=t; q_save=q; ts_save=ts; lat_save=lat; months_save=months;

%%%%%%%%%%%%% RRTMG parameters that don't change %%%%%%%%%
% transpose to match dimensions of RRTMG
p=p';
% add one level for the surface
p(end+1)=1012.5;
% switch from pressure-increasing indexing to pressure-decreasing indexing
p=fliplr(p);
play=(p(1:(end-1))+p(2:end))./2; % layer midpoint pressure (hPa)
plev=p; % interface pressure (hPa)
nlay=numel(play);

% get GHG soundings
[plev_std,play_std,tlev_std,tlay_std,~,~,gasvmr_std,~,~] = read_atm(2); % 1: standard tropical atmosphere, 2: midlatitude summer, 3: midlatitude winter
co2vmr = 330e-6*ones(size(play)); % value from Zelinka et al. 2012 (https://doi.org/10.1175/JCLI-D-11-00248.1)
ch4vmr = 1.6e-6*ones(size(play)); % value from Zelinka et al. 2012
n2ovmr = 0.28e-6*ones(size(play)); % value from Zelinka et al. 2012
o3vmr = interp1(play_std,gasvmr_std(3,:),play,'linear','extrap');
covmr = interp1(play_std,gasvmr_std(5,:),play,'linear','extrap');
o2vmr = interp1(play_std,gasvmr_std(7,:),play,'linear','extrap');
clearvars gasvmr_std plev_std play_std tlev_std tlay_std t q ts mwch4 mwco2 mwdry mwf11 mwf12 mwh2o mwn2o mwo3

% molecular weights for gases
mwdry =  28.966;
mwh2o =  18.016;

ncol=1;
icld=1; % cloud overlap assumption. 1: random overlap

adjes = 1; % if set to 1 then a constant Earth-Sun distance of 1 AU is used
dyofyr = 1; % dummy variable if adjes==1
scon=1366; % value from Zelinka et al. 2012 (W m-2)

inflgsw = 2; % value of 1 or 2 means that cloud physical properties are intput to the model. Value of 0 means that cloud optical depth, single-scattering albedo, etc. are input to the model
iceflgsw = 3; % 3 means generalized effective particle size (Dge) from Fu 1996 J. Cli.(https://doi.org/10.1175/1520-0442(1996)009<2058:AAPOTS>2.0.CO;2)
liqflgsw = 1; % 1 means liquid clouds are active in the model, 0 means inactive
iaer = 0; % 0 means aerosols are not active in the calculation

% these are dummy variables. they must be input to RRTMG but they aren't
% used in the calculation
get_rrtmg_lw_parms
get_rrtmg_sw_parms
taucld = zeros(nbndsw,ncol,nlay);
ssacld = zeros(nbndsw,ncol,nlay);
fsfcld = zeros(nbndsw,ncol,nlay);
asmcld = zeros(nbndsw,ncol,nlay);
tauaer = zeros(ncol,nlay,nbndsw);
ssaaer = zeros(ncol,nlay,nbndsw);
asmaer = zeros(ncol,nlay,nbndsw);
ecaer = zeros(ncol,nlay,naerec);

CTP=250; % assume CTP of 250 hPa for ice clouds
% find level that equal CTP
if CTP>=max(plev)
    CTP_ind=find(plev==max(plev));
else
    CTP_ind=find(plev==CTP)-1;
end
%%%%%%%%%%% end of RRTMG parameters that don't change %%%%%%%

% initialize array for kernel
SW_kernel=NaN*ones(numel(albcs),numel(lat_save),numel(Reff_mdpts_MODIS),numel(IWP_mdpts_MODIS),numel(months_save));

for albedo_ind=1:numel(albcs)
for lat_ind=1:numel(lat_save)
for month_ind=1:numel(months_save)

t=squeeze(t_save(lat_ind,:,month_ind));
q=squeeze(q_save(lat_ind,:,month_ind));
ts=squeeze(ts_save(lat_ind,month_ind));
lat=lat_save(lat_ind);
months=months_save(month_ind);
asdir=albcs(albedo_ind);
aldir=albcs(albedo_ind);
asdif=albcs(albedo_ind);
aldif=albcs(albedo_ind);

if months==2
    day_of_yr=datenum(2000,months,14)-datenum(2000,1,1)+1;
else
    day_of_yr=datenum(2000,months,15)-datenum(2000,1,1)+1;
end
hrs=0:24;
for n=1:(numel(hrs)-1)
    dt=0.01;
    day_num=day_of_yr+((hrs(n):dt:(hrs(n+1)-dt)))./24;
    cosSZA = compute_cosSZA(lat,0,day_num);
    cosSZA(cosSZA<0)=0;
    cosSZA_mean(n)=nanmean(cosSZA);
end
clearvars day_of_yr hrs n dt day_num cosSZA

% add one level for the surface
t(end+1)=ts; q(end+1)=q(end);

% switch from pressure-increasing indexing to pressure-decreasing indexing
q=fliplr(q); t=fliplr(t);

% fill missing values where pressure levels are always below the surface (only affects Antarctica)
if any(isnan(t))
    t=fill_t_q_profile(t,p);
end

if any(isnan(q))
    q=fill_t_q_profile(q,p);
end

% RRTMG settings
tlay=(t(1:(end-1))+t(2:end))./2; % layer midpoint temperature (K)
tlev=t; % interface temperature (K)
tsfc=ts; %surface temperature (K)
h2ovmr=(mwdry/mwh2o).*(q(1:(end-1))+q(2:end))./2; % water vapor mole fraction

% calculate SW CRE
SWCRE=NaN*ones(numel(Reff_bins_MODIS),numel(IWP_bins_MODIS));

for n=1:numel(Reff_bins_MODIS)
    for m=1:numel(IWP_bins_MODIS)
        Reff=Reff_bins_MODIS(n); % um
        IWP=IWP_bins_MODIS(m); % g m-2
        
        Dge=8./(3*sqrt(3)).*Reff; % equation 3.12 of Fu 1996

        cicewp=zeros(size(play)); cliqwp=zeros(size(play)); reice=zeros(size(play)); reliq=zeros(size(play)); cldfr=zeros(size(play));
        cicewp(CTP_ind)=IWP; 
        reice(CTP_ind)=Dge; 
        clearvars Dge Reff

        cldfr(cicewp>0)=1;
        
        % radiative transfer calculations. loop over SZA for diurnal cycle
        for k=1:numel(cosSZA_mean)
            coszen=cosSZA_mean(k);
            % CALL RRTMG SW MEX INTERFACE
            [uflx_sw, ~, ~, uflxc_sw, ~, ~]  ...
              = rrtmg_sw_wrapper(ncol,nlay,icld, ...
                         play,plev,tlay,tlev,tsfc, ...
                         h2ovmr,o3vmr,co2vmr,ch4vmr, ...
                         n2ovmr,o2vmr, ...
                         asdir, asdif, aldir, aldif, ...
                         coszen, adjes, dyofyr, scon, ...
                         inflgsw,iceflgsw,liqflgsw,iaer, ...
                         cldfr,taucld,ssacld,asmcld,fsfcld, ...
                         cicewp,cliqwp,reice,reliq, ...
                         tauaer,ssaaer,asmaer,ecaer);
             uflx_sw_TOA_save(k)=uflx_sw(end);
             uflxc_sw_TOA_save(k)=uflxc_sw(end);
        end
        uflx_sw_TOA=nanmean(uflx_sw_TOA_save);
        uflxc_sw_TOA=nanmean(uflxc_sw_TOA_save);
        clearvars k uflx_sw_TOA_save uflxc_sw_TOA_save
        
        if IWP==0
            SWCRE(n,m)=0;
        else
            SWCRE(n,m)=uflxc_sw_TOA - uflx_sw_TOA;
        end
        clearvars IWP uflx_sw_TOA uflxc_sw_TOA
    end
end
clearvars n m

for n=1:numel(Reff_mdpts_MODIS)
    for m=1:numel(IWP_mdpts_MODIS)
        SWCREa=SWCRE(n:(n+1),m:(m+1));
        SW_kernel(albedo_ind,lat_ind,n,m,month_ind)=mean(SWCREa(:))./100; % divide by 100 to convert from overcast CRE to W m-2 %-1
    end
end

if lat_ind==1 & albedo_ind==1
    month_ind
end

clearvars -except albcs t_save q_save ts_save lat_save months_save p CTP_ind CTP SW_kernel IWP_bins_MODIS IWP_mdpts_MODIS ...
    Reff_bins_MODIS Reff_mdpts_MODIS lat_ind month_ind albedo_ind p play plev nlay co2vmr ch4vmr n2ovmr o3vmr covmr o2vmr mwdry mwh2o ...
    ncol icld adjes dyofyr scon inflgsw iceflgsw liqflgsw iaer taucld ssacld fsfcld asmcld tauaer ssaaer asmaer ecaer

end % end of month loop
lat_ind
end % end of lat loop
end % end of albedo loop

%% save kernel
% user needs to input the file name to save the kernel as a string:
% file_name_save='/path/to/file/filename.mat';

Reff_bins=Reff_bins_MODIS; Reff_mdpts=Reff_mdpts_MODIS;
IWP_bins=IWP_bins_MODIS; IWP_mdpts=IWP_mdpts_MODIS;
lat=lat_save;
months=months_save;

save(file_name_save,'Reff_bins','Reff_mdpts','IWP_bins','IWP_mdpts',...
    'SW_kernel','albcs','lat','months');

%% plot to check
clearvars; clc;
load(file_name_save,'Reff_bins','Reff_mdpts','IWP_bins','IWP_mdpts',...
    'SW_kernel','albcs','lat','months');
SW_kernela=squeeze(nanmean(SW_kernel,5));

for l=1:numel(lat)
    for n=1:numel(Reff_mdpts)
        for m=1:numel(IWP_mdpts)
            SW_kernel_mean(l,n,m)=interp1(albcs,squeeze(SW_kernela(:,l,n,m)),0.07);
        end
    end
end
clearvars n m SW_kernela

for n=1:numel(Reff_mdpts)
    for m=1:numel(IWP_mdpts)
        SW_kernel_meana(n,m)=sum(squeeze(SW_kernel_mean(:,n,m)).*cos(lat.*pi./180))./sum(cos(lat.*pi./180));
    end
end
SW_kernel_mean=SW_kernel_meana;
clearvars n m SW_kernel_meana

addpath('/home/cawall/scripts/general_use/DrosteEffect-BrewerMap-7fdcfcf/')
cm=flipud(brewermap(100,'blues'));
aspect_ratio=[9 7 1];

figure
imagesc(SW_kernel_mean)
colormap(cm)
set(gca,'tickdir','out','box','on','fontsize',12,'xtick',(1:numel(IWP_bins))-0.5,'xticklabel',IWP_bins,...
    'ytick',(1:numel(Reff_bins))-0.5,'yticklabel',Reff_bins)
colorbar
caxis([-2.5 0])
pbaspect(aspect_ratio)    
