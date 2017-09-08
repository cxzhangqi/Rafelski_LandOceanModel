%% Calculated airborne fraction by fitting fossil fuel emissions to 
%% atmospheric CO2 from 1959-1979

clear all

%% Define parameters

ts = 12;
start_year = 1800;
end_year = 2005.5;


%% Load the variables

[dtdelpCO2a,dpCO2a,year,dt,CO2a] = MLOinterpolate_increment2(ts,start_year,end_year); % CO2 increment in ppm/year
% dtdelpCO2a = annual increment in atmospheric CO2
% dpCO2a = atmospheric CO2 minus initial value
% year = timeseries of dates
% dt = 1/ts
% CO2a = atmospheric CO2 interpolated to monthly values

[landusemo,ff1,fas,Aoc,extratrop_landmo] = getsourcesink_scale3;

%load land_temp.mat

%% Find FF and CO2 data from 1959 to 1979
i3 = find(ff1(:,1) == 1959);
j3 = find(ff1(:,1) == 1979);

i4 = find(floor(100*dtdelpCO2a(:,1)) == floor(100*(1959+(1/24))));
j4 = find(floor(100*dtdelpCO2a(:,1)) == floor(100*(1979+(1/24))));
% 
x3 = ff1(i3:j3,2);

d = dtdelpCO2a(i4:j4,2); %CO2 record from 1959-1979: annual increment

G2 = [ones(size(x3)) x3]; % fossil fuels from 1959-1979

m3 = G2\d; % regression of ff from 1959-1979. Second value in vector m3 is the airborne fraction 
x4 = ff1(:,2);

G3 = [ones(size(x4)) x4];

%% Scale fossil fuel record based on airborne fraction

dm3(:,1) = ff1(i3:j3,1);
dm3(:,2) = G2*m3; % fossil fuels only

dm4(:,1) = ff1(:,1);
dm4(:,2) = G3*m3; % whole fossil fuel record used

%% Integrate scale fossil fuel record and atmospheric co2 increments

[ffi] = integrate_series_trap2(dm4,1,2,12); % constant fraction of fossil fuel record

[co2i] = integrate_series_trap2(dtdelpCO2a,1,2,12); % integrate CO2 record

%% Shift constant fraction of fossil fuel record so it's on the same scale
%% as CO2 (i.e. add integration constant)
a = find(floor(100*CO2a(:,1)) == floor(100*(1958.5+(1/24))));
b = CO2a(a,2);
cc = find(ffi(:,1) == 1958.5);
ffi(:,2) = ffi(:,2) - ffi(cc,2) + b; 

e = find(floor(100*co2i(:,1)) == floor(100*(1958.5+(1/24))));
co2i(:,2) = co2i(:,2) - co2i(e,2) + b; % add integration constant to integrated CO2 record


%% Save the results of nonlin_land_Qs_annotate as a .mat file, and load
%% here:
load co2_update_VHMV_11jan11.mat


%% Fit to find integration constant to add to model results of atmospheric CO2
%% to get agreement with data

i = find(datm(:,2) == 0);
datm(i,2) = NaN;

i12 = find(datm(:,1) == 1959);
j12 = find(datm(:,1) == 1979);

X = datm(i12:j12,2);

i13 = find(floor(100*co2i(:,1)) == floor(100*(1959+(1/24))));
j13 = find(floor(100*co2i(:,1)) == floor(100*(1979+(1/24))));

y6 = co2i(i13:j13,2);

beta = 288;
[betahat,resid,J] = nlinfit(X,y6,'forward_run_fit',beta);

datm(:,2) = datm(:,2) + betahat;


%% Calculate airborne fraction anomaly
airfranom(:,1) = datm(:,1);
airfranom(:,2) = datm(:,2) - ffi(1189:3104,2);  % use for 2010 cases

j14 = find(floor(100*co2i(:,1)) == floor(100*(1850+(1/24))));
airfranom_obs(:,1) = co2i(j14:4436,1);
airfranom_obs(:,2) = co2i(j14:4436,2) - ffi(1189:3104,2);

j11 = find(floor(100*co2i(:,1)) == floor(100*(1800+(1/24))));

figure
subplot('Position',[0.075 0.5 0.9 0.4])%(2,1,1)
plot(co2i(j11:end,1),co2i(j11:end,2),ffi(:,1),ffi(:,2),datm(:,1),datm(:,2))
title('forward run, CLM-CN')
legend('measured','fossil fuel','model','Location','Southeast')
xlabel('Year')
ylabel('ppm C')
set(gca,'Xlim',[1850 2010])
set(gca,'Ylim',[280 395])
set(gca,'Xminortick','on')

subplot('Position',[0.075 0.1 0.9 0.4])%(2,1,2,'align')
plot(airfranom(:,1),airfranom(:,2),'-r',airfranom_obs(:,1),airfranom_obs(:,2),'-b')
%legend('modeled','measured','Location','Northwest')
xlabel('Year')
ylabel('ppm C')
set(gca,'Xlim',[1850 2010])
set(gca,'Xminortick','on')

