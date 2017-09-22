% Change to fit ff from 1959-1980

clear all

% Define parameters

ts = 12;
start_year = 1800;
end_year = 2005.5;


% Load the variables

[dtdelpCO2a,dpCO2a,year,dt,CO2a] = MLOinterpolate_increment(ts,start_year,end_year); % CO2 increment in ppm/year

[landusemo,ff1,fas,Aoc,extratrop_landmo] = getsourcesink_scale2;

[tau,nino3] = get_climateindices;

load land_temp.mat
i2 = find(nino3(:,1) == 1958.5);
j2 = find(nino3(:,1) == 2005.5);


% i3 = find(ff1(:,1) == 1958.5);
% j3 = find(ff1(:,1) == 2005.5);

i3 = find(ff1(:,1) == 1959);
j3 = find(ff1(:,1) == 1979);

% i4 = find(dtdelpCO2a(:,1) == 1958.5);
% j4 = find(dtdelpCO2a(:,1) == 2005.5);
% i4 = find(dtdelpCO2a(:,1) == 1959);
% j4 = find(dtdelpCO2a(:,1) == 1979);
i4 = find(floor(100*dtdelpCO2a(:,1)) == floor(100*(1959+(1/24))));
j4 = find(floor(100*dtdelpCO2a(:,1)) == floor(100*(1979+(1/24))));
% 

x1 = nino3(i2-3:j2-4,2);
x2 = tau(1303:1867-1,2);
%x3 = ff1(i3:j3-1,2);
x3 = ff1(i3:j3,2);

%d = dtdelpCO2a(i4:j4-1,2);
d = dtdelpCO2a(i4:j4,2); %CO2 record from 1959-1979

%G = [ones(size(x3)) x3 x1 x2]; 

G2 = [ones(size(x3)) x3]; % fossil fuels from 1959-1979

%m = inv(G'*G)*G'*d; 

%m2 = G\d; % regression of co2 increment using all variables

m3 = G2\d; % regression of ff from 1959-1979
%m3 = [-0.0259;0.57469];
x4 = ff1(:,2);

G3 = [ones(size(x4)) x4];

% Regressions

% dm2(:,1) = ff1(i3:j3-1,1);
% dm2(:,2) = G*m2; % all variables

dm3(:,1) = ff1(i3:j3,1);
dm3(:,2) = G2*m3; % fossil fuels only

dm4(:,1) = ff1(:,1);
dm4(:,2) = G3*m3; % whole fossil fuel record used

% figure
% plot(dtdelpCO2a(i4:j4-1,1),d,dm2(:,1),dm2(:,2),dm3(:,1),dm3(:,2))

% subtract out regressed NINO3 and tau parts from CO2 record

% canom(:,1) = ff1(i3:j3-1,1);
% canom(:,2) = d - m2(3)*x1 - m2(4)*x2;

% subtract out all in regression

% canom2(:,1) = ff1(i3:j3-1,1);
% canom2(:,2) = d - dm2(:,2);

% figure
% plot(dtdelpCO2a(i4:j4-1,1),d,canom(:,1),canom(:,2))

% Integrals

[ffi] = integrate_series_trap(dm4,1,2,12);

%[alli] = integrate_series_trap(dm2,1,2,12);

%[co2i] = integrate_series_trap(dtdelpCO2a(i4:j4-1,:),1,2,12); 

[co2i] = integrate_series_trap(dtdelpCO2a,1,2,12); % integrate whole record

%[anomi] = integrate_series_trap(canom,1,2,12); %integrate anomaly


%a = find(floor(100*CO2a(:,1)) == floor(100*(1958.5)))
a = find(floor(100*CO2a(:,1)) == floor(100*(1958.5+(1/24))))
b = CO2a(a,2)
cc = find(ffi(:,1) == 1958.5);
ffi(:,2) = ffi(:,2) - ffi(cc,2) + b;

% dd = find(alli(:,1) == 1958.5);
% alli(:,2) = alli(:,2) - alli(dd,2) + b;
% 
%e = find(floor(100*co2i(:,1)) == floor(100*(1958.5)));
e = find(floor(100*co2i(:,1)) == floor(100*(1958.5+(1/24))));
co2i(:,2) = co2i(:,2) - co2i(e,2) + b;

% 
%f = find(anomi(:,1) == 1958.5);
% anomi(:,2) = anomi(:,2) - anomi(f,2) + b;

% Subtract from CO2 record to get anomaly

% i6 = find(CO2a(:,1) == 1958.5);
% j6 = find(CO2a(:,1) == 2005.5);
% 
% CO2_anom(:,1) = CO2a(i6:j6-1,1);
% CO2_anom(:,2) = CO2a(i6:j6-1,2) - ffi(i3:j3-1,2);
% 
% CO2_anoma(:,1) = CO2a(i6:j6-1,1);
% CO2_anoma(:,2) = CO2a(i6:j6-1,2) - alli(:,2);
% 
% co2ia(:,1) = co2i(i4:j4-1,1);
% co2ia(:,2) = co2i(i4:j4-1,2) - ffi(i3:j3-1,2); % integrated CO2 - ff

% co2iaa(:,1) = co2i(i4:j4-1,1);
% %co2iaa(1:i4-1,2) = co2i(:,2);
% co2iaa(:,2) = co2i(i4:j4-1,2) - alli(:,2); % integrated CO2 minus integrated regression

% co2iaai(:,1) = co2iaa(:,1);
% co2iaai(:,2) = co2iaa(:,2) + ffi(i3:j3-1,2); % integrated anomaly plus ff

% figure
% plot(CO2a(i6:j6-1,1),CO2a(i6:j6-1,2),ffi(:,1),ffi(:,2),alli(:,1),alli(:,2))
% title('real record')
% 
% figure
% plot(CO2_anom(:,1),CO2_anom(:,2),CO2_anoma(:,1),CO2_anoma(:,2))
% title('real record')
% 
% figure
% plot(co2i(:,1),co2i(:,2),ffi(:,1),ffi(:,2),alli(:,1),alli(:,2))
% title('calc record')
% % 
% figure
% plot(co2ia(:,1),co2ia(:,2),co2iaa(:,1),co2iaa(:,2))
% title('calc record')
% 
% figure
% plot(CO2a(i6:j6-1,1),CO2a(i6:j6-1,2),co2i(:,1),co2i(:,2))

% boxcar average the temperature

[avg_temp] = l_boxcar(tland4,1,12,1,2483,1,2);

avg_temp(1:6,2) = avg_temp(7,2);

[avg_temp2] = l_boxcar(tland4,10,12,1,2483,1,2);

% Calculate "predicted" anomaly

year2 = year';
[delC1dt2,delC2dt2,delCdt2, delC12,delC22] = bioboxtwo_sub7(0.8,1,1,1,ts,year2,dpCO2a,avg_temp);


% figure
% plot(delC1dt2(:,1),-delC1dt2(:,2))
% title('flux to atm from box 1')
% ylabel('ppm/yr')
% 
% figure
% plot(delC2dt2(:,1),-delC2dt2(:,2))
% title('flux to atm from box 2')
% ylabel('ppm/yr')

C12(:,2) = 110 + 2.12*delC12(:,2);

C22(:,2) = 1690+ 2.12*delC22(:,2);

% figure
% plot(delC12(:,1),C12(:,2))
% title('box 1 size')
% ylabel('GtC')
% 
% figure
% plot(delC22(:,1),C22(:,2))
% title('box 2 size')
% ylabel('GtC')
% figure
% plot(delCdt2(:,1),-delCdt2(:,2))

co2p(:,1) = year(1,601:2467);
co2p(:,2) = ff1(1189:3055,2) + landusemo(1:1867,2) -  Aoc*fas(601:2467,2) - delCdt2(601:2467,2);

% figure
% plot(ff1(1189:3055,1),ff1(1189:3055,2),landusemo(1:1867,1),landusemo(1:1867,2),fas(601:2467,1),-Aoc*fas(601:2467,2),delCdt2(601:2467,1),-delCdt2(601:2467,2))
% 
% figure
% plot(co2p(:,1),co2p(:,2),dtdelpCO2a(i4:j4-1,1),d)

[predict_int] = integrate_series_trap(co2p,1,2,12);

i8 = find(predict_int(:,1) == 1958.5);
j8 = find(predict_int(:,1) == 2005.5);

predict_int(:,2) = predict_int(:,2) - predict_int(i8,2) + b;

% panom(:,1) = predict_int(i8:j8-1,1);
% panom(:,2) = predict_int(i8:j8-1,2) - ffi(i3:j3-1,2);

%load atm_forward.mat

%load forward_run3_6land2.mat

% load forward_run3_8.mat

%load forward_run3_21_lowT_2.mat

%load forward_run4_25_noT.mat

%load forward_eps9.mat

%load forward_varT_lowLU_10.mat

%load forward_CVHD_9.mat

%load meure_CHMV_1900_1976.mat

load CHMC_2box.mat

%load undecon_CHMC_meure_10apr.mat

%load forward_run4_25.mat

%%-------------------------------------%%
% %redo for modern ff
% clear x4 G3 dm4 ffi
% x4 = ff1(:,2);
% 
% G3 = [ones(size(x4)) x4];
% dm4(:,1) = ff1(:,1);
% dm4(:,2) = G3*m3; % whole fossil fuel record used
% 
% [ffi] = integrate_series_trap(dm4,1,2,12);
% a = find(floor(100*CO2a(:,1)) == floor(100*(1958.5+(1/24))))
% b = CO2a(a,2)
% cc = find(ffi(:,1) == 1958.5);
% ffi(:,2) = ffi(:,2) - ffi(cc,2) + b;
%%-------------------------------------%%

i = find(datm(:,2) == 0);
datm(i,2) = NaN;

i9 = find(datm(:,1) == 1958.5);
j9 = find(datm(:,1) == 1978.5);

i12 = find(datm(:,1) == 1959);
j12 = find(datm(:,1) == 1979);

X = datm(i12:j12,2);


% i13 = find(floor(100*co2i(:,1)) == floor(100*(1959)));
% j13 = find(floor(100*co2i(:,1)) == floor(100*(1979)));
i13 = find(floor(100*co2i(:,1)) == floor(100*(1959+(1/24))));
j13 = find(floor(100*co2i(:,1)) == floor(100*(1979+(1/24))));

y6 = co2i(i13:j13,2);

beta = 288;

[betahat,resid,J] = nlinfit(X,y6,'forward_run_fit',beta);

betahat

%p = polyfit(X,y6,0)

datm3(:,1) = datm(:,1);
datm3(:,2) = datm(:,2) - datm(i9,2) + b;

%G10 = [ones(size(x6)) x6];

%m6 = G10\y6;

datm(:,2) = datm(:,2) + betahat;

%datm(:,2) = datm(:,2) + 290;

x4 = datm(i9:j9,2);
x5 = datm(i9:end,2);

clear e

j10 = find(floor(100*co2i(:,1)) == floor(100*(1978.5+(1/24))));
e = find(floor(100*co2i(:,1)) == floor(100*(1958.5+(1/24))));
% j10 = find(floor(100*co2i(:,1)) == floor(100*(1978.5)));
% e = find(floor(100*co2i(:,1)) == floor(100*(1958.5)));
y2 = co2i(e:j10,2);

G3 = [ones(size(x4)) x4];
G4 = [ones(size(x5)) x5];

m4 = G3\y2;

datm2(:,1) = datm(i9:end,1);
datm2(:,2) = G4*m4;

j11 = find(floor(100*co2i(:,1)) == floor(100*(1800+(1/24))));
%j11 = find(floor(100*co2i(:,1)) == floor(100*(1800)));
for i = 1201:2423
    
   a = find(floor(100*datm(:,1)) == floor(100*avg_temp2(i,1)));
   %b = find(floor(100*co2i(:,1)) == floor(100*avg_temp2(i,1)));
   b = find(floor(100*co2i(:,1)) == floor(100*(avg_temp2(i,1)+(1/24))));
   temp_diff(i-1200,1) = avg_temp2(i,2)-8.5;
   temp_diff(i-1200,2) = co2i(b,2)-datm(a,2);
     
end

j14 = find(floor(100*co2i(:,1)) == floor(100*(1850+(1/24))));
%j14 = find(floor(100*co2i(:,1)) == floor(100*(1850)));

[p,s,mu] = polyfit(temp_diff(:,1),temp_diff(:,2),1);
model = polyval(p,temp_diff(:,1),s,mu);

e = model - temp_diff(:,2);
misfit = e'*e/length(model)

% 
airfranom(:,1) = datm(:,1);
airfranom(:,2) = datm(:,2) - ffi(1189:3055,2);  % use for normal cases
%airfranom(:,2) = datm(:,2) - ffi(1189:3025,2);  % use for var SST cases
%airfranom(:,2) = datm(:,2) - ffi(1189:3078,2);  % use for updated (2007) cases

% airfranom(:,1) = datm(:,1);
% airfranom(:,2) = datm(:,2) - ffi(1189:3049,2);

airfranom_obs(:,1) = co2i(j14:4387,1);
airfranom_obs(:,2) = co2i(j14:4387,2) - ffi(1189:3055,2);

% airfranom_obs(:,1) = co2i(j14:4410,1);
% airfranom_obs(:,2) = co2i(j14:4410,2) - ffi(1189:3078,2);


for i = 1:365%((ts/2)+1):(length(years)-(ts/2))
    n = i*12;
    annaveMLOSPO(i,1) = CO2a(n,1);
    annaveMLOSPO(i,2) = sum(CO2a(n-(ts/2):n+(ts/2),2))/13;
end  

airfranom_year(:,1) = annaveMLOSPO(112:end,1);
for i = 112:length(annaveMLOSPO)
    j = find(ffi(:,1) == annaveMLOSPO(i,1)-(annaveMLOSPO(112,1)-ffi(12,1)));
    airfranom_year(i-111,2) = annaveMLOSPO(i,2) - ffi(j,2);
end

%load meure_data.mat

%airfranom_data(:,1) = icecore_data(:,1);

% for i = 56:length(icecore_data)
%     if(icecore_data(i,1) - floor(icecore_data(i,1)) == 0)
%         q = find(ffi(:,1) == icecore_data(i,1));
%         airfranom_data(i,2) = icecore_data(i,2) - ffi(q,2);
%     end
% end

% figure
% plot(temp_diff(:,1),temp_diff(:,2),'.',temp_diff(:,1),model,'-')
% xlabel('Land temperature (degrees C)')
% ylabel('ppm CO2')
% title('Difference between CO2 record and temperature independent model')

% figure
% plot(co2i(:,1),co2i(:,2),ffi(:,1),ffi(:,2),alli(:,1),alli(:,2),predict_int(:,1),predict_int(:,2))
% % title('calc record')
% 
% figure
% plot(co2i(j11:end,1),co2i(j11:end,2),ffi(:,1),ffi(:,2),alli(:,1),alli(:,2),datm(:,1),datm(:,2))
% title('forward run: slow box has no temp effects')
% legend('measured','fossil fuel','regression','model','Location','Northwest')
% xlabel('Year')
% ylabel('ppm C')

% figure
% plot(co2i(j11:end,1),co2i(j11:end,2),datm(:,1),datm(:,2),'r')
% title('forward run: fast box has no temp effects')
% legend('measured','model','Location','Northwest')
% xlabel('Year')
% ylabel('ppm C')

figure
% plot(co2i(j11:end,1),co2i(j11:end,2),ffi(:,1),ffi(:,2),co2iaai(:,1),co2iaai(:,2),datm(:,1),datm(:,2))
% title('forward run, low land use, temp')
% legend('measured','fossil fuel','detrended','model','Location','Northwest')
% xlabel('Year')
% ylabel('ppm C')

subplot('Position',[0.075 0.5 0.9 0.4])%(2,1,1)
plot(co2i(j11:end,1),co2i(j11:end,2),ffi(:,1),ffi(:,2),datm(:,1),datm(:,2))
title('forward run, CLM-CN')
legend('measured','fossil fuel','model','Location','Southeast')
xlabel('Year')
ylabel('ppm C')
set(gca,'Xlim',[1850 2010])
set(gca,'Ylim',[280 380])
set(gca,'Xminortick','on')

subplot('Position',[0.075 0.1 0.9 0.4])%(2,1,2,'align')
plot(airfranom(:,1),airfranom(:,2),'-r',airfranom_obs(:,1),airfranom_obs(:,2),'-b')
%legend('modeled','measured','Location','Northwest')
xlabel('Year')
ylabel('ppm C')
set(gca,'Xlim',[1850 2010])
set(gca,'Xminortick','on')

% subplot('Position',[0.2 0.6 0.2 0.2])
% plot(datm(:,1),datm(:,2))
% 

% figure
% plot(co2i(j11:end,1),co2i(j11:end,2),ffi(:,1),ffi(:,2),co2iaai(:,1),co2iaai(:,2),datm3(:,1),datm3(:,2))
% title('forward run, low land use, temp')
% legend('measured','fossil fuel','detrended','model','Location','Northwest')
% xlabel('Year')
% ylabel('ppm C')
% 
% figure
% plot(co2i(e:end,1),co2i(e:end,2),ffi(cc:end,1),ffi(cc:end,2),co2iaai(:,1),co2iaai(:,2),datm3(i9:end,1),datm3(i9:end,2))
% title('forward run, low land use, temp')
% legend('measured','fossil fuel','detrended','model','Location','Northwest')
% xlabel('Year')
% ylabel('ppm C')

% figure
% plot(co2i(:,1),co2i(:,2),ffi(:,1),ffi(:,2),alli(:,1),alli(:,2),datm2(:,1),datm2(:,2))
% title('forward run shifted')
% legend('measured','fossil fuel','regression','model','Location','Northwest')
% xlabel('Year')
% ylabel('ppm C')
% 
% figure
% plot(co2i(:,1),co2i(:,2),ffi(:,1),ffi(:,2),co2iaai(:,1),co2iaai(:,2),datm2(:,1),datm2(:,2))
% title('forward run shifted')
% legend('measured','fossil fuel','detrended','model','Location','Northwest')
% xlabel('Year')
% ylabel('ppm C')

% figure
% plot(co2ia(:,1),co2ia(:,2),co2iaa(:,1),co2iaa(:,2),panom(:,1),panom(:,2))
% title('Anomalies')
% legend('using fossil fuels','using ff, nino3 and tau','using model')
% xlabel('Year')
% ylabel('ppm C')

% Write to files

% fid = fopen('co2_full_record_meure_shift.txt','wt')
% fprintf(fid,'%9.4f\t %8.4f\n',co2i')
% 
% fid = fopen('co2_detrend2.txt','wt')
% fprintf(fid,'%9.4f\t %8.4f\n',co2iaai')
% 
% fid = fopen('ff_fit_all_meure.txt','wt')
% fprintf(fid,'%9.4f\t %8.4f\n',ffi')
% %  
% fid = fopen('co2_forward_atm.txt','wt')
% fprintf(fid,'%9.4f\t %8.4f\n',datm2')
%
% fid = fopen('CVLDN_undecon_meure_shift_fit.txt','wt')
% fprintf(fid,'%9.4f\t %8.4f\n',datm')
% % 
% fid = fopen('CHMC_2box_diff_shift_2007.txt','wt')
% fprintf(fid,'%9.4f\t %8.4f\n',airfranom')
% % 
% fid = fopen('obs_shift_diff_2007.txt','wt')
% fprintf(fid,'%9.4f\t %8.4f\n',airfranom_obs')

% fid = fopen('co2_forward10_no_bigT.txt','wt')
% fprintf(fid,'%9.4f\t %8.4f\n',datm')
% % 
% fid = fopen('co2_forward_CHMC_2box.txt','wt')
% fprintf(fid,'%9.4f\t %8.4f\n',datm')

% % 
% fid = fopen('CCHD_plusT.txt','wt')
% fprintf(fid,'%7.4f\t %7.4f\n',dtland10')

% fid = fopen('co2_minusff.txt','wt')
% fprintf(fid,'%7.4f\t %7.4f\n',CO2_anom')
% 
% fid = fopen('co2_anom_minusff.txt','wt')
% fprintf(fid,'%7.4f\t %7.4f\n',CO2_anoma')

% fid = fopen('atm_increment.txt','wt')
% fprintf(fid,'%7.4f\t %7.4f\n',dtdelpCO2a')
% 
% fid = fopen('model_increment.txt','wt')
% fprintf(fid,'%7.4f\t %7.4f\n',datmdt')
% 
% fid = fopen('anom_increment.txt','wt')
% fprintf(fid,'%7.4f\t %7.4f\n',canom')

% status = fclose('all');