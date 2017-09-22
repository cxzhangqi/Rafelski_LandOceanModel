%% Modifications
%% 11/27/2007: added NPP weighted land temperature datasets: npp_T.mat
%% 11/27/2007: added land weighted land T datasets: landwt_T.mat
%% 12/18/2007: added option to do fit with unfiltered data
%% 3/13/2008: changed data filter so that only 1957 to present is filtered
%% 4/2/2008: change index numbers for dtdelpCO2a to use with new dataset.
%% change 2053 to 2521, change 3919 to 4387
%% 3/31/2009: change index numbers to do run to 6/2007 for mode 1

% forward run model with temperature dependence in land

clear all

for s = 1
    for n = 1
mode = s; % 1 = standard cases; 2 = high ocean uptake; 3 = low ocean uptake; 4 = SST

LU = n; %1 = high land use; 2 = low land use

nitrogen = 0; % 1 = yes, 0 = no

load land_temp.mat

load npp_T.mat

load landwt_T.mat

%load landwt_T_2009.mat

%load joos_rkDIC_varT_water.mat

load joos_rkDIC_varT_water_3.mat

if (mode == 2 | mode == 3)
load joos_hilda_bigsm_scale.mat
elseif (mode == 4)
load joos_rkDIC_varT_water_3.mat
end
% 
[landusemo,ff1,fas,Aoc,extratrop_landmo] = getsourcesink_scale2;

%clear ff1 fas Aoc

%load joos_hilda_2007.mat


 clear year start_year end_year ts

ts = 12;
start_year = 1850;

if (mode == 1)
   end_year = 2005.5;
%end_year = 2007+(5/12);%.5; %from 2005.5
elseif (mode == 2 | mode ==3)
end_year = 2005; % high and low sabine
elseif (mode == 4)
end_year = 2003;
end

beta = [0.5;2];
%beta = [0.5];

filter = 1; % 1 = 10 year filter; 2 = unfiltered

[dtdelpCO2a,dpCO2a,year,dt,CO2a] = MLOinterpolate_increment(ts,start_year,end_year);



% clear ff1
% load ff3.mat
% 
% ff1 = ff4;

% landusemo(1874:1890,1) = year(1874:1890);
% landusemo(1874:1890,2) = landusemo(1873,2);
% 
extratrop_landmo(1802:1873,1) = landusemo(1802:1873,1);
extratrop_landmo(1802:1873,2) = 0;

% extratrop_landmo(1802:1890,1) = landusemo(1802:1890,1);
% extratrop_landmo(1802:1890,2) = 0;

if (mode == 1)
    
% residual(:,1) = year(1,1:1890);
% residual(:,2) = dtdelpCO2a(2521:4410,2) - ff1(1189:3078,2)....
% + Aoc*fas(601:2490,2) - landusemo(1:1890,2);
% 
% residual2(:,1) = year(1,1:1890);
% residual2(:,2) = dtdelpCO2a(2521:4410,2) - ff1(1189:3078,2)....
% + Aoc*fas(601:2490,2) - extratrop_landmo(1:1890,2);
    
residual(:,1) = year(1,1:1867);
residual(:,2) = dtdelpCO2a(2521:4387,2) - ff1(1189:3055,2)....
+ Aoc*fas(601:2467,2) - landusemo(1:1867,2);

residual2(:,1) = year(1,1:1867);
residual2(:,2) = dtdelpCO2a(2521:4387,2) - ff1(1189:3055,2)....
+ Aoc*fas(601:2467,2) - extratrop_landmo(1:1867,2);

elseif (mode == 2)
    
residual(:,1) = year(1,1:1861);
residual(:,2) = dtdelpCO2a(2521:4381,2) - ff1(1189:3049,2)....
+ AocH*fas_h_big(601:2461,2) - landusemo(1:1861,2);  % 2521 to 2053; 4387 to 3919

residual2(:,1) = year(1,1:1861);
residual2(:,2) = dtdelpCO2a(2521:4381,2) - ff1(1189:3049,2)....
+ AocH*fas_h_big(601:2461,2) - extratrop_landmo(1:1861,2);

elseif (mode == 3)
 
residual(:,1) = year(1,1:1861);
residual(:,2) = dtdelpCO2a(2521:4381,2) - ff1(1189:3049,2)....
+ AocH*fas_h_sm(601:2461,2) - landusemo(1:1861,2);

residual2(:,1) = year(1,1:1861);
residual2(:,2) = dtdelpCO2a(2521:4381,2) - ff1(1189:3049,2)....
+ AocH*fas_h_sm(601:2461,2) - extratrop_landmo(1:1861,2);

elseif (mode == 4)
% residual(:,1) = year(1,1:1837);
% residual(:,2) = dtdelpCO2a(2521:4387,2) - ff1(1189:3025,2)....
% + Aoc_HvarT*fas_HvarT_wat(601:2437,2) - landusemo(1:1837,2);
% 
% residual2(:,1) = year(1,1:1837);
% residual2(:,2) = dtdelpCO2a(2521:4387,2) - ff1(1189:3025,2)....
% + Aoc_HvarT*fas_HvarT_wat(601:2437,2) - extratrop_landmo(1:1837,2);

residual(:,1) = year(1,1:1837);
residual(:,2) = dtdelpCO2a(2521:4357,2) - ff1(1189:3025,2)....
+ Aoc_HvarT*fas_HvarT(601:2437,2) - landusemo(1:1837,2);

residual2(:,1) = year(1,1:1837);
residual2(:,2) = dtdelpCO2a(2521:4357,2) - ff1(1189:3025,2)....
+ Aoc_HvarT*fas_HvarT(601:2437,2) - extratrop_landmo(1:1837,2);

end

[avg_temp] = l_boxcar(tland4,1,12,1,2483,1,2);

avg_temp(1:6,2) = avg_temp(7,2);

% [avg_temp] = l_boxcar(tland4,10,12,1,2483,1,2);
% 
% avg_temp(1:60,2) = avg_temp(61,2);

%%----------------------------------------------------------------------%%
%% Pick the temperature record to use
%%----------------------------------------------------------------------%%


% temp_anom(:,1) = avg_temp(:,1);  % Starts at the year 1800
% temp_anom(:,2) = avg_temp(:,2) - 8.5;

% X = temp_anom(601:2467,:); % for regular 1yr run, start at 1850

% X = temp_anom(601:2461,:); % for high and low sabine

% X = temp_anom(601:2437,:); % for variable SST

%%----------------------------------------------------------------------%%
% 
%  temp_anom(1:6,1) =  avg_temp(601:606,1);
%  temp_anom(1:360,2) = avg_temp(1,2)-8.5;%npptglob(355,2); %355 instead of 1, 360 instead of 6
%  
%  temp_anom(7:1867,1) = npptglob(1:1861,1); % Starts at the year 1850.5
%  temp_anom(361:1867,2) = npptglob(355:1861,2); % 361 instead of 7, 355 instead of 1
%  
%  temp_anom(1:6,1) =  avg_temp(601:606,1);
%  temp_anom(1:6,2) = npptno(1,2);
% 
%  temp_anom(7:1867,1) = npptno(1:1861,1); % Starts at the year 1850.5
%  temp_anom(7:1867,2) = npptno(1:1861,2);
 
%  temp_anom(1:6,1) =  avg_temp(601:606,1);
%  temp_anom(1:15,2) = npptso(10,2);
% 
%  temp_anom(7:1867,1) = npptso(1:1861,1); % Starts at the year 1850.5
%  temp_anom(16:1867,2) = npptso(10:1861,2);
 
%  temp_anom(1:6,1) =  avg_temp(601:606,1);
%  temp_anom(1:18,2) = nppttro(13,2);
% 
%  temp_anom(7:1867,1) = nppttro(1:1861,1); % Starts at the year 1850.5
%  temp_anom(19:1867,2) = nppttro(13:1861,2);

%%----------------------------------------------------------------------%%
% 
% ***Use these for various T tests***

 temp_anom(1:6,1) =  avg_temp(601:606,1);
 temp_anom(1:6,2) = landtglob(1,2); %355 instead of 1, 360 instead of 6
 
 temp_anom(7:1867,1) = landtglob(1:1861,1); % Starts at the year 1850.5. 1867, 1861 = 2005.5, 1891, 1885 = 2007.5
 temp_anom(7:1867,2) = landtglob(1:1861,2); % 361 instead of 7, 355 instead of 1
 
%  [p1,s1,mu1] = polyfit(temp_anom(601:end,1) , temp_anom(601:end,2),1);
%  tfit = polyval(p1,temp_anom(:,1),s1,mu1);
%  
%   temp_anom(:,2) = temp_anom(:,2) - tfit;
%  
 
%  temp_anom(1:6,1) =  avg_temp(601:606,1);
%  temp_anom(1:6,2) = npptglob(1,2); %355 instead of 1, 360 instead of 6
%  
%  temp_anom(7:1867,1) = npptglob(1:1861,1); % Starts at the year 1850.5
%  temp_anom(7:1867,2) = npptglob(1:1861,2); % 361 instead of 7, 355 instead of 1
%  
%  temp_anom(1:6,1) =  avg_temp(601:606,1);
%  temp_anom(1:6,2) = npptno(1,2);
% 
%  temp_anom(7:1867,1) = npptno(1:1861,1); % Starts at the year 1850.5
%  temp_anom(7:1867,2) = npptno(1:1861,2);
%  
%  temp_anom(1:6,1) =  avg_temp(601:606,1);
%  temp_anom(1:15,2) = npptso(10,2);
% 
%  temp_anom(7:1867,1) = npptso(1:1861,1); % Starts at the year 1850.5
%  temp_anom(16:1867,2) = npptso(10:1861,2);
%  
%  temp_anom(1:6,1) =  avg_temp(601:606,1);
%  temp_anom(1:18,2) = nppttro(13,2);
% 
%  temp_anom(7:1867,1) = nppttro(1:1861,1); % Starts at the year 1850.5
%  temp_anom(19:1867,2) = nppttro(13:1861,2);

%  temp_anom(:,1) = T(601:2467,1);
%  temp_anom(:,2) = T(601:2467,2);
%  
X = temp_anom(:,:);
 
%X = temp_anom(1:1861,:); % for high and low sabine

%X = temp_anom(1:1837,:); % for variable SST

%X = temp_anom(601:end,2);

%%----------------------------------------------------------------------%%

if(LU==1)
%[residual10] = l_boxcar(residual,10,12,1,length(residual),1,2);
[residual10a] = l_boxcar(residual,10,12,1225,length(residual),1,2);
residual10(1:1284,:) = residual(1:1284,:);
residual10(1285:(length(residual10a)),:) = residual10a(1285:end,:);

decon = residual;
elseif(LU ==2)
%[residual10] = l_boxcar(residual2,10,12,1,length(residual2),1,2);
[residual10a] = l_boxcar(residual2,10,12,1225,length(residual2),1,2);
residual10(1:1284,:) = residual2(1:1284,:);
residual10(1285:(length(residual10a)),:) = residual10a(1285:end,:);

decon = residual2;
end



if(filter == 1)

    [betahat,resid,J] = nlinfit(X,residual10(601:end,2),'land_fit_Qs',beta); %change 601:end to 1081:1513; change 601 to 1297

elseif(filter == 2)

    [betahat,resid,J] = nlinfit(X,decon(1057:1513,2),'land_fit_Qs',beta); %change 601:end to 1081:1513. After 1958: 1297; was on 1309:end

end

[N,P] = size(X);
% 
covariance = inv(J(1:1177,:)'*J(1:1177,:))*sum(resid(1:1177,:).^2)/(N-P) %1206 to 1200 to 1177

[isize,jsize] = size(covariance);
for i=1:isize
    for j=1:jsize
correlation(i,j) = covariance(i,j)/sqrt(covariance(i,i)*covariance(j,j));
    end
end

correlation

%[betahat,resid,J] = nlinfit(X,landint(1:1823,2),'land_fit',beta);

ci = nlparci(betahat,resid,J);

% 
if(nitrogen == 1)
epsilon = 0;% betahat(2)
gamma = betahat(1)
Q1 = betahat(2)
Q2 = 1;%betahat(3)
%C1 = betahat(3)
else
epsilon = betahat(1)
Q1 = betahat(2)
Q2 = 1;%betahat(2)
%C1 = betahat(3)
end


% epsilon = betahat(1)
% gamma = betahat(2)
% Q1 = betahat(3)
% Q2 = betahat(4)

year2 = year';

%temp_model = avg_temp(601:end,:);

% [C1dt,delCdt,delC1] = bioboxone_sub10(epsilon,Q1,ts,year2,dpCO2a,X);

 [C1dt,C2dt,delCdt,delC1,delC2] = bioboxtwo_sub10(epsilon,Q1,Q2,ts,year2,dpCO2a,X); 
 
%  [C1dt,delCdt,delC1] = bioboxone_boxsize(epsilon,Q1,C1,ts,year2,dpCO2a,X);
  %  [C1dt,delCdt,delC1] = bioboxone_boxsize(epsilon,Q1,gamma,C1,ff1,ts,year2,dpCO2a,X);
  
%[C1dt,C2dt,delCdt,delC1,delC2] = bioboxtwo_boxsize(epsilon,Q1,Q2,C1,ts,year2,dpCO2a,X);
  

%%% Nitrogen%%%

% [C1dt,delCdt,delC1] = bioboxone_subN(epsilon,Q1,gamma,ff1,ts,year2,dpCO2a,X);

% [C1dt,C2dt,delCdt,delC1,delC2] = bioboxtwo_subN(epsilon,Q1,Q2,gamma,ff1,ts,year2,dpCO2a,X);
   
%[C1dt,C2dt,delCdt,delC1,delC2] = bioboxtwo_Nsize(epsilon,Q1,Q2,gamma,C1,ff1,ts,year2,dpCO2a,X);
   


   % calculate land uptake, temperature input doesn't matter because don't use it
  
delCdt(:,2) = -delCdt(:,2);

%[landint_noT] = integrate_series_trap(delCdt,1,2,12);

    [delC10] = l_boxcar(delCdt,10,12,1,length(delCdt),1,2);
%     [delC10a] = l_boxcar(delCdt,10,12,1225,length(delCdt),1,2);
%     delC10(1:1284,:) = delCdt(1:1284,:);
%     delC10(1285:length(delC10a),:) = delC10a(1285:end,:);

   
    
% Use filtered or unfiltered yhat
%----------------------------------------------------------------%
%
% Filtered
%
%------------------------------------------------

if(filter == 1)

yhat2 = delC10(:,2);

%yhat2 = yhat2 + (0 - yhat2(1));

if (mode == 1)
    
    
    e = delC10(61:1806,2) - residual10(61:1806,2); % look at MSE 
    
    misfit = e'*e/length(delC10(61:1806,2))  % chi square?
    
    e2 = delC10(601:1806,2) - residual10(601:1806,2); % look at MSE. change 1806 to 1081. change 601 to 1297
    misfit2 = e2'*e2/length(delC10(601:1806,2))  % chi square?

%     e2 = delC10(601:1513,2) - residual10(601:1513,2); % look at MSE. change 1806 to 1081. change 601 to 1297
%     misfit2 = e2'*e2/length(delC10(601:1513,2))  % chi square?
    
elseif(mode == 2 | mode == 3)
    
    e = delC10(61:1800,2) - residual10(61:1800,2); % look at MSE 
    
    misfit = e'*e/length(delC10(61:1800,2))  % chi square?
    
    e2 = delC10(601:1800,2) - residual10(601:1800,2); % look at MSE 
    
    misfit2 = e2'*e2/length(delC10(601:1800,2))  % chi square?
  
    
elseif(mode == 4)
    e = delC10(61:1776,2) - residual10(61:1776,2); % look at MSE % For var SST
    
    misfit = e'*e/length(delC10(61:1776,2))  % chi square?

    e2 = delC10(601:1776,2) - residual10(601:1776,2); % look at MSE 
    
    misfit2 = e2'*e2/length(delC10(601:1776,2))  % chi square?
    
end


   error1= betahat(1)-ci(1)
%   error2= betahat(2)-ci(2)

%     error3= betahat(3)-ci(3)
 
 %   yhat2 = delCdt(:,2);    
%    
%      e3 = delCdt(1297:(end-1),2) - decon(1297:(end-1),2); % look at MSE from 1958-present
%     misfit3 = e3'*e3/length(delCdt(1297:(end-1),2))  % chi square?    
%     
%      C = cov(decon(1297:(end-1),2),yhat2(1297:(end-1),1))
%    
%  [R,P,RLO,RUP] = corrcoef(yhat2(1297:(end-1),1),decon(1297:(end-1),2))
%   
% R2 = R(1,2)^2

  
%   figure
% plot(decon(:,1),decon(:,2),delCdt(:,1),yhat2)
% xlabel('year')
% ylabel('ppm CO2/year')
% title('land uptake')
% legend('Residual uptake','land uptake without T effects','land uptake with T effects')

   
    C = cov(residual10(601:(end-1),2),yhat2(601:(end-1),1))
%    
   [R,P,RLO,RUP] = corrcoef(yhat2(601:(end-1),1),residual10(601:(end-1),2));
%   
   R(1,2)^2
   
figure
plot(residual10(:,1),residual10(:,2),delC10(:,1),yhat2)
xlabel('year')
ylabel('ppm CO2/year')
title('land uptake')
legend('Residual uptake','land uptake without T effects','land uptake with T effects')
set(gca,'Xlim',[1850 2010])   
%      

%----------------------------------------------------------------%
%
% Unfiltered
%
%----------------------------------------------------------------%
elseif(filter == 2)

yhat2 = delCdt(:,2);    
  
      e = delCdt(61:1806,2) - decon(61:1806,2); % look at MSE 
    
   % misfit = e'*e/length(delCdt(61:1806,2))  % chi square?
    
    e2 = delCdt(601:1806,2) - decon(601:1806,2); % look at MSE. change 1806 to 1081
    
   % misfit2 = e2'*e2/length(delCdt(601:1806,2))  % chi square?    

    e3 = delCdt(1309:(end-1),2) - decon(1309:(end-1),2); % look at MSE from 1958-present
    misfit3 = e3'*e3/length(delCdt(1309:(end-1),2))  % chi square?    
    
      C = cov(decon(1309:(end-1),2),yhat2(1309:(end-1),1))
   
  [R,P,RLO,RUP] = corrcoef(yhat2(1297:(end-1),1),decon(1297:(end-1),2));
  
  %R(1,2)^2

  error1= betahat(1)-ci(1)
    
figure
plot(decon(:,1),decon(:,2),delCdt(:,1),yhat2)
xlabel('year')
ylabel('ppm CO2/year')
title('land uptake')
legend('Residual uptake','land uptake without T effects','land uptake with T effects')

end

% %------------------------------------------------------------------%
% 
% fid = fopen('CHMV_unfilt_meure_10yr.txt','wt');
% fprintf(fid,'%7.4f\t %7.4f\n', delCdt');
% fclose(fid);

% 
% fid = fopen('CHMV_unfilt_resid_meure.txt','wt');
% fprintf(fid,'%7.4f\t %7.4f\n', decon');
% fclose(fid);
% 
% % 
% fid = fopen('CLMC_meure_2box.txt','wt');
% fprintf(fid,'%7.4f\t %7.4f\n', delC10');
% fclose(fid);
% 
% fid = fopen('CLM_meure_decon_shift.txt','wt');
% fprintf(fid,'%7.4f\t %7.4f\n', residual10');
% fclose(fid);

    end
end

if(LU==1)
newat(:,1) = year(1,1:1867);
newat(:,2) =  ff1(1189:3055,2)....
- Aoc*fas(601:2467,2) + landusemo(1:1867,2) + delCdt(:,2) ;

% newat(:,1) = year(1,1:1890);
% newat(:,2) =  ff1(1189:3078,2)....
% - Aoc*fas(601:2490,2) + landusemo(1:1890,2) + delCdt(:,2) ;
elseif(LU==2)
newat(:,1) = year(1,1:1867);
newat(:,2) =  ff1(1189:3055,2)....
- Aoc*fas(601:2467,2) + extratrop_landmo(1:1867,2) + delCdt(:,2) ;
end


% if(LU==1)
% newat(:,1) = year(1,1:1837);
% newat(:,2) =  ff1(1189:3025,2)....
% - Aoc_HvarT*fas_HvarT(601:2437,2) + landusemo(1:1837,2) + delCdt(:,2) ;
% elseif(LU==2)
% newat(:,1) = year(1,1:1837);
% newat(:,2) =  ff1(1189:3025,2)....
% - Aoc_HvarT*fas_HvarT(601:2437,2) + extratrop_landmo(1:1837,2) + delCdt(:,2) ;
% end

[datm] = integrate_series_trap(newat,1,2,12);
% 
figure
plot(datm(:,1),datm(:,2)+295,CO2a(:,1),CO2a(:,2))
