%% Modifications
%% 11/27/2007: added NPP weighted land temperature datasets: npp_T.mat
%% 11/27/2007: added land weighted land T datasets: landwt_T.mat
%% 12/18/2007: added option to do fit with unfiltered data
%% 3/13/2008: changed data filter so that only 1957 to present is filtered
%% 4/2/2008: change index numbers for dtdelpCO2a to use with new dataset.
%% change 2053 to 2521, change 3919 to 4387
%% 3/31/2009: change index numbers to do run to 6/2007 for mode 1
%% 1/10/2011: add in joos_hilda_2011.mat
%% 8/23/2012: annotate code better

% forward run model with temperature dependence in land

clear all

for n = 1
%% define what kind of run you want to do

LU = n; %1 = high land use scenario; 2 = low land use scenario

nitrogen = 0; % 1 = yes, 0 = no; account for nitrogen fertilization?

filter = 1; % filter the data? 1 = 10 year filter; 2 = unfiltered

%%
load land_temp.mat % land temperature records

load npp_T.mat % NPP-weighted temperature record

load landwt_T_2011.mat % land temperature anomaly

%% load CO2 sources and sinks
%% landusemo: land use emissions, in ppm/year, interpolated to monthly
%% resolution
%% ff1: fossil fuel emissions
%% fas: ocean flux per m^2
%% Aoc: ocean surface area
%% extratrop_landmo: extratropical land use emissions

[landusemo,ff1,fas,Aoc,extratrop_landmo] = getsourcesink_scale3; 

 clear year start_year end_year ts

ts = 12; % timesteps per year
start_year = 1850;
end_year = 2009+(7/12); 

beta = [0.5;2]; % initial guesses for model fit

[dtdelpCO2a,dpCO2a,year,dt,CO2a] = MLOinterpolate_increment2(ts,start_year,end_year); 

% Extend land use record by making recent emissions equal to last
% record
landusemo(1874:1916,1) = year(1874:1916);
landusemo(1874:1916,2) = landusemo(1873,2);
% % 
%Extend extratropical emissions by assuming emissions are zero
extratrop_landmo(1802:1916,1) = landusemo(1802:1916,1);
extratrop_landmo(1802:1916,2) = 0;

%% Calculate residual land uptake
% run to 8/2009
% using high land use emissions
residual(:,1) = year(1,1:1916);
residual(:,2) = dtdelpCO2a(2521:4436,2) - ff1(1189:3104,2)....
+ Aoc*fas(601:2516,2) - landusemo(1:1916,2);

% using extratropical emissions only
residual2(:,1) = year(1,1:1916);
residual2(:,2) = dtdelpCO2a(2521:4436,2) - ff1(1189:3104,2)....
+ Aoc*fas(601:2516,2) - extratrop_landmo(1:1916,2);

%% tland4: the temperature record started at 1880. tland4 extends the
%% record back to 1800 by using the mean temperature for the first year
%% 
%% do a moving boxcar average of the land temperature: 1 year average
%% note: in this case the box length (1; second term in l_boxcar) is in
%% units of years. dt (12, third term) is the number of points per year
%% first column of avg_temp gives the date, second column gives the moving
%% average of the land temperature
[avg_temp] = l_boxcar(tland4,1,12,1,2483,1,2); 

avg_temp(1:6,2) = avg_temp(7,2); % make the first 6 points 

%% 10 year moving boxcar average of land temperature
% [avg_temp] = l_boxcar(tland4,10,12,1,2483,1,2);
% 
% avg_temp(1:60,2) = avg_temp(61,2);

%%----------------------------------------------------------------------%%
%% Pick the temperature record to use
%%----------------------------------------------------------------------%%

%%----------------------------------------------------------------------%%
% 
% ***Use these for various T tests***

 temp_anom(1:6,1) =  avg_temp(601:606,1); %Jan 1850-May 1850
 temp_anom(1:6,2) = landtglob(1,2); %355 instead of 1, 360 instead of 6
  
 temp_anom(7:1916,1) = landtglob(1:1910,1); % Starts at the year 1850.5. 
 temp_anom(7:1916,2) = landtglob(1:1910,2); % 
 
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
 

%%----------------------------------------------------------------------%%

%% calculate a 10-year running boxcar average of the residual land uptake
%% don't use 10 year mean before 1957, because data are already smoothed
%% (ice core)
if(LU==1) %high land use
%[residual10] = l_boxcar(residual,10,12,1,length(residual),1,2);
[residual10a] = l_boxcar(residual,10,12,1225,length(residual),1,2);
residual10(1:1284,:) = residual(1:1284,:);
residual10(1285:(length(residual10a)),:) = residual10a(1285:end,:);

decon = residual;
elseif(LU ==2) % low land use
%[residual10] = l_boxcar(residual2,10,12,1,length(residual2),1,2);
[residual10a] = l_boxcar(residual2,10,12,1225,length(residual2),1,2);
residual10(1:1284,:) = residual2(1:1284,:);
residual10(1285:(length(residual10a)),:) = residual10a(1285:end,:);

decon = residual2;
end

%% find model fit using a nonlinear regression
if(filter == 1) % fit to 10-year filtered record

    [betahat,resid,J] = nlinfit(X,residual10(601:end,2),'land_fit_Qs_annotate',beta); %change 601:end to 1081:1513; change 601 to 1297

elseif(filter == 2) % fit to unfiltered record

    [betahat,resid,J] = nlinfit(X,decon(1057:1513,2),'land_fit_Qs_annotate',beta); %change 601:end to 1081:1513. After 1958: 1297; was on 1309:end

end

%% Look at covariances and correlations between model result and calculated land uptake 
[N,P] = size(X);
% 
covariance = inv(J(1:1177,:)'*J(1:1177,:))*sum(resid(1:1177,:).^2)/(N-P) 

[isize,jsize] = size(covariance);
for i=1:isize
    for j=1:jsize
correlation(i,j) = covariance(i,j)/sqrt(covariance(i,i)*covariance(j,j));
    end
end

%% Get uncertainties of best fit values
ci = nlparci(betahat,resid,J);

%% Redefine values of epsilon, gamma and Q1
if(nitrogen == 1)
epsilon = 0;% betahat(2)
gamma = betahat(1)
Q1 = betahat(2)
Q2 = 1;%betahat(3)
else
epsilon = betahat(1)%0.79;%
Q1 = betahat(2)%4.91;
Q2 = 1;%betahat(2)
end

year2 = year';

%% Run the best fit values in the model again to plot
 [C1dt,C2dt,delCdt,delC1,delC2] = bioboxtwo_sub10_annotate(epsilon,Q1,Q2,ts,year2,dpCO2a,X); 
 
%%% Nitrogen%%%

% [C1dt,C2dt,delCdt,delC1,delC2] = bioboxtwo_subN_annotate(epsilon,Q1,Q2,gamma,ff1,ts,year2,dpCO2a,X);
   
    
delCdt(:,2) = -delCdt(:,2);

%% 10 year moving boxcar average of model result
[delC10] = l_boxcar(delCdt,10,12,1,length(delCdt),1,2);

%% Quantify how good the fit is using the mean squared error (MSE)       
% Use filtered or unfiltered yhat
%----------------------------------------------------------------%
%
% Filtered
%
%------------------------------------------------

if(filter == 1)

yhat2 = delC10(:,2);

%yhat2 = yhat2 + (0 - yhat2(1));
    
    e = delC10(61:1806,2) - residual10(61:1806,2); % look at MSE for 1855 to 2000
    
    misfit = e'*e/length(delC10(61:1806,2));  
    
    e2 = delC10(601:1806,2) - residual10(601:1806,2); % look at MSE for 1900 to 2000
    misfit2 = e2'*e2/length(delC10(601:1806,2));  


   error1 = betahat(1)-ci(1);
   
%   error2= betahat(2)-ci(2)

%   error3= betahat(3)-ci(3)
   
%   figure
% plot(decon(:,1),decon(:,2),delCdt(:,1),yhat2)
% xlabel('year')
% ylabel('ppm CO2/year')
% title('land uptake')
% legend('Residual uptake','land uptake without T effects','land uptake with T effects')

   
    C = cov(residual10(601:(end-1),2),yhat2(601:(end-1),1));
%    
   [R,P,RLO,RUP] = corrcoef(yhat2(601:(end-1),1),residual10(601:(end-1),2));
%   
   R(1,2)^2;
   
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
  
      e = delCdt(61:1806,2) - decon(61:1806,2); % look at MSE for 1855-2000
    
   % misfit = e'*e/length(delCdt(61:1806,2))  
    
    e2 = delCdt(601:1806,2) - decon(601:1806,2); % look at MSE for 1900-2000
    
   % misfit2 = e2'*e2/length(delCdt(601:1806,2))   

    e3 = delCdt(1309:(end-1),2) - decon(1309:(end-1),2); % look at MSE from 1958-present
    misfit3 = e3'*e3/length(delCdt(1309:(end-1),2));  
    
      C = cov(decon(1309:(end-1),2),yhat2(1309:(end-1),1));
   
  [R,P,RLO,RUP] = corrcoef(yhat2(1297:(end-1),1),decon(1297:(end-1),2));
  
  %R(1,2)^2

  error1= betahat(1)-ci(1);
    
figure
plot(decon(:,1),decon(:,2),delCdt(:,1),yhat2)
xlabel('year')
ylabel('ppm CO2/year')
title('land uptake')
legend('Residual uptake','land uptake without T effects','land uptake with T effects')

end

    
end

%% Do "reverse deconvolution" to calculate modeled atmospheric change in
%% CO2
if(LU==1)
newat(:,1) = year(1,1:1916);
newat(:,2) =  ff1(1189:3104,2)....
- Aoc*fas(601:2516,2) + landusemo(1:1916,2) + delCdt(:,2) ;

elseif(LU==2)
newat(:,1) = year(1,1:1867);
newat(:,2) =  ff1(1189:3055,2)....
- Aoc*fas(601:2467,2) + extratrop_landmo(1:1867,2) + delCdt(:,2) ;
end

[datm] = integrate_series_trap(newat,1,2,12);
