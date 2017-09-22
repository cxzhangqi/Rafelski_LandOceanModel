% Things to check: end year chosen? lines 12-14
% beta variables defined correctly? lines 43-54
% correct model chosen? lines 58-69
% correct yhat? lines 75-77

function yhat = land_fit_Qs(beta,X)

% Initial conditions - set according to Joos et al

ts = 12;
start_year = 1850;
%end_year = 2007+(5/12);%.5;
end_year = 2005.5;
%end_year = 2005; % high and low sabine
%end_year = 2003;

% Get parameters

[dtdelpCO2a,dpCO2a,year,dt,CO2a] = MLOinterpolate_increment(ts,start_year,end_year);

[landusemo,ff1,fas,Aoc,extratrop_landmo] = getsourcesink_scale2;

% Set arrays to zero

% load land_temp.mat
% 
% [avg_temp2] = l_boxcar(tlanfigured4,10,12,1,2483,1,2);
% 
% temp_anom2(:,1) = avg_temp2(:,1);
% temp_anom2(:,2) = avg_temp2(:,2) - 8.5;
% 
% [avg_temp] = l_boxcar(tland4,1,12,1,2483,1,2);
% 
% temp_anom(:,1) = avg_temp(:,1);
% temp_anom(:,2) = avg_temp(:,2) - 8.5;
% 
% temp_model = avg_temp(601:end,:);

% X = temp_anom(601:2461,2);
% 
% beta = [0.9;18.2; 1];

epsilon = beta(1);
Q1 = beta(2);
Q2 = 1;%beta(2);
%C1 = beta(3);
% 
% epsilon = 0;
% gamma = beta(1);
% Q1 = beta(2);
% Q2 = 1;%beta(3);
% % C1 = beta(3);

% epsilon = beta(1);
% gamma = beta(2);
% Q1 = beta(3);
% Q2 = beta(4);

year2 = year';

%[C1dt,delCdt,delC1] = bioboxone_sub10(epsilon,Q1,ts,year2,dpCO2a,X);
  
[C1dt,C2dt,delCdt,delC1,delC2] = bioboxtwo_sub10(epsilon,Q1,Q2,ts,year2,dpCO2a,X); 

 % [C1dt,delCdt,delC1] = bioboxone_boxsize(epsilon,Q1,C1,ts,year2,dpCO2a,X);
  
% [C1dt,delCdt,delC1] = bioboxone_boxsize(epsilon,Q1,gamma,C1,ff1,ts,year2,dpCO2a,X);
   
% [C1dt,C2dt,delCdt,delC1,delC2] = bioboxtwo_boxsize(epsilon,Q1,Q2,C1,ts,year2,dpCO2a,X);

%%% Nitrogen %%%

% [C1dt,delCdt,delC1] = bioboxone_subN(epsilon,Q1,gamma,ff1,ts,year2,dpCO2a,X);

%[C1dt,C2dt,delCdt,delC1,delC2] = bioboxtwo_subN(epsilon,Q1,Q2,gamma,ff1(601:end,:),ts,year2,dpCO2a,X);
   % calculate land uptake, temperature input doesn't matter because don't use it
 

%[C1dt,C2dt,delCdt,delC1,delC2] = bioboxtwo_Nsize(epsilon,Q1,Q2,gamma,C1,ff1(601:end,:),ts,year2,dpCO2a,X);
   
delCdt(:,2) = -delCdt(:,2);
   
    [delC10] = l_boxcar(delCdt,10,12,1,length(delCdt),1,2);
%     [delC10a] = l_boxcar(delCdt,10,12,1225,length(delCdt),1,2);
%     delC10(1:1284,:) = delCdt(1:1284,:);
%     delC10(1285:length(delC10a),:) = delC10a(1285:end,:);

yhat = delC10(601:end,2);

%yhat = delC10(1081:1513,2);

%yhat = delCdt(1057:1513,2);

%yhat = delC10(601:1777,2);

% figure
% plot(CO2a(2053:3875,1),yhat,CO2a(2053:3875,1),yhat2,CO2a(2053:3875,1),CO2a(2053:3875,2))

% figure
% plot(CO2a(2053:2653,1),CO2a(2053:2653,2))