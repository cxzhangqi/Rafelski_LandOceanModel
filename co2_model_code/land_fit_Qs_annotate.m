% Things to check: end year defined correctly? line 12
% beta variables defined correctly? lines 22-32
% correct model chosen? lines 36-43
% correct time period for yhat? line 55

function yhat = land_fit_Qs_annotate(beta,X)

% Initial conditions - set according to Joos et al

ts = 12;
start_year = 1850;
end_year = 2009+(7/12);

% Get parameters

[dtdelpCO2a,dpCO2a,year,dt,CO2a] = MLOinterpolate_increment2(ts,start_year,end_year);

[landusemo,ff1,fas,Aoc,extratrop_landmo] = getsourcesink_scale3;

%% To make temperature-independent: set Q1 and Q2 to 1

%% For CO2 fertilization model
epsilon = beta(1);
Q1 = beta(2);
Q2 = 1;%beta(2);

%% For N fertilization model. Set epsilon = 0 to not allow CO2
%% fertilization
% epsilon = 0;
% gamma = beta(1);
% Q1 = beta(2);
% Q2 = 1;%beta(3);

year2 = year';

%% CO2 fertilization

[C1dt,C2dt,delCdt,delC1,delC2] = bioboxtwo_sub10_annotate(epsilon,Q1,Q2,ts,year2,dpCO2a,X); 


%% Nitrogen fertilization

%[C1dt,C2dt,delCdt,delC1,delC2] = bioboxtwo_subN(epsilon,Q1,Q2,gamma,ff1(601:end,:),ts,year2,dpCO2a,X);
    
delCdt(:,2) = -delCdt(:,2);

%% 10-year moving boxcar average of model
   [delC10] = l_boxcar(delCdt,10,12,1,length(delCdt),1,2);


%% yhat is the term that is compared to the residual flux in nlinfit. 
%% Change the index numbers here and in nonlin_land_Qs_annotate (e.g. line
%% 158) to fit to a different time period

yhat = delC10(601:end,2);
