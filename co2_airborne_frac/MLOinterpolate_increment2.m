%% 3/13/08: changed CO2 dataset to spline fit with sigma =0.6 to capture
%% 1940s plateau
%% 4/2/08: changed CO2 dataset to updated ice core data, sigma=0.6
%% 4/24/08: add in an option to decrease ice core data by 2 ppm; change to
%% linear interpolation for this option
%% 1/10/11: add updated CO2 dataset through 2010

function [annincMLOSPO,dpCO2a,year,dt,MLOSPOiceinterp] = MLOinterpolate_increment2(ts,start_year,end_year)

load co2_2011_2.mat

dt = 1/ts; 

years0 = mlospo_meure(:,1);

CO2 = mlospo_meure(:,2);

%% Create new time array
years = [years0(1):1/ts:years0(length(years0))];

%% Do interpolation
CO22 = interp1(years0,CO2,years,'spline');

MLOSPOiceinterp(:,1) = years;
MLOSPOiceinterp(:,2) = CO22;

%% CO2 increment with monthly resolution, in ppm/year
% n = 7 is 7/1958, last value is 7/2005

for n = ((ts/2)+1):(length(years)-(ts/2))
    annincMLOSPO(n,1) = MLOSPOiceinterp(n,1);
    annincMLOSPO(n,2) = MLOSPOiceinterp(n+(ts/2),2) - MLOSPOiceinterp(n-(ts/2),2);
end

year = start_year:dt:end_year;

%% Calculate change in atmospheric concentration
i = find(floor(100*MLOSPOiceinterp(:,1)) == floor(100*(start_year+(1/24))));
dpCO2a(:,1) = MLOSPOiceinterp(i:length(MLOSPOiceinterp),1); 
dpCO2a(:,2) = MLOSPOiceinterp(i:length(MLOSPOiceinterp),2)-MLOSPOiceinterp(i,2);
