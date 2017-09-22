%% 3/13/08: changed CO2 dataset to spline fit with sigma =0.6 to capture
%% 1940s plateau
%% 4/2/08: changed CO2 dataset to updated ice core data, sigma=0.6
%% 4/24/08: add in an option to decrease ice core data by 2 ppm; change to
%% linear interpolation for this option

function [annincMLOSPO,dpCO2a,year,dt,MLOSPOiceinterp] = MLOinterpolate_increment(ts,start_year,end_year)

%load MLOSPO.mat

%load MLOSPO_6.mat

% load co2_2008.mat
%load co2_2008_2.mat
load co2_2008_3.mat

dt = 1/ts;

%% Interpolate ice core CO2 data
% years0 = MLOSPOice(:,1);
% 
% CO2 = MLOSPOice(:,2);
% 
% years0 = co2sig6(:,1);
% 
% CO2 = co2sig6(:,2);

% years0 = co2_2008(:,1);
% 
% CO2 = co2_2008(:,2);

% Decrease ice core data by 2 ppm; cut out ice core data after 1955
% 
% lowice(1:3780,1) = mlospo_meure(1:3780,1);
% lowice(1:3780,2) = mlospo_meure(1:3780,2) - 2;
% 
% mlospo_meure_2(1:3780,:) = lowice(:,:);
% mlospo_meure_2(3781:4378,:) = mlospo_meure(3807:end,:);


%% This one:
years0 = mlospo_meure(:,1);

CO2 = mlospo_meure(:,2);

%Create new time array
years = [years0(1):1/ts:years0(length(years0))];

%Do interpolation
%CO22 = interp1(years0,CO2,years,'linear');
CO22 = interp1(years0,CO2,years,'spline');

MLOSPOiceinterp(:,1) = years;
MLOSPOiceinterp(:,2) = CO22;

% CO2 increment with monthly resolution, in ppm/year
% n = 7 is 7/1958, last value is 7/2005

for n = ((ts/2)+1):(length(years)-(ts/2))
    annincMLOSPO(n,1) = MLOSPOiceinterp(n,1);
    annincMLOSPO(n,2) = MLOSPOiceinterp(n+(ts/2),2) - MLOSPOiceinterp(n-(ts/2),2);
end

year = start_year:dt:end_year;

%i = find(MLOSPOiceinterp(:,1) == start_year);
i = find(floor(100*MLOSPOiceinterp(:,1)) == floor(100*(start_year+(1/24))));
dpCO2a(:,1) = MLOSPOiceinterp(i:length(MLOSPOiceinterp),1);
dpCO2a(:,2) = MLOSPOiceinterp(i:length(MLOSPOiceinterp),2)-MLOSPOiceinterp(i,2);
