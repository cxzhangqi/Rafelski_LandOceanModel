% 3/23/09: Change input file to BP_extrap_CDIAC_data_2007.xls
% 1/10/11: Change input file to BP_extrap_CDIAC_data_2009.xls

function [ff1] = load_fossil(ts)

%% load fossil fuel data

%Read in the data to be interpolated

data1=xlsread('BP_extrap_CDIAC_data_2009.xls');
yr0 = data1(:,1);
fos0 = data1(:,2);

%Get rid of any NaN values
fos0(isnan(yr0)) = [];
yr0(isnan(yr0)) = [];

%Computed cumulative flux
cf = 0;
fos1 = fos0;
yr1 = yr0;
for I = 1:length(yr0),
    cf = cf + fos0(I);
    fos1(I) = cf;
    yr1(I) = yr0(I)+1;  %add 1 because integral is valid at the end of the calendar year
end

%Create new time array
yr2 = [yr0(1):(1/ts):yr0(length(yr0))+1];

%Do interpolation
fos2 = interp1(yr1,fos1,yr2,'spline');

%Computed fluxes by taking discrete time derivative of integrated fluxes
fos = fos2;  %creates dummy arrays
yr = yr2;
fos(length(yr2)) = [];
yr(length(yr2)) = [];
for I = 1:(length(yr2)-1),
    fos(I) = ts.*(fos2(I+1)-fos2(I));
    yr(I) = yr2(I);
end

% convert to ppm

fosppm = fos/2.12;

ff1(:,2) = fosppm(1,:);
ff1(:,1) = yr(1,:);
