% Get relevant land use and Joos model results
%
% Updates
%
% 4/13/07 - Lauren Elmegreen - use land use change 
% emissions that are extrapolated to the present by keeping values 
% constant at the last value, instead of increasing by 1.4%

function [landusemo,ff1,fas,Aoc,extratrop_landmo] = getsourcesink;

load jooshildascale2.mat
load fossil.mat

% interpolate land use changes

month = 1850:(1/12):2006;
landmonth = interp1(landnowppm(:,1),landnowppm(:,12),month);

month2 = 1850:(1/12):2000;
extralandmonth = interp1(extratrop_landppm(:,1),extratrop_landppm(:,2),month2);

landusemo(:,1) = month;
landusemo(:,2) = landmonth; % value in ppm

extratrop_landmo(:,1) = month2;
extratrop_landmo(:,2) = extralandmonth; % value in ppm

% other variables should just be loaded and passed