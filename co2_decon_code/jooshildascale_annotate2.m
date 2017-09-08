%% January 10, 2011: extended data to 2010

clear all

ts = 12; % number of data points/year
start_year = 1800;
end_year = 2010;
Aoc = 3.62E14; % surface area of ocean, m^2, from Joos 1996
c = 1.722E17; % unit converter, umol m^3 ppm^-1 kg^-1, from Joos 1996
h = 75; % mixed layer depth, m, from Joos 1996
T = 18.2; % surface temperature, deg C, from Joos 1996
kg = 1/9.06; % gas exchange rate, yr^-1, from Joos 1996

[dtdelpCO2a,dpCO2a,year,dt] = MLOinterpolate_increment2(ts,start_year,end_year); % get atmospheric CO2 record

[ff1] = load_fossil2(ts); % get fossil fuel emissions


% Response function to calculate ocean uptake
for n = 1:length(year)
     t(n,1) = year(n) - year(1);
     r(n,1) = t(n,1);
    %Calculate response function based on HILDA equation in Joos 1996
     if t(n,1) == 0
         r(n,2) = 1;
     elseif t(n,1) <= 2
         r(n,2)= (1/0.95873)*(0.12935+0.21898*exp(-t(n,1)/0.034569)+0.17003*exp(-t(n,1)/0.26936)...
             +0.24071*exp(-t(n,1)/0.96083)+0.24093*exp(-t(n,1)/4.9792));
     else
         r(n,2) = (1/0.95873)*(0.022936+0.24278*exp(-t(n,1)/1.2679)+0.13963*exp(-t(n,1)/5.2528)...
                +0.089318*exp(-t(n,1)/18.601)+0.037820*exp(-t(n,1)/68.736)...
                +0.035549*exp(-t(n,1)/232.3));
     end
end

%% calculate ocean uptake
[fas,dpCO2s] = joos_general_fast_annotate2(year,dpCO2a,c,h,kg,T,Aoc,r,dt); 

%% Calculate land flux using fossil fuel sources and ocean sink (ocean flux*area of ocean)
 
for p = 1:(length(year)-7);%1:(length(year)-1)
    q = find(ff1(:,1) == year(1,p));
    landflux(p,1) = year(p);
    landflux(p,2) = dtdelpCO2a(p+((1800-1640)*ts),2) - ff1(q,2) + fas(p,2)*Aoc; 
end


x = 0*year;

figure
plot(ff1(:,1),ff1(:,2),'-k',dtdelpCO2a(:,1),dtdelpCO2a(:,2),'-r',fas(:,1),-Aoc*fas(:,2),'-b',landflux(:,1),landflux(:,2),'-g',year(1,:),x,'--k')
axis([1800 2010 -10 10])
legend('fossil fuel','atmosphere','ocean','land','Location','SouthWest')
title('Sources and sinks from Joos response function ')
xlabel('Year ')
ylabel('ppm/year  Positive = source, negative = sink ')