% New subroutine friendly version of jooshilda scaled

clear all

ts = 12;
start_year = 1800;
end_year = 2008;%2008-(1/12);
Aoc = 3.62E14;
c = 1.722E17;
h = 75;
T = 18.2;
kg = 1/9.06;

[dtdelpCO2a,dpCO2a,year,dt] = MLOinterpolate_increment(ts,start_year,end_year);

[ff1] = load_fossil(ts);


% Response function
for n = 1:length(year)
     t(n,1) = year(n) - year(1);
     r(n,1) = t(n,1);
    %Calculate response function based on HILDA equation
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

%[fas] = joos_general2(year,dpCO2a,c,h,kg,T,Aoc,r,dt);
[fas,dpCO2s] = joos_general_fast(year(1:2496),dpCO2a,c,h,kg,T,Aoc,r,dt); 
% year(1:2496) is because have to go to 2008 to get month values to agree, 
%but only have data until 2007+(11/12)

% Calculate land flux using fossil fuel sources and ocean sink (ocean
% flux*area of ocean)


for p = 1:(length(year)-7);%1:(length(year)-1)
    q = find(ff1(:,1) == year(1,p));
    landflux(p,1) = year(p);
    landflux(p,2) = dtdelpCO2a(p+((1800-1640)*ts),2) - ff1(q,2) + fas(p,2)*Aoc; % changed 1679 to 1640: new CO2 record goes farther back
end


x = 0*year;

figure
plot(ff1(:,1),ff1(:,2),'-k',dtdelpCO2a(:,1),dtdelpCO2a(:,2),'-r',fas(:,1),-Aoc*fas(:,2),'-b',landflux(:,1),landflux(:,2),'-g',year(1,:),x,'--k')
axis([1800 2010 -10 10])
legend('fossil fuel','atmosphere','ocean','land','Location','SouthWest')
title('Sources and sinks from Joos 10 yr response function ')
xlabel('Year ')
ylabel('ppm/year  Positive = source, negative = sink ')