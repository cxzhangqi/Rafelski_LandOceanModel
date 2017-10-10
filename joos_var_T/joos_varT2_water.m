function [fas]= joos_varT2_water(year,dpCO2a,c,h,kg,T,Aoc,r,dt)

dpCO2s = zeros(length(dpCO2a),2);
dpCO2s(:,1) = dpCO2a(:,1);
dpCO2s_2(:,1) = dpCO2a(:,1);
fas = zeros(length(year),2);
integral = zeros(length(year),length(year));
delDIC = zeros(length(year),2);

DIC0 = 2007.3559;
pCO2s0 = 280;

for m = 1:(length(year)-1)
    %Calculate flux 
   fas(m,1) = year(m);
   fas(m,2) = (kg/Aoc)*(dpCO2a(m,2) - dpCO2s(m,2));
   


    for k = 1:m
        i = find(r(:,1) == year(m-k+1) - year(1));
        integral(k,m) = fas(k,2)*r(i,2)*dt;
    end

    % Calculate delDIC 
    delDIC(m+1,1) = year(m+1);
    delDIC(m+1,2) = (c/h)*sum(integral(1:m,m));
    
    i = find(floor(100*T(:,1)) == floor(100*year(m+1)));
    
    %Calculate dpCO2s
    
%     z0(m+1,1) = 1.7561-0.031618*T(i,2)+0.0004444*(T(i,2))^2;
%     z1(m+1,1) = 0.004096-7.7086E-5*T(i,2)+6.10E-7*(T(i,2))^2;
%     
%     dpCO2s(m+1,2) = z0(m+1,1)*delDIC(m+1,2)/(1-z1(m+1,1)*delDIC(m+1,2));
    
    DIC(m+1,1) = delDIC(m+1,2) + DIC0;
    
    [pCO2s(m+1,2)] = pCO2_calc_water(DIC(m+1,1),T(i,2));
    
    dpCO2s(m+1,2) = pCO2s(m+1,2) - pCO2s0;
    
    dpCO2s_2(m+1,2) = (1.5568 - (1.3993E-2)*T(i,2))*delDIC(m+1,2) + (7.4706-0.20207*T(i,2))*10^(-3)*...
        (delDIC(m+1,2))^2 - (1.2748-0.12015*T(i,2))*10^(-5)*(delDIC(m+1,2))^3 + (2.4491-0.12639*T(i,2))...
        *10^(-7)*(delDIC(m+1,2))^4 - (1.5468-0.15326*T(i,2))*10^(-10)*(delDIC(m+1,2))^5;
%     
%         pCO2s(m+1,2) = (1.5568 - (1.3993E-2)*T(i,2))*(delDIC(m+1,2)+DIC0) + (7.4706-0.20207*T(i,2))*10^(-3)*...
%         (delDIC(m+1,2)+DIC0)^2 - (1.2748-0.12015*T(i,2))*10^(-5)*(delDIC(m+1,2)+DIC0)^3 + (2.4491-0.12639*T(i,2))...
%         *10^(-7)*(delDIC(m+1,2)+DIC0)^4 - (1.5468-0.15326*T(i,2))*10^(-10)*(delDIC(m+1,2)+DIC0)^5;
%     
%     dpCO2s_2(m+1,2) = pCO2s(m+1,2) - pCO2s0;
%     
end

figure
plot(dpCO2s(:,1),dpCO2s(:,2),dpCO2s_2(:,1),dpCO2s_2(:,2))
title('dpco2')

zero1 = find(delDIC(:,1) == 0);
delDIC(zero1,2) = NaN;
delDIC(zero1,1) = NaN;

zero2 = find(dpCO2s(:,2) == 0);
dpCO2s(zero2,2)  = NaN;

zero3 = find(fas(:,1) == 0);
fas(zero3,1:2) = NaN;