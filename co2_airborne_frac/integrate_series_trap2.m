function [func_int] = integrate_series_trap2(func,timecol,numcol,dt)

% Integrate data based on trapezoid rule. Time points and data points
% should be in the same array, with time points in the column "timecol" and
% data points in the column "numcol"

% dt = number of points per year

func_int = zeros(length(func),numcol);
func_int(:,1) = func(:,timecol);

for i = 2:length(func(:,timecol))
    func_int(i,2) = func_int(i-1,2) + 0.5*(func(i-1,numcol)+func(i,numcol))/dt;
end