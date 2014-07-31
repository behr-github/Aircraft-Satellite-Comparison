function [ avg_slope, avg_R2 ] = bootstrapping_analysis( x, y, percent_of_data, number_of_trials, varargin )
%bootstrapping_analysis Randomly selects a portion of the data and does a
%fit, repeated n times.
%   Detailed explanation goes here

if isrow(x); x = x'; end
if isrow(y); y = y'; end
if numel(x) ~= numel(y); error('bootstrap:xy_unequal','x and y must have equal lengths'); end
if percent_of_data > 1; error('bootstrap:percent_of_data','Percent of data must be a number between 0 and 1'); end

nans = isnan(x) | isnan(y);
x(nans) = []; y(nans) = [];
M = [x,y];
count = ceil(percent_of_data*numel(x));

slopes = zeros(number_of_trials,1);
Rsquareds = zeros(number_of_trials,1);

for a=1:number_of_trials
    x_set = zeros(count,1); y_set = zeros(count,1);
    M_set = M;
    for b=1:count
        ri = randi([1, size(M_set,1)],1);
        x_set(b) = M_set(ri,1); y_set(b) = M_set(ri,2);
        M_set(b,:) = [];
    end
    
    [m,~,r] = lsqfitma(x_set,y_set);
    slopes(a) = m;
    Rsquareds(a) = r^2;
end

avg_slope = mean(slopes);
avg_R2 = mean(Rsquareds);
end

