function [ bl_height ] = find_bdy_layer_height( vals, altitude, varargin )
%find_bdy_layer_height: Finds the height of the boundary layer based on criteria specific to the data field in question.
%   This function calculates the boundary layer height.  For most
%   measurments

p = inputParser;
p.addRequired('vals',@isnumeric);
p.addRequired('altitude', @isnumeric);
p.addOptional('type','',@isstr);

p.parse(vals, altitude, varargin{:});
pout = p.Results;
vals = pout.vals;
altitude = pout.altitude;
find_type = pout.type;

if ~isempty(regexpi(find_type,'theta'))
    delta = diff(vals)./(diff(altitude));
    xx = find(delta > 3);
    bl_height = -1;
    for a=1:numel(xx)
        ind_low = xx(a)+1;
        %end_count = 1 / median(diff(altitude));
        %ind_high = min(xx(a)+end_count,numel(delta));
        ind_high = min(xx(a)+4,numel(delta));
        if all(delta(ind_low:ind_high) > 1);
            bl_height = altitude(xx(a));
            break
        end
    end
else
    xx = find(abs(diff(vals)) == max(abs(diff(vals(2:end)))));
    bl_height = mean(altitude(xx:xx+1));
end
end

