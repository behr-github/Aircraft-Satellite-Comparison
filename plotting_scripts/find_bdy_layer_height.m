function [ bl_height ] = find_bdy_layer_height( vals, altitude, varargin )
%find_bdy_layer_height: Finds the height of the boundary layer based on criteria specific to the data field in question.
%   Returns: boundary layer height in same units as "altitude"
%
%   This function calculates the boundary layer height.  For most
%   measurments, this will simply be the largest difference between two
%   bins, this works well for well-mixed boundary layers.
%
%
%   For exponential-shaped boundary layers, pass 'exp', 'expo', or
%   'exponential' as the optional third argument. This will then find the
%   boundary layer height as the altitude where the data value has
%   decreased to 1/e of its maximum value.
%
%   For potential temperature, pass 'theta' as the optional third argument.
%   This function then looks for the first instance where the difference between
%   two bins is sufficiently positive and a certain number of the following
%   differences are also sufficiently positive, though less so than the
%   first.
%
%   To plot the vertical profile and mark the boundary layer, pass the
%   parameter "debug" a value of 1

p = inputParser;
p.addRequired('vals',@isnumeric);
p.addRequired('altitude', @isnumeric);
p.addOptional('type','',@isstr);
p.addParamValue('debug',0,@isscalar);

p.parse(vals, altitude, varargin{:});
pout = p.Results;
vals = pout.vals;
altitude = pout.altitude;
find_type = pout.type;
DEBUG = pout.debug;

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
    
% A special method for potential temperature intended for use when a user
% has manually identified a range in which the d(theta)/dz > 0 transition
% occurs
elseif strcmpi(find_type,'thetalite'); 
    delta = diff(vals)./(diff(altitude));
    xx = find(delta > 1);
    bl_height = -1;
    for a=1:numel(xx)
        ind_low = xx(a)+1;
        ind_high = min(xx(a)+4,numel(delta));
        if all(delta(ind_low:ind_high) > 0);
            bl_height = altitude(xx(a));
            break
        end
    end
    
elseif any(strcmpi(find_type,{'exp','expo','exponential'}));
    % Sort and average same x values so that altitude is a monotonically increasing vector
    [altitude, vals] = average_same_x(altitude, vals);
    if ~isrow(altitude); altitude = altitude'; end
    if ~isrow(vals); vals = vals'; end
    efold = exp(-1) * max(vals(:));
    interp_alt = min(altitude):0.01:max(altitude);
    interp_vals = interp1(altitude,vals,interp_alt);
    
    % Find all concentrations that are within 10% of the efolding value and
    % have a negative change in NO2 with altitude
    dx_dz = [diff(interp_vals)./diff(interp_alt), 0];
    efold_logical = abs(interp_vals - efold)/efold < 0.1 & dx_dz < 0;
    bl_height = median(interp_alt(efold_logical));
else
    binwidth = nanmean(diff(altitude));
    dx_dz = diff(vals)./diff(altitude);
    mag_dval = abs(diff(vals));
    if ~iscolumn(altitude); altitude = altitude'; end
    if ~iscolumn(mag_dval); mag_dval = mag_dval'; end
    if ~iscolumn(dx_dz); dx_dz = dx_dz'; end
    S = sortrows([[mag_dval; 0], [(dx_dz<0);0], altitude]); S = flipud(S);
    xx = find(S(:,2),1,'first'); % Find the largest differnce in the input value that has a negative slope w.r.t. altitude
    %xx = find(abs(diff(vals)) == max(abs(diff(vals(2:end)))) & dx_dz < 0);
    bl_height = S(xx,3) + 0.5*binwidth;
end

if DEBUG;
    nfig = nextfig; figure(nfig);
    plot(vals, altitude);
    line([0, max(vals)],[bl_height, bl_height],'color','k','linestyle','--','linewidth',2)
end
end

