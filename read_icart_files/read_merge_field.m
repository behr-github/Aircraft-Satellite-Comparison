function [ varargout ] = read_merge_field( Merge, fieldname )
%read_merge_field(Merge, fieldname): Read in a merge field along with UTC, altitude, and lat/lon data, replacing fills with nans

ulod = Merge.metadata.upper_lod_flag;
llod = Merge.metadata.lower_lod_flag;
fill = Merge.Data.(fieldname).Fill;

yy = Merge.Data.(fieldname).Values == fill | Merge.Data.(fieldname).Values == ulod | Merge.Data.(fieldname).Values == llod;

data = Merge.Data.(fieldname).Values; 
data(yy) = NaN;
varargout{1} = data;

if nargout >= 2
    utc = Merge.Data.UTC.Values; 
    utc(yy) = NaN;
    varargout{2} = utc;
end
if nargout >= 3
    alt = Merge.Data.ALTP.Values; 
    alt(yy) = NaN;
    varargout{3} = alt;
end
if nargout >= 4
    lon = Merge.Data.LONGITUDE.Values - 360; 
    lon(yy) = NaN;
    varargout{4} = lon;
end
if nargout >= 5
    lat = Merge.Data.LATITUDE.Values; 
    lat(yy) = NaN;
    varargout{5} = lat;
end

end

