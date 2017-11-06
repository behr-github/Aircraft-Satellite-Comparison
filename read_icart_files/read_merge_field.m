function [ varargout ] = read_merge_field( Merge, fieldname, varargin )
%READ_MERGE_FIELD Read a fill from a Merge structure, replacing fills with NaNs
%
%   [ VALUE ] = READ_MERGE_FIELD( MERGE, FIELDNAME ) will read the Value
%   field from MERGE.Data.(FIELDNAME), replace any fill, upper LOD, or
%   lower LOD values with NaNs, and return that as VALUE.
%
%   [ VALUE, UTC, ALT, LON, LAT ] = READ_MERGE_FIELD( MERGE, FIELDNAME )
%   will simultaneously return the UTC time, altitude, longitude, and
%   latitude, assuming that their field names are "UTC", "ALTP",
%   "LONGITUDE", and "LATITUDE".
%
%   [ ___ ] = READ_MERGE_FIELD( ___, 'unit', UNIT ) will use CONVERT_UNITS
%   to convert the data given in FIELDNAME from the units listed in the
%   Merge structure to UNIT.

p = inputParser;
p.addParameter('unit', '');

p.parse(varargin{:});
pout = p.Results;

target_unit = pout.unit;

ulod = Merge.metadata.upper_lod_flag;
llod = Merge.metadata.lower_lod_flag;
fill = Merge.Data.(fieldname).Fill;

yy = Merge.Data.(fieldname).Values == fill | Merge.Data.(fieldname).Values == ulod | Merge.Data.(fieldname).Values == llod;

data = Merge.Data.(fieldname).Values; 
data(yy) = NaN;
if ~isempty(target_unit)
    data = convert_units(data, Merge.Data.(fieldname).Unit, target_unit);
end
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

