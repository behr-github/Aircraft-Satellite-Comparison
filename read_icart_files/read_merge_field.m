function [ data, utc, alt, lon, lat ] = read_merge_field( Merge, fieldname )
%read_merge_field(Merge, fieldname): Read in a merge field along with UTC, altitude, and lat/lon data, replacing fills with nans

ulod = Merge.metadata.upper_lod_flag;
llod = Merge.metadata.lower_lod_flag;
fill = eval(sprintf('Merge.Data.%s.Fill',fieldname));

yy = eval(sprintf('Merge.Data.%s.Values',fieldname)) == fill | eval(sprintf('Merge.Data.%s.Values',fieldname)) == ulod | eval(sprintf('Merge.Data.%s.Values',fieldname)) == llod;

data = eval(sprintf('Merge.Data.%s.Values',fieldname)); data(yy) = NaN;
utc = Merge.Data.UTC.Values; utc(yy) = NaN;
alt = Merge.Data.ALTP.Values; alt(yy) = NaN;
lon = Merge.Data.LONGITUDE.Values - 360; lon(yy) = NaN;
lat = Merge.Data.LATITUDE.Values; lat(yy) = NaN;

end

