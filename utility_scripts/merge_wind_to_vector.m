function [ u, v ] = merge_wind_to_vector( WND, WNS )
%MERGE_WIND_TO_VECTOR( WND, WNS ) Convert ICARTT merge wind direction to u & v vectors
%   As far as I know, wind directions in ICARTT files are given as degrees
%   clockwise from north, meaning that 0 deg indication the wind is blowing
%   from the south to the north, 90 deg means from the west to the east,
%   and so on. This is highly inconvinent for converting to wind vectors,
%   as the usual definition of U and V wind vectors relies on wind blowing
%   to due east being defined as 0, so that the winds are described on a
%   traditional cartesian grid. (Some campaigns may have U and V winds
%   directly included in the merge, in which case this is unnecessary.)
%
%   The transformation between the two coordinate systems is actually just
%   a straightforward y = x transform:
%
%                               | cos(0) = sin(90)
%                               | sin(0) = cos(90)
%                               |                   cos(90) = sin(0)
%                               |                   sin(90) = cos(0)
%   ----------------------------------------------------------
%   cos(270) = sin(180)         |
%   sin(270) = cos(180)         |
%                               | cos(180) = sin(270)
%                               | sin(180) = cos(270)
%
%   The quantities on the left side of the = are in the original (deg. CW
%   from north) coordinate system, the quantities on the right are with the
%   angles in the new (deg CCW from east) coordinate system.
%
%   The first input is the wind direction from the merge files. The second
%   is optional, it it the wind speed from the merge files.  If given, then
%   the u and v output vectors will reflect the magnitude of the wind
%   vector. If not given, the magnitude will be assumed to be 1.
%
%   Josh Laughner <joshlaugh5@gmail.com> 31 Jul 2015


if nargin < 2
    WNS = ones(size(WND));
end

v = cosd(WND).*WNS;
u = sind(WND).*WNS;



end

