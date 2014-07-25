function varargout = quality_flag_filter(bad_flags, qflag_in, varargin)

% quality_flag_filter - filters an arbitrary number of data sets based on
% their bitarray quality flags.
%
%   This function takes as arguments 1) a number representing the bits that
%   should NOT be set, 2) a list of the quality flags that correspond the
%   the data in, and 3+) each of the data sets that you want filtered.  
%
%   This function works by comparing all quality flags associated with a
%   particular measurement against the bad_flags input with bitand.  If any
%   of the flags match, then that measurement is rejected.  As an example,
%   consider data with a 4-bit flag.  If bad_flags is set to 2 (i.e. 0010),
%   then any measurement with a quality flag where the second bit (**1*) is
%   1 will be removed.
%
%   This returns each data set with offending entries removed and the
%   quality flags as the final argument.

% At least one data set must be passed in or this is a silly function
if numel(varargin) < 1;
    error('quality_flag_filter:numargs','Must pass at least one data set');
end

% Copy varargin 
data = varargin;

% Check the qflags passed in; if not in a cell array, make them so
if ~iscell(qflag_in);
    qflags = mat2cell(qflag_in);
else
    qflags = qflag_in;
end

% Prep a logical index variable that will keep track of which measurements
% to remove
xx = true(size(qflags));

for a=1:numel(qflags)
    if any(bitand(bad_flags,qflags{a}))
        xx(a) = false;
    end
end

qflags_out = qflags(xx);
data_out = cell(size(data));
for v=1:numel(data);
    data_out{v} = data{v}(xx);
end

varargout = data_out;
varargout{end+1} = qflags_out;

end
