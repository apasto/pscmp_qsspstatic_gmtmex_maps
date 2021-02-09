function [lonRange, latRange, lonStep, latStep] = PSCMPsnapshotCheckGrid(snapshotAsTable, varargin)
%PSCMPsnapshotCheckGrid given a 'snapshot' table, read with 'PSCMPsnapshot2table'
%                  from a 'snapshot' file obtained with PSCMP
%                  assert: is it a regular grid? Then obtain the extents
%                  and sampling intervals, arguments required by gmt xyz2grd
%
%   Input arguments:
%      - snapshotAsTable : table, obtained with snapshot2table by reading snapshot file
%      - (optional tolStep) : tolerance used in extracting unique
%                             values, when asserting for constant sampling
%                             'tol' argument of 'uniquetol'
%                             it is an 'epsilon' in a comparison between floats
%                             default: 1e-4 (empirically found)
%                             too small: results in false positives
%                             too large: may miss uneven sampling
%
%   Output arguments:
%      - lonRange : min and max longitude, length 2 vector
%      - latRange : min and max latitude, length 2 vector
%      - lonStep  : longitude sampling interval, scalar
%      - latStep  : latitude sampling interval, scalar
%
%
%   Note: these tests do not assert univocally for a regular grid!
%   They do just enough to avoid patently wrong files
%   (e.g. profiles, sparse points), for any other check we can assert if
%   there are NaNs after xyz2grd, or look for xyz2grd warnings.
%
% 2021-01-07 AP

narginchk(1,2)
nargoutchk(4,4)

assert(istable(snapshotAsTable),...
    'input provided to snapshotCheckGrid is not a table')

% uniquetol tolerance, optional argument
if nargin == 2
    tolStep = varargin{1};
else
    tolStep = 1e-4;
end

% common preamble of error messages, if assertions fail (followed by: reason)
message_notGrid = 'data is not a regular lon by lat grid: ';

% get the extents of data (whatever its form)
lonRange = [min(snapshotAsTable.Londeg), max(snapshotAsTable.Londeg)];
latRange = [min(snapshotAsTable.Latdeg), max(snapshotAsTable.Latdeg)];

% assert that we are not dealing with constant-longitude data (min==max)
assert(lonRange(1)~=lonRange(2),...
    [message_notGrid, 'all points have the same latitude']);

% assert that we are not dealing with constant-latitude data (min==max)
assert(latRange(1)~=latRange(2),...
    [message_notGrid, 'all points have the same longitude']);

% finding the last row with lon = lon0, extract the 'vector of latitudes'
% (an 'axis vector')
% by doing so we get: the number of samples along lat and the lat sampling interval

% index of first lon 'jump' (lon0 to lon1), minus 1 == length of lat vector
latLength = find(snapshotAsTable.Londeg > snapshotAsTable.Londeg(1),1,'first') - 1;

% lat vector
latV = snapshotAsTable.Latdeg(1:latLength);

% assert that the observed range in latV is the same of latRange
assert(all([min(latV), max(latV)] == latRange),...
    [message_notGrid, 'an axis vector of latitude values was extracted, '...
    'but its range is different from the range of latitudes in the whole file'])

% assert that sampling interval is positive and constant
% since we are dealing with floats, we must use 'uniquetol' (with a tolerance)
% otherwise we end up with seemingly identical, non-unique numerical values
latStep = uniquetol(diff(latV), tolStep);
assert(all(latStep > 0),... % using 'all' since we are asserting for constant step in the next test
    [message_notGrid, 'negative sampling interval along latitude detected']);
if numel(latStep) ~= 1
    disp('ERROR: multiple sampling interval along latitude found:')
    if numel(latStep) < 10
        disp(latStep)
    else
        disp('(the first 10 are printed)')
        disp(latStep)
    end
    disp('Note: if values seem apparently identical,')
    disp('      a larger uniquetol tolerance should be tried')
    disp('      (''tolStep'' argument:, an epsilon used in float comparison)')
    disp('      tolStep = ', num2str(tolStep));
    error([message_notGrid, 'sampling interval along latitude is not constant'])
end

% extrapolate the number of along longitude samples
% knowing the number of along latitude samples (in latV)
% and the total number of elements
lonLength = size(snapshotAsTable,1) / size(latV,1);
% extract the along-longitude sampling interval
lonStep = (lonRange(2) - lonRange(1)) / (lonLength - 1);

% the conditions we just asserted seem enough to provide the required
% arguments to gmt xyz2grd: -R range and -I interval
% note: they are not enough to determine univocally that the provided data
%       is a regular grid. A test that we omitted:
%       'find all the unique lon and lat values: are latV and lonV
%       equal to those unique values?' (note that we do not build lonV,
%       while we could)
% if despite these checks the coverage is uneven, xyz2grd will leave NaNs

end

