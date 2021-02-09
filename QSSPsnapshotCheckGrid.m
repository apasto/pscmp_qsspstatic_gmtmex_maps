function [lonRange, latRange, lonStep, latStep] = QSSPsnapshotCheckGrid(snapshotAsTable, varargin)
%QSSPsnapshotCheckGrid given a 'snapshot' table, read with 'QSSPsnapshot2table'
%                  from a 'snapshot' file obtained with QSSP
%                  assert: is it a regular grid? Then obtain the extents
%                  and sampling intervals, arguments required by gmt xyz2grd
%                  Note: while PSCMP can produce regular grid, if asked
%                        QSSP needs them to be externally defined.
%                        Therefore we cannot make the same assumptions!
%                        
%
%   Input arguments:
%      - snapshotAsTable : table, obtained with snapshot2table by reading snapshot file
%      - (optional tolStep) : tolerance used in extracting unique values
%                             'tol' argument of 'uniquetol'
%                             it is an 'epsilon' in a comparison between floats
%                             default: 1e-12 (double precision)
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
% 2021-01-23 AP

narginchk(1,2)
nargoutchk(4,4)

assert(istable(snapshotAsTable),...
    'input provided to snapshotCheckGrid is not a table')

% uniquetol tolerance, optional argument
if nargin == 2
    tolStep = varargin{1};
else
    tolStep = 1e-8;
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

% extract the sampling step along longitude and latitude
% by finding the "smallest positive non-zero interval" between subsequent points
uniqueLonSteps = uniquetol(diff(snapshotAsTable.Londeg), tolStep);
uniqueLatSteps = uniquetol(diff(snapshotAsTable.Latdeg), tolStep);
% get rid of negatives and zero (tolerance included)
uniqueLonSteps(uniqueLonSteps<=(0 + tolStep)) = [];
uniqueLatSteps(uniqueLatSteps<=(0 + tolStep)) = [];
% usually we should get either: one negative value (jump when cycle of lon or lats
% begins) and interval, or zero and interval
% depends on how points are ordered: along lat-first or lon-first
lonStep = min(uniqueLonSteps);
latStep = min(uniqueLatSteps);

end

