function gridsInStruct = snapshotTable2gmtMexGrid(...
    snapshotAsTable, lonStep, latStep, lonRange, latRange, varargin)
%snapshotTable2gmtMexGrid using gmt xyz2grd, convert a loaded PSCMP
%                         snapshot to a structure containing
%                         a gmt mex grid structure for each variable
%                         in the snapshot.
%
%   Input arguments:
%      - snapshotAsTable : table, obtained with snapshot2table by reading snapshot file
%      (the following arguments are obtained with snapshotCheckGrid)
%      - lonRange : min and max longitude, length 2 vector
%      - latRange : min and max latitude, length 2 vector
%      - lonStep  : longitude sampling interval, scalar
%      - latStep  : latitude sampling interval, scalar
%      - optional, verboseFlag : print progress and provide '-V' option to xyz2grd
%                                default: false (no progress printout and use gmt defaults)
%
%   Output arguments:
%      - gridsInStruct : scalar struct, each field a gmt mex grid structure
%
% 2021-01-13 AP

narginchk(5,6)
nargoutchk(0,1)

if nargin==6
    verboseFlag = varargin{1};
else
    verboseFlag = false;
end

if verboseFlag
    verboseOption = ' -V';
else
    verboseOption = '';
end
    

assert(istable(snapshotAsTable),...
    'input provided to snapshotTable2gmtMexGrid is not a table')

% get variable (columns) in snapshot
VariableNames = snapshotAsTable.Properties.VariableNames;
% assert that there are two columns namend 'Latdeg' and 'Londeg'
% (otherwise the wrong table was provided or something changed)
assert(any(strcmp('Latdeg', VariableNames)),...
    'table provided to snapshotTable2gmtMexGrid has no ''Latdeg'' column')
assert(any(strcmp('Londeg', VariableNames)),...
    'table provided to snapshotTable2gmtMexGrid has no ''Londeg'' column')
% remove them from VariableNames
VariableNames(strcmp('Londeg', VariableNames)) = [];
VariableNames(strcmp('Latdeg', VariableNames)) = [];

% perform conversion with gmt xyz2grd
for n=1:size(VariableNames, 2)
    if verboseFlag
    fprintf(['\nxyz2grd: called on ', VariableNames{n}, '\n'])
    end
    gridsInStruct.(VariableNames{n}) = gmt(...
        ['xyz2grd', verboseOption,...
        ' -I', num2str(lonStep), '/', num2str(latStep),...
        ' -R', num2str(lonRange(1)), '/', num2str(lonRange(2)), '/',...
        num2str(latRange(1)), '/', num2str(latRange(2))],...
        horzcat(...
            snapshotAsTable.Londeg,...
            snapshotAsTable.Latdeg,...
            snapshotAsTable.(VariableNames{n})));
    if verboseFlag
        fprintf(['xyz2grd: done with ', VariableNames{n}, '\n'])
    end
end

end

