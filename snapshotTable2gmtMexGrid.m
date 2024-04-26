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
%         either a table or a 2-cell cell array, with a table in each cell
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
    
table_was_provided = istable(snapshotAsTable);
cell_was_provided = iscell(snapshotAsTable);
if cell_was_provided
    snapshotAsTable = snapshotAsTable(:);
    % assert correct size of cell
    assert(all(size(snapshotAsTable) == [2, 1]), ...
        'cell array provided to snapshotTable2gmtMexGrid is not a 2-cell cell array');
    % assert table in each cell
    assert(istable(snapshotAsTable{1}) && istable(snapshotAsTable{2}), ...
        'cell array provided to snapshotTable2gmtMexGrid must contain  a table in each cell')
end

assert(table_was_provided || cell_was_provided,...
    'input provided to snapshotTable2gmtMexGrid is not a table or a cell')

% get variable (columns) in snapshot
if table_was_provided
    VariableNames = snapshotAsTable.Properties.VariableNames;
else
    VariableNames = snapshotAsTable{1}.Properties.VariableNames;
    VariableNames_b = snapshotAsTable{2}.Properties.VariableNames;
    % TO DO: variable names may be the same, but with different order
    % should assert using set and contains / does not contain
    assert(all(strcmp(VariableNames, VariableNames_b)), ...
        'Different variable names in the two tables provided')
end
% assert that there are two columns namend 'Latdeg' and 'Londeg'
% (otherwise the wrong table was provided or something changed)
assert(any(strcmp('Latdeg', VariableNames)),...
    'table provided to snapshotTable2gmtMexGrid has no ''Latdeg'' column')
assert(any(strcmp('Londeg', VariableNames)),...
    'table provided to snapshotTable2gmtMexGrid has no ''Londeg'' column')
% remove them from VariableNames
VariableNames(strcmp('Londeg', VariableNames)) = [];
VariableNames(strcmp('Latdeg', VariableNames)) = [];

% two tables provided: compute difference FIRST
% second minus one (this is reflected in filenames in 'snapshotGrid2gmtMap.m'
% this avoids precision issues (data in gmt / gmt mex is single precision)
if cell_was_provided
    snapshotData = snapshotAsTable{2};
    for n=1:size(VariableNames, 2)
        snapshotData.(VariableNames{n}) = ...
            snapshotData.(VariableNames{n}) - snapshotAsTable{1}.(VariableNames{n});
    end
else
    snapshotData = snapshotAsTable;
end

% perform conversion with gmt xyz2grd
for n=1:size(VariableNames, 2)
    if verboseFlag
    fprintf(['\nxyz2grd: called on ', VariableNames{n}, '\n'])
    end
    gridsInStruct.(VariableNames{n}) = gmt(...
        ['xyz2grd', verboseOption,...
        ' -fg -I', num2str(lonStep), '/', num2str(latStep),...
        ' -R', num2str(lonRange(1)), '/', num2str(lonRange(2)), '/',...
        num2str(latRange(1)), '/', num2str(latRange(2))],...
        horzcat(...
            snapshotData.Londeg,...
            snapshotData.Latdeg,...
            snapshotData.(VariableNames{n})));
    if verboseFlag
        fprintf(['xyz2grd: done with ', VariableNames{n}, '\n'])
    end
end

end
