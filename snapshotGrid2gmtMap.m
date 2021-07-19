function snapshotGrid2gmtMap(dataPath, snapshotFilename, PSCMPQSSPswitch, ExtentsBuffer, varargin)
% snapshotGrid2gmtMap still a work in progress! TO DO: documentation
%    (sub-functions called from here are already documented)
%PSCMPsnapshotCheckGrid given a 'snapshot' table, read with 'PSCMPsnapshot2table'
%                  from a 'snapshot' file obtained with PSCMP
%                  assert: is it a regular grid? Then obtain the extents
%                  and sampling intervals, arguments required by gmt xyz2grd
%
%   Automated parameters:
%    - major and minor z-ticks (to do: documentation)
%    - grid line step for lon and lat: if range in lon or lat (or both)
%                                      is less than 2 deg
%                                      interval is 0.5 deg,
%                                      else: 1 deg
%    - country borders ('gmt pscoast -N') are plotted only if range
%      in lon or lat (or both) is more than 6 deg
%
%   Input arguments:
%      - dataPath : string, path to the directory with the snapshot file
%                   trailing slash included!
%      - snapshotFilename : string or cell array of strings, filename(s) of the snapshot file(s)
%      - PSCMPQSSPswitch : string, switch between PSCMP and QSSP output files
%      - ExtentsBuffer : trim data on edges or buffer around edges
%                           either a:
%                            - 4 element vector [deg] with forced extents
%                                 [lon0, lon1, lat0, lat1]
%                            - scalar float [deg], tighter (negative) or larger (positive)
%                              border around data (e.g. a negative value trims data and
%                              zooms on center of covered extents)
%                            - 2 element vector: [deg], same behaviour as scalar
%                              but is separated in [lon bugger, lat buffer]
%                            If argument 'sourceLonLat' is NOT provided,
%                            the trim/buffer values refer to the extents
%                            of the provided grid.
%                            If 'sourceLonLat' is provided, they refer to
%                            the provided center point or to the
%                            average coordinates of the 4 fault patch vertices.
%
%   Optional input arguments:
%      - (optional gsPath) : string, path to gswin64c or gswin32c
%                            (including gswin64c/gswin32c exe file)
%                            required if it cannot be found in your PATH
%                            (to test it: call either gswin64c or gswin32c
%                            from a terminal prompt: does it succeed?)
%                            default: 'gswin64c'
%      - (optional doContours) : plot contours? (default: true)
%      - (optional cptQuantiles) : quantiles for colormap min/max, 2 element vector
%                                  (default: 0.001, 0.999)
%      - (optional sourceLonLat) : plot the epicenter position, 2 element vector
%                                  (default: is not plotted)
%      - (optional noPlots) : only grd conversion, skip plotting
%                             (default: false)
%      - (optional extendedTitle) : filename in map title
%                                   (default: false)
%
%   Output arguments:
%      none
%
% 2021-01-24 AP

narginchk(4,10)
nargoutchk(0,0)

% gspath optional argument
gsPath = 'gswin64c'; % default: assume it can be found in PATH
if nargin>4 && ~isempty(varargin{1})
    assert((ischar(varargin{1}) || isstring(varargin{2})),...
        'gsPath argument must be type char or string.');
    gsPath = varargin{1};
end

% doContours optional argument
doContours = true; % default
if nargin>5 && ~isempty(varargin{2})
    doContours = logical(varargin{2});
end

% cpt quantiles optional argument
MinQuantile = 0.001; % default
MaxQuantile = 0.999; % default
if nargin>6 && ~isempty(varargin{3})
    assert(isequal(size(varargin{3}(:)), [2, 1]),...
        'CPT quantiles argument must be a 2 element vector');
    MinQuantile = varargin{3}(1);
    MaxQuantile = varargin{3}(2);
end

% source point optional argument
SourceType = 'none';
sourceLonLat = [];
if nargin>7 && ~isempty(varargin{4})
    % TO DO: this could allow for drawing multiple rectangles or multiple
    %        point sources, e.g. by providing a cell array.
    %        Also: may provide focal mechanism to relevant gmt module.
    sourceLonLat = varargin{4};
    if isequal(size(sourceLonLat), [1, 2]) % one point [lon, lat]
        SourceType = 'point';
    elseif isequal(size(sourceLonLat), [4, 2]) % 4-point rectangle
        SourceType = 'rectangle';
    else
        error('Source Lon Lat argument must be a 2 element vector (point) or a 4 by 2 array (rectangle).')
    end
end

% optional argument: only grd conversion, no plotting
if nargin>8 && ~isempty(varargin{5})
    noPlotsFlag = logical(varargin{5});
else
    noPlotsFlag = false;
end

% optional argument: extended titles with filename
if nargin>9 && ~isempty(varargin{6})
    extendedTitle = logical(varargin{6});
else
    extendedTitle = false;
end

% delete leftovers from gmtset
if exist('./gmt.conf', 'file')
    delete('./gmt.conf') 
end

% PSCMPQSSPswitch = 'QSSP'; % 'PSCMP' / 'QSSP', case insensitive
assert(ischar(PSCMPQSSPswitch) || isstring(PSCMPQSSPswitch),...
    'Switch argument between PSCMP and QSSP must be type char or string.');
assert(any(strcmpi(PSCMPQSSPswitch, {'PSCMP', 'QSSP'})),...
    ['Switch between PSCMP and QSSP only. Provided argument: ', PSCMPQSSPswitch])

% manage edge buffer/trim or forced extents
if isscalar(ExtentsBuffer) % only edge buffer/trim
    ForceExtents = false;
    edgeBufferLon = ExtentsBuffer;
    edgeBufferLat = ExtentsBuffer;
elseif isequal(size(ExtentsBuffer(:)), [2 ,1]) % [lon buffer, lat buffer]
    ForceExtents = false;
    edgeBufferLon = ExtentsBuffer(1);
    edgeBufferLat = ExtentsBuffer(2);
elseif isequal(size(ExtentsBuffer(:)), [4 ,1]) % [lon0, lon1, lat0, lat1] forced extents
    ForceExtents = true;
    % will use ExtentsBuffer as extents
else
    error('Size of ExtentsBuffer argument is not valid.')
end

% static parameters for plotting fine/coarse grid line interval
GridLines_FineIntervalExtentsThreshold = 2.25; % [deg], if smaller go from 'coarse' to 'fine' interval
GridLines_FineInterval = 0.5; % [deg]
GridLines_CoarseInterval = 1.0; % [deg]

% static parameters for plotting country borders or not
PlotCountries_ExtentsThreshold = 6.0; % deg

%% allow for single file or couple of files (difference, 2nd minus 1st)
if ischar(snapshotFilename) || isstring(snapshotFilename)
    % one single file provided (as char array or string)
    IsDifference = false; % => we are plotting a single-snapshot map
    % convert to cell, to use brace indexing
    snapshotFilename_temp = snapshotFilename;
    snapshotFilename = cell(1,1);
    snapshotFilename{1} = snapshotFilename_temp;
    clear snapshotFilename_temp
elseif iscell(snapshotFilename) && isequal(size(snapshotFilename(:)), [2, 1])
    % two files provided, as a (2, 1)-sized cell-array
    assert(... % assert that the elements are type char or string
        or(ischar(snapshotFilename{1}), isstring(snapshotFilename{1})) && ...
        or(ischar(snapshotFilename{2}), isstring(snapshotFilename{2})),...
        'Contents of snapshot filenames cell array must be type char or string.')
    snapshotFilename = snapshotFilename(:); % ensure consistent size: 2x1
    IsDifference = true; % => we are plotting a snapshot-couple difference map
else % cases not accounted for: throw error
    error('The snapshot filename type is not allowed.')
end

% if a .mat file is provided, set a flag to skip import function
if ~IsDifference
    if strcmpi(snapshotFilename{1}(end-3:end),'.mat')
        isMATfile = true;
    else
        isMATfile = false;
    end
else
    if strcmpi(snapshotFilename{1}(end-3:end),'.mat')
        assert(strcmpi(snapshotFilename{2}(end-3:end),'.mat'),...
            'If the first filename is a mat file, the second one must be in the same format.')
        isMATfile = true;
    else
        isMATfile = false;
    end
end
       
%% filename and output directories

% complete path to snapshot to be read
if ~IsDifference
    filename = {[dataPath, snapshotFilename{1}]};
else
    filename = {...
        [dataPath, snapshotFilename{1}], ...
        [dataPath, snapshotFilename{2}]};
end

% include snapshot filename in name of output directories
% remove extension from snapshot filename: find LAST period (if any)
snapshotFilename_periodIndex = cell(size(snapshotFilename)); % preallocate
snapshotFilename_noext = cell(size(snapshotFilename)); % preallocate
for f=1:size(snapshotFilename, 1) % allows for single and couple-difference
    snapshotFilename_periodIndex{f} = find(snapshotFilename{f}== '.', 1, 'last');
    if ~isempty(snapshotFilename_periodIndex{f})
        snapshotFilename_noext{f} = snapshotFilename{f}(1:snapshotFilename_periodIndex{f} - 1);
    else % there may be no extension, it would be still fine
        snapshotFilename_noext{f} = snapshotFilename{f};
    end
end

% output directories names
if ~IsDifference
    % prepend filename of snapshot
    snapshotFilename_noext_forOutputDir = snapshotFilename_noext{1};
else
    % prepend filename of 2nd snapshot 'minus' filename of 1st snapshot
    snapshotFilename_noext_forOutputDir = [...
        snapshotFilename_noext{2},...
        '_minus_',...
        snapshotFilename_noext{1}];
end
outputGridPath = [dataPath, snapshotFilename_noext_forOutputDir, '_out_grids/'];
outputPSPath = [dataPath, snapshotFilename_noext_forOutputDir, '_out_figures/'];

% output directory for converted ps files
outputConvFigPath = [outputPSPath,'psconverted/'];

%% output format for psconvert: not working right now, using gs as fallaback

% do_psconvert = false;  % TO DO: conversion is failing, call to gs fails!
% temporary fix: gswin64c is called directly

% output format
% OutputFormat = 'png'; % pdf, png supported
% switch OutputFormat
%     case 'pdf'
%         PsConvertTypeString = '-Tf';
%     case 'png'
%         PsConvertTypeString = '-Tg';
%     otherwise
%         error('Unsupported file format');
% end

% TODO: event position and/or rupture area projection

%% set GMT parameters
gmt('gmtset MAP_FRAME_PEN thick,black')
gmt('gmtset MAP_GRID_PEN_PRIMARY thinnest,gray')
gmt('gmtset FONT_ANNOT_PRIMARY 12p')
gmt('gmtset FONT_LABEL 10p')
gmt('gmtset PS_PAGE_ORIENTATION portrait')
gmt('gmtset MAP_TITLE_OFFSET 36p') % avoids overlap with tick labels

% small font size for title, if filenames are included
if extendedTitle % (else: leave GMT defaults)
    gmt('gmtset FONT_TITLE 12p')
end

%% import data
% before loading: do the snapshot file(s) exist?
for f=1:size(filename, 2)
    assert(logical(exist(filename{f}, 'file')),... % exist() returns 2 if true
        ['Provided snapshot filename: ''', filename{f}, ''' does not exist.']);
end

if ~isMATfile
    switch lower(PSCMPQSSPswitch)
        case 'pscmp'
            if ~IsDifference
                snapshotData = PSCMPsnapshot2table(filename{1});
            else
                snapshotData{1} = PSCMPsnapshot2table(filename{1});
                snapshotData{2} = PSCMPsnapshot2table(filename{2});
            end
        case 'qssp'
            if ~IsDifference
                snapshotData = QSSPsnapshot2table(filename{1});
            else
                snapshotData{1} = QSSPsnapshot2table(filename{1});
                snapshotData{2} = QSSPsnapshot2table(filename{2});
            end
    end
else
    if ~IsDifference
        snapshotData = load(filename{1});
        % extract field from structure (expecting 1 field only)
        loaded_fieldnames = fieldnames(snapshotData);
        snapshotData = snapshotData.(loaded_fieldnames{1});
    else
        snapshotData{1} = load(filename{1});
        % extract field from structure (expecting 1 field only)
        loaded_fieldnames = fieldnames(snapshotData{1});
        snapshotData{1} = snapshotData{1}.(loaded_fieldnames{1});
        snapshotData{2} = load(filename{2});
        % extract field from structure (expecting 1 field only)
        loaded_fieldnames = fieldnames(snapshotData{2});
        snapshotData{2} = snapshotData{2}.(loaded_fieldnames{1});
    end
end

%% assert if we are dealing with a grid, extract range and interval
% we need to:
%    1) avoid providing scattered data to gmt xyz2grd
%        (if data is scattered, item 2 would make no sense)
%    2) extract the grid interval (which is a required argument for xyz2grd)
%    3) extract the grid extents (again, a xyz2grd required argument)
% row order:
%    empirically, we know that data in PSCMP snapshots is ordered
%    along LATITUDE first: [ lat0, lon0, ... ; lat1, lon 0, ... ; ... ]
%    the same assumption cannot be made for qssp snapshot files

switch lower(PSCMPQSSPswitch)
    case 'pscmp'
        if ~IsDifference
            [lonRange, latRange, lonStep, latStep] = PSCMPsnapshotCheckGrid(snapshotData);
        else
            [lonRange(1, :), latRange(1, :), lonStep(1), latStep(1)] = PSCMPsnapshotCheckGrid(snapshotData{1});
            [lonRange(2, :), latRange(2, :), lonStep(2), latStep(2)] = PSCMPsnapshotCheckGrid(snapshotData{2});
        end
    case 'qssp'
        if ~IsDifference
            [lonRange, latRange, lonStep, latStep] = QSSPsnapshotCheckGrid(snapshotData);
        else
            [lonRange(1, :), latRange(1, :), lonStep(1), latStep(1)] = QSSPsnapshotCheckGrid(snapshotData{1});
            [lonRange(2, :), latRange(2, :), lonStep(2), latStep(2)] = QSSPsnapshotCheckGrid(snapshotData{2});
        end
end

if IsDifference % assert the two snapshot for same grid
    assert(all([...
        all(lonRange(1, :)==lonRange(2, :)),...
        all(latRange(1, :)==latRange(2, :)),...
        abs(lonStep(1)-lonStep(2)) < 1e4*eps(min(abs(lonStep(1)),abs(lonStep(2)))),... % use tolerance when comparing steps
        abs(latStep(1)-latStep(2)) < 1e4*eps(min(abs(latStep(1)),abs(latStep(2))))]),...
        'The provided snapshot have different grids. Check extents and step.')
    lonRange = lonRange(1, :);
    latRange = latRange(1, :);
    lonStep = lonStep(1);
    latStep = latStep(1);
end

%% call xyz2grd
if ~IsDifference
    snapshotGrids = snapshotTable2gmtMexGrid(...
        snapshotData, lonStep, latStep, lonRange, latRange);
    % store names of grids
    ObservableNames = fieldnames(snapshotGrids);
else
    snapshotGrids{1} = snapshotTable2gmtMexGrid(...
        snapshotData{1}, lonStep, latStep, lonRange, latRange);
    snapshotGrids{2} = snapshotTable2gmtMexGrid(...
        snapshotData{2}, lonStep, latStep, lonRange, latRange);
    % store names of grids and assert for equality
    ObservableNames = fieldnames(snapshotGrids{1});
    assert(all(strcmp(ObservableNames, fieldnames(snapshotGrids{2}))),...
        'Observable names in the two provided snapshots are different.')
end

%% compute free-air-correction-removal
% psgrn/pscmp and qssp apply a free-air correction to model
% the movement of a gravimeter placed on the ground (which is deformed)
% to use the gravity data for fixed height modelling (e.g. sat altitude)
% we need to REMOVE the applied correction.
%
% two cases: 1) we are dealing with pscmp, 2) we are dealing with qssp
%
% if among the observables we have 'Gravity': pscmp snapshot
%    use 'Disp_down' displacement, reverse sign as positive upwards
% if among the observables we have 'Grav': qsspstatic snapshot
%    use 'U_z' displacement, positive upwards as it is

% we have to counter-correct the applied correction
% therefore we SUBTRACT it: points that have risen get a POSITIVE counter-correction

% behave differently for pscmp and qssp data
switch lower(PSCMPQSSPswitch)
    case 'pscmp'
        assert(any(strcmp(ObservableNames,'Disp_down')),...
            'Expecting ''Disp_down'' field, but none was found.');
        if ~IsDifference
            snapshotGrids.Gravity.z = countercorrectFreeAir_pscmp_qssp(...
                snapshotGrids.Gravity.z, -snapshotGrids.Disp_down.z);
        else
            snapshotGrids{1}.Gravity.z = countercorrectFreeAir_pscmp_qssp(...
                snapshotGrids{1}.Gravity.z, -snapshotGrids{1}.Disp_down.z);
            snapshotGrids{2}.Gravity.z = countercorrectFreeAir_pscmp_qssp(...
                snapshotGrids{2}.Gravity.z, -snapshotGrids{2}.Disp_down.z);
        end
    case 'qssp'
        assert(any(strcmp(ObservableNames,'U_z')),...
            'Expecting ''U_z'' field, but none was found.');
        if ~IsDifference
            snapshotGrids.Grav.z = countercorrectFreeAir_pscmp_qssp(...
                snapshotGrids.Grav.z, snapshotGrids.U_z.z);
        else
            snapshotGrids{1}.Grav.z = countercorrectFreeAir_pscmp_qssp(...
                snapshotGrids{1}.Grav.z, snapshotGrids{1}.U_z.z);
            snapshotGrids{2}.Grav.z = countercorrectFreeAir_pscmp_qssp(...
                snapshotGrids{2}.Grav.z, snapshotGrids{2}.U_z.z);
        end
end

% update z range in gmt-mex grid struct
% TO DO: do switch once, set name of 'gravity' field and sign of displacement
if ~IsDifference
    switch lower(PSCMPQSSPswitch)
        case 'pscmp'
            snapshotGrids.Gravity.range(5:6) = ...
                [min(snapshotGrids.Gravity.z(:)),...
                max(snapshotGrids.Gravity.z(:))];
        case 'qssp'
            snapshotGrids.Grav.range(5:6) = ...
                [min(snapshotGrids.Grav.z(:)),...
                max(snapshotGrids.Grav.z(:))];
    end
end
% if IsDifference: no need to update z range in gmt-mex grid struct
% since we are computing differences (and new z-ranges) soon

%% define units and conversion factor for each observable

% TO DO: using minor tick interval computation, determine multiplier

% TO DO: units and conversion factor in same struct
%        maybe even in grid structure

% this is done manually (at the moment no way to do otherwise)
% then equality with existing observables is asserted

% 'all units are SI unless specified' (source: example input files)

% common units
commonUnits.Displacement = 'mm';
commonUnits.Stress = 'kPa';
commonUnits.Tilt = '@~m@~rad'; % microrad
commonUnits.Rotation = commonUnits.Tilt;
commonUnits.Geoid = '@~m@~m'; % micrometers
commonUnits.Gravity = '@~m@~Gal'; % microGal

% common conversion factors
commonConv.Displacement = 1e3;
commonConv.Stress = 1e-3;
commonConv.Tilt = 1e6;
commonConv.Rotation = commonConv.Tilt;
commonConv.Geoid = 1e6;
commonConv.Gravity = 1e8; % m/s2 to microGal

switch lower(PSCMPQSSPswitch)
    case 'pscmp'
        % units strings
        Units.Disp_north = commonUnits.Displacement;
        Units.Disp_east = commonUnits.Displacement;
        Units.Disp_down = commonUnits.Displacement;
        Units.Stress_nn = commonUnits.Stress;
        Units.Stress_ee = commonUnits.Stress;
        Units.Stress_dd = commonUnits.Stress;
        Units.Stress_ne = commonUnits.Stress;
        Units.Stress_ed = commonUnits.Stress;
        Units.Stress_dn = commonUnits.Stress;
        Units.Tilt_n = commonUnits.Tilt;
        Units.Tilt_e = commonUnits.Tilt;
        Units.Rotation = commonUnits.Rotation;
        Units.Geoid = commonUnits.Geoid;
        Units.Gravity = commonUnits.Gravity;
        Units.Disp_LOS = commonUnits.Displacement;
        % conversion factors
        Conv.Disp_north = commonConv.Displacement;
        Conv.Disp_east = commonConv.Displacement;
        Conv.Disp_down = commonConv.Displacement;
        Conv.Stress_nn = commonConv.Stress;
        Conv.Stress_ee = commonConv.Stress;
        Conv.Stress_dd = commonConv.Stress;
        Conv.Stress_ne = commonConv.Stress;
        Conv.Stress_ed = commonConv.Stress;
        Conv.Stress_dn = commonConv.Stress;
        Conv.Tilt_n = commonConv.Tilt;
        Conv.Tilt_e = commonConv.Tilt;
        Conv.Rotation = commonConv.Rotation;
        Conv.Geoid = commonConv.Geoid;
        Conv.Gravity = commonConv.Gravity;
        Conv.Disp_LOS = commonConv.Displacement;
    case 'qssp'
        % units strings
        Units.Station = '';
        Units.U_n = commonUnits.Displacement;
        Units.U_e = commonUnits.Displacement;
        Units.U_z = commonUnits.Displacement;
        Units.Vstrain = '10@+-6@+'; % adminesional
        Units.Grav = commonUnits.Gravity;
        Units.Geoid = commonUnits.Geoid;
        Units.Tilt_n = commonUnits.Tilt;
        Units.Tilt_e = commonUnits.Tilt;
        % conversion factors
        Conv.Station = 1;
        Conv.U_n = commonConv.Displacement;
        Conv.U_e = commonConv.Displacement;
        Conv.U_z = commonConv.Displacement;
        Conv.Vstrain = 1e6; % m3 to cm3
        Conv.Grav = commonConv.Gravity;
        Conv.Geoid = commonConv.Geoid;
        Conv.Tilt_n = commonConv.Tilt;
        Conv.Tilt_e = commonConv.Tilt;
end

% TO DO: assert equality of defined fields in Units
% note: strcmp(fieldnames, fieldnames) imposes same order
% which is not required - the requirement is that Units CONTAINS
% all ObservableNames

%% compute snapshot-couple difference
if IsDifference
    snapshotGridsDifferences = snapshotGrids{1}; % to get the same grid structures
    for n=1:size(ObservableNames, 1)
        snapshotGridsDifferences.(ObservableNames{n}).z = ...
            snapshotGrids{2}.(ObservableNames{n}).z - snapshotGrids{1}.(ObservableNames{n}).z;
        snapshotGridsDifferences.(ObservableNames{n}).range(5:6) = ...
            [min(snapshotGridsDifferences.(ObservableNames{n}).z(:)), ...
            max(snapshotGridsDifferences.(ObservableNames{n}).z(:))];
        snapshotGridsDifferences.(ObservableNames{n}).z_unit = Units.(ObservableNames{n});
    end
    % from now on: even if IsDifference is true,
    % snapshotGrids is a scalar grid structure
    % as if we had only one snapshot do plot
    snapshotGrids = snapshotGridsDifferences;
    
end

%% perform conversion (e.g. m to mm, m/s2 to mGal, ...)
% multiplication by a factor AND update gmt mex grid structure min/max
% with a function!
for n=1:size(ObservableNames, 1)
    snapshotGrids.(ObservableNames{n}).z = ...
        snapshotGrids.(ObservableNames{n}).z * Conv.(ObservableNames{n});
    snapshotGrids.(ObservableNames{n}).range(5:6) = ...
        [min(snapshotGrids.(ObservableNames{n}).z(:)), ...
        max(snapshotGrids.(ObservableNames{n}).z(:))];
    snapshotGrids.(ObservableNames{n}).z_unit = Units.(ObservableNames{n});
end

%% create output folders, if they do not exist yet
if ~exist(outputGridPath, 'dir')
    mkdir(outputGridPath)
end
if ~exist(outputPSPath, 'dir')
    mkdir(outputPSPath)
end
if ~exist(outputConvFigPath, 'dir')
    mkdir(outputConvFigPath)
end

%% write grd files
for n=1:size(ObservableNames, 1)
    gmtMexGrid2grd(...
        snapshotGrids.(ObservableNames{n}),...
        [outputGridPath, (ObservableNames{n}), '.grd'])
end

%% return if noPlotsFlag is true
if noPlotsFlag
    disp([...
        snapshotFilename_noext_forOutputDir,...
        ' : noPlotsFlag is True, skipping plotting of maps'])
    return
end

%% colorscales
CPT_colorScale = 'roma -E -Z'; % E, Z options: linear continous colorscale

% quantiles to clip min/max limits of colorscales at (to avoid outliers)
%MinQuantile = 0.001;
%MaxQuantile = 0.999;

for n=1:size(ObservableNames, 1)
    % skip zero-z-range grid, to avoid grd2cpt errors
    if ~((snapshotGrids.(ObservableNames{n}).range(6) - snapshotGrids.(ObservableNames{n}).range(5)) == 0)
        % compute 'optimal' minor tick interval (TO DO: could be better!)
        MinorTickInt.(ObservableNames{n}) = ...
            RoundToStep(0.1, ...
                (snapshotGrids.(ObservableNames{n}).range(6) - ...
                snapshotGrids.(ObservableNames{n}).range(5)) / 25,...
                'ceil');
        if MinorTickInt.(ObservableNames{n}) > 10
            MinorTickInt.(ObservableNames{n}) = RoundToStep(...
                5, MinorTickInt.(ObservableNames{n}), 'ceil');
        end
        MajorTickInt.(ObservableNames{n}) = MinorTickInt.(ObservableNames{n})*5;
        if MajorTickInt.(ObservableNames{n}) > 5
            MajorTickInt.(ObservableNames{n}) = RoundToStep(...
                5, MajorTickInt.(ObservableNames{n}), 'ceil');
            MinorTickInt.(ObservableNames{n}) = MajorTickInt.(ObservableNames{n}) / 5;
        end
        CPTs.(ObservableNames{n}) = gmt([...
            'grd2cpt -C', CPT_colorScale,...
            ' -L',...
            num2str(RoundToStep(...
                MinorTickInt.(ObservableNames{n}),...
                quantile(snapshotGrids.(ObservableNames{n}).z(:),MinQuantile),...
                'floor')),'/',...
            num2str(RoundToStep(...
                MinorTickInt.(ObservableNames{n}),...
                quantile(snapshotGrids.(ObservableNames{n}).z(:),MaxQuantile),...
                'ceil')),...
            ' -Di'],...
            snapshotGrids.(ObservableNames{n}));
        % avoid case in which the colorscale covers less than a 2 major ticks range
        l = 0; % re-set loop counter
        while (CPTs.(ObservableNames{n}).minmax(2) - CPTs.(ObservableNames{n}).minmax(1)) < (MajorTickInt.(ObservableNames{n}) * 2)
            MajorTickInt.(ObservableNames{n}) = MinorTickInt.(ObservableNames{n}) / 2;
            MinorTickInt.(ObservableNames{n}) = MinorTickInt.(ObservableNames{n}) / 10;
            % re-call grd2cpt
            CPTs.(ObservableNames{n}) = gmt([...
                'grd2cpt -C', CPT_colorScale,...
                ' -L',...
                num2str(RoundToStep(...
                MinorTickInt.(ObservableNames{n}),...
                quantile(snapshotGrids.(ObservableNames{n}).z(:),MinQuantile),...
                'floor')),'/',...
                num2str(RoundToStep(...
                MinorTickInt.(ObservableNames{n}),...
                quantile(snapshotGrids.(ObservableNames{n}).z(:),MaxQuantile),...
                'ceil')),...
                ' -Di'],...
                snapshotGrids.(ObservableNames{n}));
            l = l+1;
            if l>4
                break
            end
        end
        % avoid case in which there are 'too many' divisions in the cpt range
        l = 0; % re-set loop counter
        while ((CPTs.(ObservableNames{n}).minmax(2) - CPTs.(ObservableNames{n}).minmax(1)) / MinorTickInt.(ObservableNames{n})) > 35
            MajorTickInt.(ObservableNames{n}) = MajorTickInt.(ObservableNames{n}) * 2;
            MinorTickInt.(ObservableNames{n}) = MinorTickInt.(ObservableNames{n}) * 2;
            l = l+1;
            if l>4
                break
            end
        end
    end
end

%% define region, -R and -J arguments (extents and projection)
RegionExtents = [lonRange, latRange]; % data extents [lon0 lon1 lat0 lat1]
if ForceExtents % [lon0 lon1 lat0 lat1] were provided, use those
    MapRegionExtents = ExtentsBuffer;
else % edge buffer/trim was provided
    switch SourceType % should we refer them to the grid coverage or around the centre point / fault plane?
        case 'none'
            % referring the edge buffer/trim to the grid extents 'RegionExtents'
            MapRegionBuffer = [...
                -edgeBufferLon, edgeBufferLon,...
                -edgeBufferLat, edgeBufferLat]; % buffer/clip around data [lon0 lon1 lat0 lat1]
            MapRegionExtents = RegionExtents + MapRegionBuffer; % map corners [lon0 lon1 lat0 lat1]
            MapRegion_isPointReferred = false;
        case 'point'
            % referring the buffer around the centre point
            MapRegion_isPointReferred = true;
            MapRegion_sourceLonLat = sourceLonLat;
        case 'rectangle'
            % referring the buffer around the fault rectangle
            MapRegion_isPointReferred = true;
            % compute the rectangle center
            %    = along column average: [average lon, average lat]
            MapRegion_sourceLonLat = mean(sourceLonLat, 1);
    end
    % avoid repetitions: the following apply if we are 'referring to point' (either point or center of rectangle)
    if MapRegion_isPointReferred
        % negative 'trim' is nonsense in this case
        edgeBufferLon = abs(edgeBufferLon);
        edgeBufferLat = abs(edgeBufferLat);
        MapRegionExtents = [...
            MapRegion_sourceLonLat(1) - edgeBufferLon, MapRegion_sourceLonLat(1) + edgeBufferLon,...
            MapRegion_sourceLonLat(2) - edgeBufferLat, MapRegion_sourceLonLat(2) + edgeBufferLat];
    end
end

% grid lines: if extents (in any direction) are less than 'GridLines_FineIntervalExtentsThreshold'
%             use finer grid line interval
if (abs(MapRegionExtents(2) - MapRegionExtents(1)) <=GridLines_FineIntervalExtentsThreshold || ...
        abs(MapRegionExtents(4) - MapRegionExtents(3)) <=GridLines_FineIntervalExtentsThreshold)
    Map_LonLatGridInterval = GridLines_FineInterval;
else
    Map_LonLatGridInterval = GridLines_CoarseInterval;
end

% plot countries only if the lon/lat extents are larger than 'PlotCountries_ExtentsThreshold'
if (abs(MapRegionExtents(2) - MapRegionExtents(1)) <=PlotCountries_ExtentsThreshold || ...
        abs(MapRegionExtents(4) - MapRegionExtents(3)) <=PlotCountries_ExtentsThreshold)
    Map_PlotCountries = false;
else
    Map_PlotCountries = true;
end

Rstring = [... % extents, r option means lower left and upper right coords
    '-R',num2str(MapRegionExtents(1)),'/',num2str(MapRegionExtents(3)),'/',...
    num2str(MapRegionExtents(2)),'/',num2str(MapRegionExtents(4)),'r ']; % extents
Jstring = ['-JA',...
    num2str(RegionExtents(1)+(RegionExtents(2)-RegionExtents(1))/2),'/',... % lon0 center
    num2str(RegionExtents(3)+(RegionExtents(4)-RegionExtents(3))/2),...     % lat0 center
    '/6i']; % projection (-JA Lambert lon0/lat0)
RJstring = [Rstring, Jstring]; % extents and projection, concatenated

%% map titles: observable name only or optional 'extended title'
% the extended title includes the snapshotFilename(s)

if ~extendedTitle
    ObservableTitles = ObservableNames;
else
    % uses the same scheme used for paths (single file and couple difference)
    ObservableTitles = cell(size(ObservableNames));
    % format: '${observable} (${file or couple difference name})'
    for n=1:size(ObservableNames, 1)
        ObservableTitles{n} = ...
            [ObservableNames{n},...
            ' (', snapshotFilename_noext_forOutputDir, ')'];
    end
end


%% figure (TO DO: in function, useful also for non-snapshot files)

% if do_psconvert==true
% disp('Note that GMT psconvert fails if the target file is open in another program!')
% end

CommonGMT_start = tic;
for m=1:size(ObservableNames, 1)
    GMT_progress = ['[',num2str(m,'%02.0f'),'/',...
        num2str(size(ObservableNames, 1),'%02.0f'),']'];
    GMT_start = tic;
    MapFilename = [outputPSPath,ObservableNames{m},'.ps'];
    MapGrid = snapshotGrids.(ObservableNames{m});
    if (snapshotGrids.(ObservableNames{m}).range(6) - snapshotGrids.(ObservableNames{m}).range(5)) == 0
        disp([GMT_progress,' data has zero range, skipping: ', ObservableNames{m}]);
    else
    
        MapCPT = CPTs.(ObservableNames{m});
        MapZUnit = Units.(ObservableNames{m});

        MapScaleWidth = '5.0i';
        MapScaleTriangles = '+e';
        
        MapMinorTickInt = MinorTickInt.(ObservableNames{m});
        MapMajorTickInt = MapMinorTickInt * 5;

        % grdimage and title
        gmt(['grdimage ',RJstring,...
            ' -E -Q -nc+c -C -K -P -Xc -Yc -B+t"',ObservableTitles{m},'" > ',...
            MapFilename],...
            MapGrid,MapCPT);

        % countries
        if Map_PlotCountries % if extents larger thant PlotCountries_ExtentsThreshold
            gmt(['pscoast -R -J',...
                ' -P -N1/0.01c,Black',...
                ' -O -K -Xc -Yc >> ',...
                MapFilename])
        end
        % coasts
        gmt(['pscoast -R -J',...
            ' -P -Di+ -A20000',...
            ' -W0.03c,Black -O -K -Xc -Yc >> ',...
            MapFilename])
        
        % contours
        if doContours
            gmt(['grdcontour -R -J',...
                ' -C',num2str(MapMinorTickInt)...
                ' -A',num2str(MapMajorTickInt),'+f13p'...
                ' -S',num2str(8),...
                ' -Gn1+r15p',... % +r: enforce minimum separation radius bw labels
                ' -t50 -O -K -P -Xc -Yc -Wathin -Wcfaint >> ',...
                MapFilename],...
                MapGrid);
        end

        % semi-trasparent grid
        gmt(['psbasemap -R -J',...
            ' -B', num2str(Map_LonLatGridInterval), 'g' , num2str(Map_LonLatGridInterval),... % e.g. -B1g1
            ' -t10 -BWesN -O -K -Xc -Yc >> ',...
            MapFilename]);

        % non-trasparent frame and ticks
        gmt(['psbasemap -R -J',...
            ' -B', num2str(Map_LonLatGridInterval), ' -BWesN -O -K -Xc -Yc >> ',...
            MapFilename]);
        
        switch SourceType % source dot or fault area rectangle (hor projection)
            case 'point'
                gmt(['psxy -R -J',...
                    ' -Sc10p -Gpurple',...
                    ' -O -K -Xc -Yc >> ',...
                    MapFilename], sourceLonLat);
            case 'rectangle'
                gmt(['psxy -R -J',...
                    ' -A -L -Wthick -Wpurple',... % -A straight-line segments
                    ' -O -K -Xc -Yc >> ',...
                    MapFilename], [sourceLonLat; sourceLonLat(1,:)]); % repeat first point, close rectangle
        end
        
        % min/max range text
        MapStatsFontSize = '12p';
        MapStatsXshift = '0i';
        MapStatsYshift = '-3.60i';
        StatsText = [...
            'min = ',num2str(MapGrid.range(5), '%6.3f'),...
            ', max = ',num2str(MapGrid.range(6), '%6.3f')];
        gmt(['pstext ',...
            ' -R1/10/1/10 -JX10 -F+f',MapStatsFontSize,'+cCM',... % +cCM : justification
            ' -P -O -K -Xc',MapStatsXshift,' -Yc',MapStatsYshift,...
            ' >> ',MapFilename], StatsText)
        
        % scale
        if isempty(MapZUnit)
            ByArgs = ' -By';
        else
            ByArgs = ' -By+l';
        end
        MapScaleYcshift = '-0.5i';
        gmt(['psscale ',...
            '-DJBC+o0/0.7i+w',MapScaleWidth,MapScaleTriangles,' -R -JX',MapScaleWidth,' -O -P ',...
            ' -Bx',num2str(MapMajorTickInt),'f',num2str(MapMinorTickInt),...
            ByArgs,MapZUnit,' -Xc -Yc',MapScaleYcshift,' >> ',...
            MapFilename], MapCPT)

        GMT_elap = toc(GMT_start);
        disp([GMT_progress,' GMT: ''',MapFilename,...
            ''' written in ',num2str(GMT_elap,'%.3f'),' s']);
%        GMT_start = tic;
%         if do_psconvert==true
%             gmt(['psconvert ',MapFilename,' -A -P ',PsConvertTypeString,' -D',outputConvFigPath])
%             GMT_elap = toc(GMT_start);
%             disp([GMT_progress,' GMT: ''',MapFilename,...
%                 ''' psconverted to ',OutputFormat,' in ',num2str(GMT_elap,'%.3f'),' s']);
%         end
        % temporary workaround: perform conversion with gs directly
        % note that by doing so we lose some gmt psconvert specifics
        % such as crop-to-bounding-box
        [gsCallStatus, gsCallCmdOut] = system([...
            gsPath, ' -dNOSAFER -dNOPAUSE -dBATCH -DPSL_no_pagefill ',...
            '-sDEVICE=png16m -dUseCropBox -r200 -sOutputFile=',...
            outputConvFigPath,ObservableNames{m},'.png ', MapFilename]);
        if gsCallStatus==0
            disp([GMT_progress,' GMT: ''',MapFilename,...
                ''' converted to png.']);
        else
            disp([GMT_progress,' GMT: ''',MapFilename,...
                ''' conversion returned non-zero exit status.']);
            disp(gsCallStatus)
            disp(gsCallCmdOut)
        end
    end
end
disp(['[DONE!] GMT call done in ',num2str(toc(CommonGMT_start),'%.3f'),' s'])

%% nested function: round to nearest value according to step
    function out = RoundToStep(step,in,direction)
        %RoundToStep round(/floor/ceil/fix) input vector to nearest arbitrary 'step' unit
        % digit-counting code comes from this answer by user Jaymin on matlabcentral:
        % https://mathworks.com/matlabcentral/answers/10795-counting-the-number-of-digits#answer_68482
        %
        % Syntax: out = RoundToStep(step,in,[direction])
        %
        % Inputs:
        %    step        : scalar, unit value
        %    in          : vector of values to be rounded
        %    [direction] : char vector, optional, direction of rounding
        %                  'round', 'floor', 'ceil', 'fix'
        %                  'round' is default behaviour
        %                  case insensitive
        %
        % Outputs:
        %    out         : rounded 'in' vector
        %

        
        % check and manage input
        narginchk(2,3)
        if nargin==2
            direction='round';
        end
        assert(isscalar(step),'Argument ''step'' must be a scalar.')
        assert(isvector(in),'Argument ''in'' must be a vector.')
        % check if 'direction' is an allowed round-like function
        allowed_directions = {'round','floor','ceil','fix'};
        assert(any(strcmpi(allowed_directions,direction)),...
            ['''',direction,''' is not an allowed rounding function.'])
        
        % count the number of digits
        stepDIG = numel(num2str(step))-numel(strfind(num2str(step),'-'))-numel(strfind(num2str(step),'.'));
        
        % number of zsteps in next power of 10
        stepdiv = (10^(stepDIG))/step;
        
        % round/floor/ceil/fix to step
        roundfunc = str2func(lower(direction)); % create handle to selected rounding function
        out = ((roundfunc((in/(10^(stepDIG)))*stepdiv))/stepdiv)*(10^(stepDIG));
        
    end

gmt('destroy') % housekeeping, free memory

end
