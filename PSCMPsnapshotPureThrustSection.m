function [Figures, Axes, Plots] = PSCMPsnapshotPureThrustSection(...
    snapshotFile, profileLon, faultOriginLat, faultOriginDepth, faultDip, faultWidth, plotLatWidth, varargin)
%PSCMPsnapshotPureThrustSection ad-hoc N-S section plot for a Pure Thrust
% NOTE: this is NOT written for a general case, which would require a general
%       reconstruction of the rupture plane using strike and dip.
%       The provided faultDip is ASSUMED SOUTH-DIPPING and with 90° strike.
%
%   Input arguments:
%      - snapshotFile :
%      - profileLon :
%      - faultOriginLat :
%      - faultOriginDepth : [km]
%      - faultDip : [deg]
%      - faultWidth : [km]
%      - plotLatWidth : [deg] symmetrical plot extents around the fault
%                       midpoint along latitude
%                       (i.e. how wide should we plot the section, around
%                       the fault origin)
%      - saveFilename : path to figure to be saved (optional)
%                       complete filename will be:
%                       saveFilename + '_' + faultDip + 'deg_' + faultWidth + 'km.png'
%
%   Output arguments:
%      - figure, axex, plots handles (optional)
%
% 2021-01-24 AP

narginchk(7,8)
nargoutchk(0,3)

% optional argument: save figure to file
if nargin>7
    PlotFigure = true;
    saveFilename = varargin{1};
else
    PlotFigure = false;
end

% load the snapshot
snapshotData = PSCMPsnapshot2table(snapshotFile);

% find the closest longitude to profileLon
uniqueLons = uniquetol(snapshotData.Londeg);
[minLonDiff, minLonDiff_index] = min(abs(uniqueLons - profileLon)); %#ok<ASGLU>
profileLonActual = uniqueLons(minLonDiff_index);

% extract indices of rows at the actual profile longitude
% to do: include tolerance? Equality may be dangerous (precision!)
profileIndices = find(snapshotData.Londeg==profileLonActual);

% extract latitudes along profile
profileLats = snapshotData.Latdeg(profileIndices);
% latitudes are not guaranteed to be sorted! (depends on ordering in file)
[profileLatsSorted, profileLatsSorted_index] = sort(profileLats);
% use the obtained sorting to arrenge the profileIndices
% which we are then using to extract the data (gravity, displacement) along profile
profileIndicesSorted = profileIndices(profileLatsSorted_index);

% extract the data along profile
profileGravity = snapshotData.Gravity(profileIndicesSorted);
profileUz = -snapshotData.Disp_down(profileIndicesSorted); % change sign: positive upwards

% correct for the applied Free Air correction (i.e. counter-correct it)
profileGravity = countercorrectFreeAir_pscmp_qssp(...
    profileGravity, profileUz);

% convert: gravity to microGal, displacement to mm
profileGravity = profileGravity * 1e8;
profileUz = profileUz * 1e3;

% rupture plane: construct line using origin, dip, along-dip width
faultBottomDepth = faultOriginDepth + faultWidth * sin(deg2rad(faultDip));
% latitude of fault bottom - note assumptions: south dipping
% get Earth radius for length to arc-along-latitude conversion
R = referenceSphere('Earth').MeanRadius * 1e-3; % [km]
faultHorizontalWidth = faultWidth * cos(deg2rad(faultDip)); % projection along horizontal
faultHorizontalLatArc = rad2deg(faultHorizontalWidth / R);
faultBottomLat = faultOriginLat - faultHorizontalLatArc;
faultMidpointLat = faultOriginLat - (faultHorizontalLatArc / 2);

% figure
PlotLineWidth = 1.5;
Figures = figure();
Figures.Position = [400, 150, 900, 700];

Axes.Gravity = subplot(5, 1, [1, 2], 'Parent', Figures);
Axes.Displacement = subplot(5, 1, [3, 4], 'Parent', Figures);
Axes.Fault = subplot(5, 1, 5, 'Parent', Figures);

Axes_names = fieldnames(Axes);
for n=1:size(Axes_names, 1)
    Axes.(Axes_names{n}).Box = 'on';
    Axes.(Axes_names{n}).XGrid = 'on';
    Axes.(Axes_names{n}).YGrid = 'on';
    Axes.(Axes_names{n}).Layer = 'top';
    hold(Axes.(Axes_names{n}), 'on');
end

Axes.Gravity.XAxisLocation = 'top';
Axes.Displacement.XTickLabel = [];

Axes.Gravity.XLabel.String = 'Latitude [deg]';
Axes.Fault.XLabel.String = 'Latitude [deg]';

Axes.Gravity.YLabel.String = '\it{g} \rm{change} [\muGal]';
Axes.Displacement.YLabel.String = 'Vert. displacement (up) [mm]';
Axes.Fault.YLabel.String = 'Depth [km]';

Axes.Fault.YDir = 'reverse'; % depth : downwards positive

Plots.Gravity = plot(profileLatsSorted, profileGravity,...
    'Parent', Axes.Gravity,...
    'LineWidth', PlotLineWidth);
Plots.Displacement = plot(profileLatsSorted, profileUz,...
    'Parent', Axes.Displacement,...
    'LineWidth', PlotLineWidth);
Plots.Fault = plot(...
    [faultOriginLat, faultBottomLat],...
    [faultOriginDepth, faultBottomDepth],...
    'Parent', Axes.Fault,...
    'LineWidth', PlotLineWidth);

linkaxes([Axes.Gravity, Axes.Displacement, Axes.Fault], 'x');

Axes.Gravity.XLim = [...
    faultMidpointLat - plotLatWidth,...
    faultMidpointLat + plotLatWidth];

Axes.Gravity.YLimMode = 'manual'; % prevent update of Y range

Axes.Fault.YLim = [0, faultBottomDepth * 1.05];

% projection of the fault at the surface: fill area
Plots.Gravity_ProjFill = fill(...
    [faultOriginLat, faultBottomLat,...
        faultBottomLat, faultOriginLat],...
    [Axes.Gravity.YLim(1), Axes.Gravity.YLim(1),...
        Axes.Gravity.YLim(2), Axes.Gravity.YLim(2)],...
    [0.8, 0.8, 0.8],...
    'EdgeColor', 'none',...
    'Parent', Axes.Gravity);
Plots.Displacement_ProjFill = fill(...
    [faultOriginLat, faultBottomLat,...
faultBottomLat, faultOriginLat],...
    [Axes.Displacement.YLim(1), Axes.Displacement.YLim(1),...
        Axes.Displacement.YLim(2), Axes.Displacement.YLim(2)],...
    [0.8, 0.8, 0.8],...
    'EdgeColor', 'none',...
    'Parent', Axes.Displacement);
Plots.Fault_ProjFill = fill(...
    [faultOriginLat, faultBottomLat,...
        faultBottomLat, faultOriginLat],...
    [Axes.Fault.YLim(1), Axes.Fault.YLim(1),...
        Axes.Fault.YLim(2), Axes.Fault.YLim(2)],...
    [0.8, 0.8, 0.8],...
    'EdgeColor', 'none',...
    'Parent', Axes.Fault);

% bring plots back on top
uistack(Plots.Gravity, 'top');
uistack(Plots.Displacement, 'top');
uistack(Plots.Fault, 'top');

% set aspect ratio for fault plot, to keep true angles
% x axis in Lat [degrees], y axis in depth [km]
% x:y ratio should be such as (1 deg in km) : (1 km)
% or viceversa (1 deg) : (1 km in deg)
% which we can also express as (1 deg in rad) : (1 km in rad)
Lat2DepthAspectRatio = (pi/180) / (1 / R);
Axes.Fault.DataAspectRatio = [1, Lat2DepthAspectRatio, 1];
% note that if figure is set too wide, Axes.Fault will be less wide
% than the two top axes (gravity and displacement), to fit the
% imposed DataAspectRatio

% Axes.Fault fontsize changes unexpectedly, re-set it to the common size
Axes.Fault.FontSize = Axes.Gravity.FontSize;

if PlotFigure
    % normally extension is added by 'print', it is not happening here...
    % saveFilename + '_' + faultDip + 'deg_' + faultWidth + 'km.png'
    complete_saveFilename = [...
        saveFilename, '_', num2str(faultDip), 'deg_',...
        num2str(faultWidth), 'km.png'];
    print(Figures,complete_saveFilename, '-dpng', '-noui')
end

end
