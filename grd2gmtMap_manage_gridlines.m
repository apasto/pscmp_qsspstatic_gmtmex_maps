function Map_LonLatGridInterval = grd2gmtMap_manage_gridlines(MapRegionExtents, varargin)
% grd2gmtMap_manage_gridlines

% static parameters for plotting fine/coarse grid line interval
% GridLines_FineIntervalExtentsThreshold: if smaller go from 'coarse' to 'fine' interval

narginchk(1, 4)
if nargin>1 && ~isempty(varargin{1})
    GridLines_FineIntervalExtentsThreshold = varargin{1};
else
    GridLines_FineIntervalExtentsThreshold = 2.25; % [deg]
end
if nargin>2 && ~isempty(varargin{2})
    GridLines_FineInterval = varargin{2};
else
    GridLines_FineInterval = 0.5; % [deg]
end
if nargin>3 && ~isempty(varargin{3})
    GridLines_CoarseInterval = varargin{3};
else
    GridLines_CoarseInterval = 1.0; % [deg]
end

% grid lines: if extents (in any direction) are less than 'GridLines_FineIntervalExtentsThreshold'
%             use finer grid line interval
if (abs(MapRegionExtents(2) - MapRegionExtents(1)) <=GridLines_FineIntervalExtentsThreshold || ...
        abs(MapRegionExtents(4) - MapRegionExtents(3)) <=GridLines_FineIntervalExtentsThreshold)
    Map_LonLatGridInterval = GridLines_FineInterval;
else
    Map_LonLatGridInterval = GridLines_CoarseInterval;
end
end

