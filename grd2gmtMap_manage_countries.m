function Map_PlotCountries = grd2gmtMap_manage_countries(MapRegionExtents, varargin)
%grd2gmtMap_manage_countries

narginchk(1, 2)
if nargin>1 && ~isempty(varargin{1})
    PlotCountries_ExtentsThreshold = varargin{1};
else
    PlotCountries_ExtentsThreshold = 6.0; % [deg]
end

% plot countries only if the lon/lat extents are larger than 'PlotCountries_ExtentsThreshold'
if (abs(MapRegionExtents(2) - MapRegionExtents(1)) <=PlotCountries_ExtentsThreshold || ...
        abs(MapRegionExtents(4) - MapRegionExtents(3)) <=PlotCountries_ExtentsThreshold)
    Map_PlotCountries = false;
else
    Map_PlotCountries = true;
end

end

