function [MapRegionExtents, RJstring] = grd2gmtMap_create_RJ_options(...
    lonRange, latRange, ExtentsBuffer, ForceExtents,...
    SourceType, sourceLonLat,...
    edgeBufferLon, edgeBufferLat)
%grd2gmtMap_create_RJ_options

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
    % avoid repetitions: the following apply if we are 'referring to point'
    % (either point or center of rectangle)
    if MapRegion_isPointReferred
        % negative 'trim' is nonsense in this case
        edgeBufferLon = abs(edgeBufferLon);
        edgeBufferLat = abs(edgeBufferLat);
        MapRegionExtents = [...
            MapRegion_sourceLonLat(1) - edgeBufferLon, MapRegion_sourceLonLat(1) + edgeBufferLon,...
            MapRegion_sourceLonLat(2) - edgeBufferLat, MapRegion_sourceLonLat(2) + edgeBufferLat];
    end
end

Rstring = [... % extents, r option means lower left and upper right coords
    '-R',num2str(MapRegionExtents(1)),'/',num2str(MapRegionExtents(3)),'/',...
    num2str(MapRegionExtents(2)),'/',num2str(MapRegionExtents(4)),'r ']; % extents
Jstring = ['-JA',...
    num2str(RegionExtents(1)+(RegionExtents(2)-RegionExtents(1))/2),'/',... % lon0 center
    num2str(RegionExtents(3)+(RegionExtents(4)-RegionExtents(3))/2),...     % lat0 center
    '/6i']; % projection (-JA Lambert lon0/lat0)
RJstring = [Rstring, Jstring]; % extents and projection, concatenated

end