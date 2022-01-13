function [ForceExtents, edgeBufferLon, edgeBufferLat] = grd2gmtMap_manage_buffer_extents(ExtentsBuffer)
%grd2gmtMap_manage_buffer_extents

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

end