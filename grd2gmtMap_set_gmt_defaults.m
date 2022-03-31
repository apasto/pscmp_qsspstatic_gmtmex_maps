function grd2gmtMap_set_gmt_defaults(varargin)
%grd2gmtMap_set_gmt_defaults

narginchk(0, 1)

if nargin==1 && ~isempty(varargin{1})
    small_title = true;
else
    small_title = false;
end

gmt('gmtset MAP_FRAME_PEN thick,black')
gmt('gmtset MAP_GRID_PEN_PRIMARY thinnest,gray')
gmt('gmtset FONT_ANNOT_PRIMARY 12p')
gmt('gmtset FONT_LABEL 10p')
gmt('gmtset PS_PAGE_ORIENTATION portrait')
gmt('gmtset MAP_TITLE_OFFSET 36p') % avoids overlap with tick labels

% small font size for title, if filenames are included
if small_title % (else: leave GMT defaults)
    gmt('gmtset FONT_TITLE 12p')
end

end