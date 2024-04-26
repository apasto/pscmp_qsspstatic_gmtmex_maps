function grd2gmtMap_plot_map(...
    out_filename, map_grid, RJstring, cpt,...
    unit, tick_int_minor, tick_int_major, map_title,...
    plot_countries, plot_contours,...
    source_type, source_lon_lat,...
    lon_lat_grid_interval, varargin)
% grd2gmtMap_plot_map
narginchk(13, 19)

if nargin>13 && ~isempty(varargin{1})
    countries_color = varargin{1};
else
    countries_color = 'White';
end
if nargin>14 && ~isempty(varargin{2})
    coasts_color = varargin{2};
else
    coasts_color = 'White';
end
if nargin>15 && ~isempty(varargin{3})
    scale_width = varargin{3};
else
    scale_width = '5.0i';
end
if nargin>16 && ~isempty(varargin{4})
    scale_triangles = varargin{4};
else
    scale_triangles = '+e';
end
% opt argument: vertical shift of colorscale (needs to be finetuned)
if nargin>17 && ~isempty(varargin{5})
    scale_Yc_shift = varargin{5};
else
    scale_Yc_shift = '-0.5i';
end
% opt argument: vertical shift of min/max text above colorscale
if nargin>18 && ~isempty(varargin{6})
    stats_Y_shift = varargin{6};
else
    stats_Y_shift = '-3.60i';
end

if isempty(tick_int_major)
    tick_int_major = tick_int_minor * 5;
end

% save arguments (except grid) to structure
% for reproducibility
plot_map_args.RJstring = RJstring;
plot_map_args.cpt = cpt;
plot_map_args.unit = unit;
plot_map_args.tick_int_minor = tick_int_minor;
plot_map_args.tick_int_major = tick_int_major;
plot_map_args.map_title = map_title;
plot_map_args.plot_countries = plot_countries;
plot_map_args.plot_contours = plot_contours;
plot_map_args.source_type = source_type;
plot_map_args.source_lon_lat = source_lon_lat;
plot_map_args.lon_lat_grid_interval = lon_lat_grid_interval;
plot_map_args.countries_color = countries_color;
plot_map_args.coasts_color = coasts_color;
plot_map_args.scale_width = scale_width;
plot_map_args.scale_triangles = scale_triangles;
plot_map_args.scale_Yc_shift = scale_Yc_shift;
plot_map_args.stats_Y_shift = stats_Y_shift;
[out_path, out_name, ~] = fileparts(out_filename);
plot_map_args_filename = [out_path, '/', out_name, '_plot_map_args.mat'];
save(plot_map_args_filename, 'plot_map_args');

% grdimage and title
gmt(['grdimage ', RJstring,...
    ' -E -Q -nc+c -C -K -P -Xc -Yc -B+t"', map_title, '" > ',...
    out_filename],...
    map_grid, cpt);

% countries
if plot_countries
    gmt(['pscoast -R -J',...
        ' -P -N1/0.01c,' countries_color,...
        ' -O -K -Xc -Yc >> ',...
        out_filename])
end
% coasts
gmt(['pscoast -R -J',...
    ' -P -Di+ -A20000',...
    ' -W0.03c,', coasts_color,' -O -K -Xc -Yc >> ',...
    out_filename])

% source point/rectangle/rectangles
switch source_type % source dot or fault area rectangle (hor projection)
    case 'point'
        gmt(['psxy -R -J',...
            ' -Sc10p -Gpurple',...
            ' -O -K -Xc -Yc >> ',...
            out_filename], source_lon_lat);
    case 'rectangle'
        if ~iscell(source_lon_lat) % one rectangle
            gmt(['psxy -R -J',...
                ' -A -L -Wthick -Wpurple',... % -A straight-line segments
                ' -O -K -Xc -Yc >> ',...
                out_filename], [source_lon_lat; source_lon_lat(1,:)]); % repeat first point, close rectangle
        else
            for n=1:size(source_lon_lat, 2)
                gmt(['psxy -R -J',...
                    ' -A -L -Wthick -Wpurple',... % -A straight-line segments
                    ' -O -K -Xc -Yc >> ',...
                    out_filename], [source_lon_lat{n}; source_lon_lat{n}(1,:)]); % repeat first point, close rectangle
            end
        end
end

% contours
if plot_contours
    gmt(['grdcontour -R -J',...
        ' -C',num2str(tick_int_minor)...
        ' -A',num2str(tick_int_major),'+f13p'...
        ' -S',num2str(8),...
        ' -Gn1+r15p',... % +r: enforce minimum separation radius bw labels
        ' -t50 -O -K -P -Xc -Yc -Wathin -Wcfaint >> ',...
        out_filename],...
        map_grid);
end

% semi-trasparent grid
gmt(['psbasemap -R -J',...
    ' -B', num2str(lon_lat_grid_interval), 'g' , num2str(lon_lat_grid_interval),... % e.g. -B1g1
    ' -t10 -BWesN -O -K -Xc -Yc >> ',...
    out_filename]);

% non-trasparent frame and ticks
gmt(['psbasemap -R -J',...
    ' -B', num2str(lon_lat_grid_interval), ' -BWesN -O -K -Xc -Yc >> ',...
    out_filename]);

% min/max range text
stats_font_size = '12p';
stats_X_shift = '0i';
stats_str_min_format = '%6.3f';
stats_str_max_format = '%6.3f';
% scientific notation if min or max are very small
if abs(map_grid.range(5)) < 1e-3
    stats_str_min_format = '%6.3e';
end
if abs(map_grid.range(6)) < 1e-3
    stats_str_max_format = '%6.3e';
end
stats_text = [...
    'min = ',num2str(map_grid.range(5), stats_str_min_format),...
    ', max = ',num2str(map_grid.range(6), stats_str_max_format)];

gmt(['pstext ',...
    ' -R1/10/1/10 -JX10 -F+f',stats_font_size,'+cCM',... % +cCM : justification
    ' -P -O -K -Xc',stats_X_shift,' -Yc',stats_Y_shift,...
    ' >> ',out_filename], stats_text)

% scale
if isempty(unit)
    by_args = ' -By';
else
    by_args = ' -By+l';
end
gmt(['psscale ',...
    '-DJBC+o0/0.7i+w', scale_width, scale_triangles,' -R -JX', scale_width, ' -O -P ',...
    ' -Bx', num2str(tick_int_major),'f', num2str(tick_int_minor),...
    by_args, unit, ' -Xc -Yc', scale_Yc_shift, ' >> ',...
    out_filename], cpt)

end
