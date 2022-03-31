function grd2gmtMap_plot_reproduce(...
    grid_filename, plot_map_arguments_filename,...
    out_filename, new_cpt, z_factor, custom_RJstring,...
    varargin)
%grd2gmtMap_plot_reproduce plot a grd using a previously plotted map arguments
% TODO: varargin, optional arguments

plot_map_arguments = load(plot_map_arguments_filename, 'plot_map_args');
plot_map_arguments = plot_map_arguments.plot_map_args;
% TODO: assert type and contents (fieldnames) of plot_map_arguments

grd2gmtMap_set_gmt_defaults(true) % use small title font, is not saved in plot_map_arguments

% load grid
in_grid = gmt('read', ['-Tg ', grid_filename]);
if ~isempty(z_factor)
    in_grid.z = in_grid.z * z_factor;
    in_grid.range(5) = min(in_grid.z(:));
    in_grid.range(6) = max(in_grid.z(:));
end

if ~isempty(custom_RJstring)
    RJstring = custom_RJstring;
else
    RJstring = plot_map_arguments.RJstring;
end

% opt argument: override grid step
if nargin>6 && ~isempty(varargin{1})
    plot_map_arguments.lon_lat_grid_interval = varargin{1};
end
% opt argument: override scale_Yc_shift
if nargin>7 && ~isempty(varargin{2})
    plot_map_arguments.scale_Yc_shift = varargin{2};
end
% opt argument: override stats_Y_shift
if nargin>8 && ~isempty(varargin{3})
    plot_map_arguments.stats_Y_shift = varargin{3};
end
% opt argument: use fixed min/max colorscale range
if nargin>9 && ~isempty(varargin{4})
    colorscale_min_max_are_fixed = true;
    colorscale_min_max = varargin{4}(:);
    assert(all(size(colorscale_min_max) == [4, 1]), ...
        'Wrong size of provided colorscale min/max/minor_int/major_int, must be a 4 element vector');
else
    colorscale_min_max_are_fixed = false;
    colorscale_min_max = [];
end
% opt argument: override source type and lon lat (point, rectangles)
if nargin==12 && ~isempty(varargin{5}) && ~isempty(varargin{6})
    plot_map_arguments.source_type = varargin{5};
    plot_map_arguments.source_lon_lat = varargin{6};
end

if new_cpt
    % temp: use fixed colorscale
    cpt_colorscale = 'roma -E -Z';
    if ~colorscale_min_max_are_fixed
        % temp: use fixed min and max quantiles
        min_quantile = 0.001; % default
        max_quantile = 0.999; % default
        [...
            tick_int_minor,...
            tick_int_major,...
            cpt] = grd2gmtMap_call_grd2cpt(...
            in_grid,...
            cpt_colorscale,...
            min_quantile, max_quantile);
    else
        tick_int_minor = colorscale_min_max(3);
        tick_int_major = colorscale_min_max(4);
        cpt = gmt([...
        'grd2cpt -C', cpt_colorscale,...
        ' -L',...
        num2str(colorscale_min_max(1)),'/',...
        num2str(colorscale_min_max(2)),...
        ' -Di'], in_grid);
    end
else
    cpt = plot_map_arguments.cpt;
    tick_int_minor = plot_map_arguments.tick_int_minor;
    tick_int_major =  plot_map_arguments.tick_int_major;
end

grd2gmtMap_plot_map(...
            out_filename, in_grid,...
            RJstring, cpt,...
            plot_map_arguments.unit, tick_int_minor, tick_int_major,...
            plot_map_arguments.map_title,...
            plot_map_arguments.plot_countries, plot_map_arguments.plot_contours,...
            plot_map_arguments.source_type, plot_map_arguments.source_lon_lat,...
            plot_map_arguments.lon_lat_grid_interval,...
            plot_map_arguments.countries_color,...
            plot_map_arguments.coasts_color,...
            plot_map_arguments.scale_width,...
            plot_map_arguments.scale_triangles,...
            plot_map_arguments.scale_Yc_shift,...
            plot_map_arguments.stats_Y_shift);

[out_path, out_name, ~] = fileparts(out_filename);
[~, ~] = grd2gmtMap_psconvert_gs_fallback(...
    out_filename, './', [out_path, out_name, '_conv.png']);

end
