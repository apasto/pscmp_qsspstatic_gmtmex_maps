function [minor_tick_int, major_tick_int, out_cpt] = grd2gmtMap_call_grd2cpt(...
    in_grid, cpt_color_scale, varargin)
%grd2gmtMap_call_grd2cpt create a suitable colorscale for a grid structure
%   Provided a gmt-mex grid structure in_grid and a cpt string (-C option of grd2cpt)
%   attempt to create an 'optimal' cpt
%   example cpt_color_scale argument: 'roma -E -Z'

narginchk(2, 4)
nargoutchk(3, 3)

% optional arguments: quantiles
% min_quantile
if nargin>2 && ~isempty(varargin{1})
    min_quantile = varargin{1};
else
    min_quantile = 0;
end
% max_quantile
if nargin>3 && ~isempty(varargin{2})
    max_quantile = varargin{2};
else
    max_quantile = 1;
end

% skip zero-z-range grid, to avoid grd2cpt errors
if ~((in_grid.range(6) - in_grid.range(5)) == 0)
    % compute 'optimal' minor tick interval (TO DO: could be better!)
    minor_tick_int = ...
        grd2gmtMap_round_to_step(0.1, ...
            (in_grid.range(6) - ...
            in_grid.range(5)) / 25,...
            'ceil');
    if minor_tick_int > 10
        minor_tick_int = grd2gmtMap_round_to_step(...
            5, minor_tick_int, 'ceil');
    end
    major_tick_int = minor_tick_int*5;
    if major_tick_int > 5
        major_tick_int = grd2gmtMap_round_to_step(...
            5, major_tick_int, 'ceil');
        minor_tick_int = major_tick_int / 5;
    end
    out_cpt = gmt([...
        'grd2cpt -C', cpt_color_scale,...
        ' -L',...
        num2str(grd2gmtMap_round_to_step(...
            minor_tick_int,...
            quantile(in_grid.z(:),min_quantile),...
            'floor')),'/',...
        num2str(grd2gmtMap_round_to_step(...
            minor_tick_int,...
            quantile(in_grid.z(:),max_quantile),...
            'ceil')),...
        ' -Di'],...
        in_grid);
    % avoid case in which the colorscale covers less than a 2 major ticks range
    l = 0; % re-set loop counter
    while (out_cpt.minmax(2) - out_cpt.minmax(1)) < (major_tick_int * 2)
        major_tick_int = minor_tick_int / 2;
        minor_tick_int = minor_tick_int / 10;
        % re-call grd2cpt
        out_cpt = gmt([...
            'grd2cpt -C', cpt_color_scale,...
            ' -L',...
            num2str(grd2gmtMap_round_to_step(...
            minor_tick_int,...
            quantile(in_grid.z(:),min_quantile),...
            'floor')),'/',...
            num2str(grd2gmtMap_round_to_step(...
            minor_tick_int,...
            quantile(in_grid.z(:),max_quantile),...
            'ceil')),...
            ' -Di'],...
            in_grid);
        l = l+1;
        if l>4
            break
        end
    end
    % avoid case in which there are 'too many' divisions in the cpt range
    l = 0; % re-set loop counter
    while ((out_cpt.minmax(2) - out_cpt.minmax(1)) / minor_tick_int) > 35
        major_tick_int = major_tick_int * 2;
        minor_tick_int = minor_tick_int * 2;
        l = l+1;
        if l>4
            break
        end
    end
else % grid has zero z-range
    minor_tick_int = NaN;
    major_tick_int = NaN;
    out_cpt = NaN;
end
end