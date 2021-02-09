function gmtMexGrid2grd(GridStruct, outFilename, varargin)
%gmtMexGrid2grd write a gmt mex grid structure to a grid file
%
%   Input arguments:
%      - GridStruct : gtm mex grid structure
%      - outFilename : path and filename of output file
%      - optional, outFormat : format string (see gmt grdconvert documentation)
%                              defaults to 'nf', GMT netCDF format (32-bit float)
%                              which has wide compatibility (Surfer included)
%                              note: to output a LithoFlex compatible ascii grd
%                                    use 'gd:GSAG' (needs GDAL)
%
%   Output arguments:
%      none
%
% 2021-01-13 AP

narginchk(2,3)
nargoutchk(0,0)

if nargin==3
    outFormat=varargin{1};
else
    outFormat='nf';
end

disp(['Writing grid to: ', outFilename])
gmt('write', ['-Tg ', outFilename, '=', outFormat], GridStruct)

end

