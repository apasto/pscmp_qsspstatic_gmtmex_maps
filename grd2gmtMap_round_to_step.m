function out = grd2gmtMap_round_to_step(step,in,direction)
%grd2gmtMap_round_to_step round(/floor/ceil/fix) input vector to nearest arbitrary 'step' unit
% digit-counting code comes from this answer by user Jaymin on matlabcentral:
% https://mathworks.com/matlabcentral/answers/10795-counting-the-number-of-digits#answer_68482
%
% Syntax: out = RoundToStep(step,in,[direction])
%
% Inputs:
%    step        : scalar, unit value
%    in          : vector of values to be rounded
%    [direction] : char vector, optional, direction of rounding
%                  'round', 'floor', 'ceil', 'fix'
%                  'round' is default behaviour
%                  case insensitive
%
% Outputs:
%    out         : rounded 'in' vector
%


% check and manage input
narginchk(2,3)
if nargin==2
    direction='round';
end
assert(isscalar(step),'Argument ''step'' must be a scalar.')
assert(isvector(in),'Argument ''in'' must be a vector.')
% check if 'direction' is an allowed round-like function
allowed_directions = {'round','floor','ceil','fix'};
assert(any(strcmpi(allowed_directions,direction)),...
    ['''',direction,''' is not an allowed rounding function.'])

% count the number of digits
stepDIG = numel(num2str(step))-numel(strfind(num2str(step),'-'))-numel(strfind(num2str(step),'.'));

% number of zsteps in next power of 10
stepdiv = (10^(stepDIG))/step;

% round/floor/ceil/fix to step
roundfunc = str2func(lower(direction)); % create handle to selected rounding function
out = ((roundfunc((in/(10^(stepDIG)))*stepdiv))/stepdiv)*(10^(stepDIG));

end