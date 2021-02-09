function gravityOut = countercorrectFreeAir_pscmp_qssp(gravityIn, displacementIn, varargin)
%countercorrectFreeAir_pscmp_qssp remove FA correction from pscmp/qssp gravity
% psgrn/pscmp and qssp apply a free-air correction to model
% the movement of a gravimeter placed on the ground (which is deformed)
% to use the gravity data for fixed height modelling (e.g. sat altitude)
% we need to REMOVE the appplied correction.
%
%   Input arguments:
%      - gravityIn : gravity values as obtained by pscmp or qssp [m/s^2]
%      - displacementIn : upwards-positive displacement [m]
%                         note that pscmp provides downwards positive displacement
%                         while qssp displacement U_z is upwards positive
%      - optional, FreeAirGradient : provide a custom Free Air Gradient
%                                    default: approx -3.0827e-06 1/s^2
%                                    as used by psgrn / qssp
%
%   Output arguments:
%      - gravityOut : de-corrected gravity values
%
% 2021-01-23 AP

narginchk(2,3)
nargoutchk(0,1)

% optional FreeAirGradient argument
if nargin>2
    assert(isscalar(varargin{1}), 'provided Free Air Gradient is not a scalar')
    FreeAirGradient = varargin{1};
else
    % use Free Air Gradient as in psgrn and qssp original code
    FreeAirGradient = -(9.82 * 2 / 6371e3); %-g*2/R [1/s^2]
end

assert(isequal(size(gravityIn), size(displacementIn)),...
    'Size of provided gravity and displacement is not equal')

% we have to counter-correct the applied correction
% therefore we SUBTRACT it: points that have risen get a POSITIVE counter-correction
% => note the minus sign before FreeAirGradient
%   => positive displacement -> positive counter-correction

gravityOut = gravityIn + (displacementIn .* (-FreeAirGradient));

end

