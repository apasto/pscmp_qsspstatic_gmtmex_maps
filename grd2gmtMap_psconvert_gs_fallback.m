function [gsCallStatus, gsCallCmdOut] = grd2gmtMap_psconvert_gs_fallback(...
    ps_filename, out_path, out_filename, varargin)
%grd2gmtMap_psconvert_gs_fallback instead of psconvert, call gs

narginchk(3, 4)
if nargin>3 && ~isempty(varargin{1})
    assert((ischar(varargin{1}) || isstring(varargin{1})),...
        'gsPath argument must be type char or string.');
    gs_path = varargin{1};
else
    % TODO: default should be platform dependent
    gs_path = 'gswin64c';
end

[gsCallStatus, gsCallCmdOut] = system([...
    gs_path, ' -dNOSAFER -dNOPAUSE -dBATCH -DPSL_no_pagefill ',...
    '-sDEVICE=png16m -dUseCropBox -r200 -sOutputFile=',...
    out_path, out_filename, '.png ', ps_filename]);

end
