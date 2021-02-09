function out_table = PSCMPsnapshot2table(filename)
%PSCMPsnapshot2table import a 'snapshot' form filename obtained with PSCMP
% snapshot : see 'OUTPUTS' section of .inp file for PSCMP
%            snapshot files are output files which, if requested,
%            contain all requested computed data at a prescribed time
%               (for comparison: 'normal' non-snapshot files contain
%               ONE observable, points=cols, times=rows)
% WARNING:
%    no test has been implemented (yet) to assert if provided filename IS a snapshot
%
% 2020-01-07 AP

narginchk(1,1)
nargoutchk(0,1)

startRow = 2;
endRow = inf;
formatSpec = '%13f%13f%13f%13f%13f%13f%13f%13f%13f%13f%13f%13f%13f%13f%13f%13f%f%[^\n\r]';

fileID = fopen(filename,'r');

% Read columns of data according to the format.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

fclose(fileID);

out_table = table(...
    dataArray{1:end-1},...
    'VariableNames',...
    {'Latdeg','Londeg',...
    'Disp_north','Disp_east','Disp_down',...
    'Stress_nn','Stress_ee','Stress_dd',...
    'Stress_ne','Stress_ed','Stress_dn',...
    'Tilt_n','Tilt_e','Rotation',...
    'Geoid','Gravity','Disp_LOS'});
