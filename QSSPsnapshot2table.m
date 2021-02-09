function out_table = QSSPsnapshot2table(filename)
%QSSPsnapshot2table import a 'snapshot' form filename obtained with PSCMP
% snapshot : see 'OUTPUTS' section of .inp file for QSSP
%            snapshot files are output files which, if requested,
%            contain all requested computed data at a prescribed time
%               (for comparison: 'normal' non-snapshot files contain
%               ONE observable, points=cols, times=rows)
% WARNING:
%    no test has been implemented (yet) to assert if provided filename IS a snapshot
%
% 2021-01-23 AP

narginchk(1,1)
nargoutchk(0,1)

dataLines = [2, Inf];

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 11);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["Latdeg", "Londeg", "Station", "U_n", "U_e", "U_z", "Vstrain", "Grav", "Geoid", "Tilt_n", "Tilt_e"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Import the data
out_table = readtable(filename, opts);

end