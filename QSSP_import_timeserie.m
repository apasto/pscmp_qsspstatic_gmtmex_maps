function timeserie = QSSP_import_timeserie(filename)
startRow = 2;
endRow = inf;
formatSpec = '%12f%f%[^\n\r]';

fileID = fopen(filename, 'r');
dataArray = textscan( ...
    fileID, formatSpec, endRow(1)-startRow(1)+1, ...
    'Delimiter', '', 'WhiteSpace', '', ...
    'TextType', 'string', 'EmptyValue', NaN, ...
    'HeaderLines', startRow(1)-1, 'ReturnOnError', false, ...
    'EndOfLine', '\r\n');

fclose(fileID);

timeserie = table(dataArray{1:end-1}, 'VariableNames', {'day', 'data'});
end
