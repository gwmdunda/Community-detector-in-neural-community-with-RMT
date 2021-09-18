baseFileName = '140203_Slice_2_SMOOTHED_TRACES.xlsx';
folder = 'C:\Users\GERRY\Documents\Canvas\ELEC 5450\Project dataset and code\Photoperiod-cellular-data-traces-master';
fullFileName = fullfile(folder, baseFileName);
if exist(fullFileName, 'file')
  % File exists.  Read it into a table.
  ts = readtable(fullFileName);
  ts = table2array(ts);
else
  % File does not exist.  Warn the user.
  errorMessage = sprintf('Error: file not found:\n\n%s', fullFileName);
  uiwait(errordlg(errorMessage));
  return;
end