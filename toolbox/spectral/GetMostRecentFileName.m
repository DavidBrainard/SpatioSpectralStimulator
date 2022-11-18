function [fileName] = GetMostRecentFileName(testFileDir,testFileName,options)
% This function finds the file name with the most recent date among the
% files having the same name.
%
% Syntax:
%    [fileName] = GetMostRecentFileName(testFileDir,testFileName)
%
% Description:
%    This function finds the most up-to-date file name among the files with
%    the same name. This will be useful when we want to read the most
%    recent data file when using SpectralCalCheck and SpectralCalAnalyze.
%
% Inputs:
%    testFileDir -                Directory where the files are saved.
%    testFileName -               File name that appears in common for
%                                 multiple files.
%
% Outputs:
%    fileName -                   The file name created on the most recent
%                                 date.
%
% Optional key/value pairs:
%    'olderDate' -                Set this number to get older date than
%                                 the most recent one. If you set this to
%                                 1, the result will be the second to the
%                                 most recent date. Likewise, setting this
%                                 to 2 finds third to the most recent date.

% History:
%    11/10/21  smo                Started on it

%% Set parameters.
arguments
    testFileDir
    testFileName
    options.olderDate (1,1) = 0
end

%% Find the most recent date.
%
% Find the files with the same name.
listOfFiles = dir(append(fullfile(testFileDir,testFileName),'*'));

% Find the most recent date here.
%
% 'listOfFiles' will be the struct variable with the size of nx1 (n =
% number of files)
for ii = 1:size(listOfFiles,1)
    fileDateTemp = listOfFiles(ii).date;
    listOfDatenums(ii) = datenum(fileDateTemp);
end
[maxDateNum maxDateOrder] = max(listOfDatenums);

% Save out the file name.
% 
% If you want to save out older file than the most recent one, you can
% change the number in 'options.olderDate'
fileNameOnly = listOfFiles(maxDateOrder - options.olderDate).name;
fileName     = fullfile(testFileDir,fileNameOnly);

end