function tSp = readtfile2(tFilePath,warnOnMissing)
% READTFILE   Read t-file and return vector of raw spike times
%
% tSp = READTFILE(TFILEPATH) reads the t-file specified by the full or
% relative path TFILEPATH and returns a vector of spike times in raw units
% (one-tenth of a millisecond, November 2014).
%
% tSp = READTFILE(...,WARNONMISSING), where WARNONMISSING is a logical
% scalar, issues a warning if the specified t-file cannot be found or
% opened. WARNONMISSING is set to TRUE by default (cf. readtlist) since it
% expected that this function will only be run a single argument from the
% command line, where the user is expecting non-empty output. If the t-file
% cannot be found at the specified location, an empty vector is returned.
%
% Edited by: Christopher Cannova
%
% See also TFILEINFO, READTLIST, GETSPIKETIMES, TLISTINFO
if nargin<2
    warnOnMissing = 1;
end
[fID msg] = fopen(tFilePath, 'rb','b');
if fID == -1
    if warnOnMissing
        warning([ 'Could not open tfile ' tFilePath ': ' msg]);
    end
    tSp = [];
else
    ReadHeader(fID); % first read the header, then read a tfile
    tSp = fread(fID,inf,'uint32');	%read as 32 bit ints
    fclose(fID);
end
tSp = tSp*10^-4;
end