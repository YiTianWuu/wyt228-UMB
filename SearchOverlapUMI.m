function [overlaplist notoverlaplist] = SearchOverlapUMI(outfile1, outfile2)
% Find UMI overlap of MUTATIONS between outfile1 and outfile2
% outfile1 is from the source library of contamination (i.e. with high VAF)
% outfile2 is from the contaminated library

% Find mutation UMIs in outfile2 (don't check wt molecules)
data2 = readtable(outfile2);
searchlist = data2;

toDelete = strcmp(searchlist.Var3,'wt');
searchlist(toDelete,:) = [];

% Read outfile1 and search mutation UMIs from outfile2
filetext1 = fileread(outfile1);

toKeep = zeros(size(searchlist,1),1); % indices of overlap UMIs
for curline = 1:size(searchlist,1)
    if ~isempty(strfind(filetext1,searchlist{curline,1})) % can find current UMI in outfile1
        toKeep(curline) = 1;
    end
end

overlaplist = searchlist;
notoverlaplist = searchlist;

toKeep = logical(toKeep);
overlaplist(~toKeep,:) = [];
notoverlaplist(toKeep,:) = [];



