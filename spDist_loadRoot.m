% spDist_loadRoot.m
%
% returns the root directory - easier for distributing code, etc

function root = spDist_loadRoot

rootdir = 'spDist';

if ismac
    root = sprintf('/Volumes/data/%s/',rootdir);
else
    root = sprintf('/deathstar/data/%s/',rootdir);
end

return