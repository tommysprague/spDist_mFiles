% spDist_loadRoot.m
%
% returns the root directory - easier for distributing code, etc

function root = spDist_loadRoot

rootdir = 'spDist';

if ismac
    root = sprintf('/Volumes/data/%s/',rootdir);
else
    this_computer = char(java.net.InetAddress.getLocalHost.getHostName);
    if strcmpi(this_computer,'tcs-compute-1')
        root = '/mnt/LabShare/projects/nyu/spDist/';
    else % for vader.psych.nyu.edu
        root = '/deathstar/data/spDist/';
    end
end

return