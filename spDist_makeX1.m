% based on wmLazOri_makeX1.m
%
% TCS 9/18/2015
%
% adapted loosely from wmDrop_makeX1 - TCS 7/9/2015
% - instead of a pixel mask, just uses set of inline functions to compute
%   response of each filter at each orientation
% trnX = wmLazOri_makeX1(stimAngs,rfTh);


function trnX = spDist_makeX1(stimAngs,rfTh)

my_basis_fcn = @(th,mth) build_basis_polar_mat(th,mth);%cosd(0.5*min(mod(th-mth,360),mod(mth-th,360))).^(length(rfTh)-mod(length(rfTh),2));


trnX = nan(length(stimAngs),length(rfTh));

% stimOris in deg, rfTh in deg
for ii = 1:length(stimAngs)
    %for jj = 1:length(rfTh)
        trnX(ii,:) = my_basis_fcn(stimAngs(ii),rfTh);
    %end
end


return