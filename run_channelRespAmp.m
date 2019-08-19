% run_channelRespAmp.m
%
% UPDATED 2/23/2018 - now loops over a list that's generated containing all
% analyses to run, and assigns them in whatever order parfor decides.
% should be ~2x as quick, if not faster (no bottleneck on some subj, etc)

subj = {'AY','CC','EK','KD','MR','XL','SF'};
sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}}; % the set of sessions to be run on EVERYONE in subj...
ROIs = {'V1','V2','V3','V3AB','hV4','VO1','VO2','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS','iPCS'};


%ROIs = {'V1','V2','V3','V3AB','hV4','VO1','VO2','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS','iPCS','V1V2V3','IPS0IPS1','IPS2IPS3'};

which_vox = 0.1; % 10% VE thresh
shuf_iter = 1000;


% which analyses to run on everyone? 
% - 1 = spDist_channelRespAmp_fixTrn.m and shuffled

which_analyses = [1]; % specify which analyses to run on eveyrone (MGSMap)

to_run = nan(length(subj)*length(ROIs)*length(which_analyses),4);


thisidx = 1;
for ss = 1:length(subj)
    
    for rr = 1:length(ROIs)
        for aa = 1:length(which_analyses)
            to_run(thisidx,:) = [ss ss rr which_analyses(aa)];
            thisidx = thisidx+1;
        end
    end
    
end

% to_run is:
% 1) SUBJ
% 2) SESS 
% 3) ROI
% 4) ANALYSIS
ncores = feature('numcores');
mypool = parpool(ncores); % use all possible cores

parfor ii = 1:size(to_run,1)
    
    if to_run(ii,4) == 1
        spDist_channelRespAmp_fixTrn(     {subj{to_run(ii,1)}},{sess{to_run(ii,2)}},{ROIs{to_run(ii,3)}},which_vox,7:15)
        spDist_channelRespAmp_fixTrn_shuf({subj{to_run(ii,1)}},{sess{to_run(ii,2)}},{ROIs{to_run(ii,3)}},which_vox,7:15,shuf_iter)
    end
    
end

delete(mypool);
