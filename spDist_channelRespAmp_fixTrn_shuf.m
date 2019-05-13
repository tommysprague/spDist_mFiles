% spDist_channelRespAmp_fixTrn_shuf.m
% adapted from spDist_channelRespAmp_fixTrn.m and
% MGSMap_channelRespAmp_catSess_thruTime_fixTrn_shuf.m
%
% trains on all MGSMap runs found in trn_dir, using a fixed training time
% window, then reconstructs on all testing timepoints in tst_sess. voxel
% selection, if based on n_vox, is done using training data only
%
% computes reconstructions at each timepoint using a FIXED encoding model,
% estimated using specified timepoints on MGSMap data
%
% also computes using a shuffled model (shuffle trial labels for training
% set)
%
% input "sess" refers to TESTING sessions for now...
%
% TCS 5/8/2019
% tsprague@ucsb.edu
% github.com/tommysprague

       % for reference: 
       % (saved from spDist_scoreEyeData.m; via eyeData)
       
       % 1: distractor condition (1 = no, 2 = distractor)
       % 2,3: target position X, Y (dva)
       % 4,5: distractor position X,Y (or NaN; dva)
       % 6: distractor bin (-3:3; NaN); Cartesian, so + is CCW
       % 7: distractor direction (1 = CCW, 2 = CW, NaN = no dst)
       % 8: correct? (0 or 1; NaN)
       % 9: RT



function spDist_channelRespAmp_fixTrn_shuf(subj,sess,ROIs,which_vox,trn_tpts,n_shuf_iter)

tst_dir = 'spDist';
trn_dir = 'wmChoose';

trn_sess = 'MGSMap'; % files to load for training

root =  spDist_loadRoot;
trn_root = sprintf('%s/../wmChoose_scanner/',root);

if nargin < 1
    subj = {'AY','CC','EK','KD','MR','SF','XL'};
        
end
if nargin < 2
    sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}};
    
end

if nargin < 3
    ROIs = {'V1','V2','V3','V3AB','hV4','VO1','VO2','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS','iPCS'};
    
end


% analysis parameters:
n_chan = 8; % # of channels, evenly spaced around the screen
chan_centers = linspace(360/n_chan,360,n_chan);

% evaluate basis set at these
angs = linspace(-176,180,90);


% for debugging, also save out average reconstructions...
n_angs_sv = 30;
angs_sv = linspace(180+360/n_angs_sv,180,n_angs_sv);

if nargin < 4
    which_vox = 0.1; % top 1000 vox
end

if nargin < 5
    trn_tpts = [7:15]; % what we use to train model!
end


% # of times to shuffle training set (after computing intact model)
if nargin < 6
    n_shuf_iter = 1000;
end

align_to = {'targ_ang_all','dist_ang_all'};

func_suffix = 'surf';
delay_tpts = -3:26; % 0.8 s TR ---- what we want to reconstruct


% loop over subj, ROIs and load each session, concatenate, and process
for ss = 1:length(subj)
    
    for vv = 1:length(ROIs)
        
        % seed each ROI the same
        rng(spDist_randSeed);
        
        % load TESTING data from each session and concatenate
        data_tst = [];
        
        for sess_idx = 1:length(sess{ss})
            
            fn = sprintf('%s/%s_trialData/%s_%s_%s_%s_trialData.mat',root,tst_dir,subj{ss},sess{ss}{sess_idx},ROIs{vv},func_suffix);
            fprintf('loading TESTING data from %s...\n',fn);
            thisdata = load(fn);
            
            thisdata.sess = sess_idx*ones(size(thisdata.r_all));
            
            data_tst = cat_struct(data_tst,thisdata,{'rf','TR','which_TRs'}); % skip 'rf', these will be the same
            
        end
        
        
        % list training sessions for this ROI:
        
        tmp_fn_trn = dir( sprintf( '%s/%s_trialData/%s_%s*_%s_%s_trialData.mat', trn_root,trn_dir,subj{ss},trn_sess,ROIs{vv},func_suffix) );
        
        % loop over those and load (as above) - TRAINING data
        data_trn = [];
        fn_trn = cell(length(tmp_fn_trn),1);
        for ff = 1:length(tmp_fn_trn)
            fn_trn{ff} = sprintf('%s/%s',tmp_fn_trn(ff).folder,tmp_fn_trn(ff).name);
            fprintf('loading TRAINING data from %s...\n',fn_trn{ff});
            thisdata = load(fn_trn{ff});
            
            thisdata.sess = ff*ones(size(thisdata.r_map));
            
            data_trn = cat_struct(data_trn,thisdata,{'rf','TR','which_TRs'}); % skip 'rf', these will be the same
        end
        
        which_TRs_tst = data_tst.which_TRs;
        which_TRs_trn = data_trn.which_TRs;
        
        
        % because which_TRs doesn't necessarily start at 1...
        delay_idx = find(ismember(which_TRs_tst,delay_tpts));
        
        % for each item (target; distractor), n_trials x n_delay_tpts x
        % n_shuff_iter
        all_fidelity = cell(length(align_to),1);
        recons_sv = cell(length(align_to),1);
        for aa = 1:length(all_fidelity)
            all_fidelity{aa} = nan(size(data_tst.c_all,1),length(delay_tpts),n_shuf_iter);
            recons_sv{aa} = nan(size(data_tst.c_all,1),length(angs_sv),length(delay_tpts),n_shuf_iter);
        end

        
        
        
        % shuffle! ~~~~~~~~~~~~
        for shuf_iter = 1:n_shuf_iter
            
            
            % only save these for each shuffling iteration
            
            % save out recons rotated to align with target and with distractor
            % (nans for no-distractor trials)
            recons = cell(length(align_to),1);
            
            chan_resp = nan(size(data_tst.c_all,1),n_chan,length(delay_tpts));
            
            
            % voxel selection (from _genModelDecode2.m scripts) (NOTE: not
            % shuffled here!))
            
            
            % first - filter voxels based on whether there's signal (if no signal,
            % there are 0's the entire dataset it seems)
            trndata = mean(data_trn.dt_mapz(:,:,ismember(which_TRs_trn,trn_tpts)),3);
            mystd = std(trndata,[],1);
            
            % if which_vox < 1, that means we're using the RF VE to
            % constrain voxel choice (which will be the same for all CV
            % folds) - choose voxels which have data (no nan/0 std dev)
            % and VE >= threshold specified
            
            if which_vox < 1
                these_vox = mystd~=0 & ~isnan(mystd) & data_trn.rf.ve>=which_vox;
                
                % otherwise, we're using the top N voxels sorted by
                % quadrant-wise F-score, or all voxels, whichever is smaller
            else
                
                
                %trndata = trndata(:,mystd~=0 & ~isnan(mystd));
                
                %extra_vox = sum(mystd==0 | isnan(mystd)); % in case we have to remove voxels, add this many to F-socre
                
                % voxel selection (training data only)
                allF = nan(size(trndata,2),1);
                allp = nan(size(trndata,2),1);
                thisG = data_tst.c_map(trn_idx,2);
                for voxidx = 1:size(trndata,2)
                    
                    thisX = trndata(:,voxidx);
                    
                    
                    [p,tab,stats] = anova1(thisX,thisG,'off');
                    allF(voxidx) = tab{2,5};
                    allp(voxidx) = p;
                    clear thisX p tab stats;
                end
                
                % get rid of NaN's, which can in principle happen when
                % zero std dev, etc.
                f_sorted = sort(allF(~isnan(allF)),'descend');
                if which_vox <= length(allF) % handle case of small ROI
                    f_thresh = f_sorted(which_vox);
                else
                    f_thresh = f_sorted(end);
                end
                
                % deal with rare chance that there are some identical F
                % scores, which would allow incorrect # of vox
                these_vox = allF>=f_thresh;
                if sum(these_vox)>which_vox
                    allidx = find(these_vox);
                    these_vox(allidx((which_vox+1):end)) = 0;
                end
                
            end
                        
            trndata = trndata(:,these_vox);
            
            % n_trials x n_vox
            %trn = mean(data.dt_map(:,:,delay_idx),3);
            trn = trndata;
            
            %% set up design matrix; shuffle it
            
            X = spDist_makeX1(data_trn.c_map(:,1),chan_centers);
            
            
            X = X/max(X(:));
            
            X = X(randperm(size(X,1)),:); % shuffle!
            
            %% compute channel weights using only training set (wmMapSpace)
            
            %fprintf('computing channel weights\n');
            w = X\trn;
            
            % ~~~~~~~~ then, reconstruct each time point ~~~~~~~~~~
            
            for tpt_idx = 1:length(delay_tpts)
                %fprintf('TPT: %i\n',tpt_idx);
                
                
                tstdata = mean(data_tst.dt_allz(:,:,delay_idx(tpt_idx)),3);
                
                tstdata = tstdata(:,these_vox);
                tst = tstdata;
                
                %% use (optimized) design matrix to compute channel responses using testing data
                
                %fprintf('computing channel responses\n');
                
                chan_resp(:,:,tpt_idx) = (inv(w*w.')*w*tst.').';
                
                clear tst tstdata;
                
            end % end of testing/reconstruction loop 
            
            clear trn trn_idx tst_idx trn_runs tst_runs w trndata allF f_thresh thisG mystd;
            
            
            
            %% after all channel resp computed, coregister each trial
            
            
            for aa = 1:length(align_to)
                recons{aa} = nan(size(data_tst.c_all,1),length(angs),length(delay_tpts));
                for tt = 1:size(chan_resp,1)
                    
                    % we have to build a basis set for each trial, each target
                    
                    % remove the polar angle of the aligned target
                    rot_by = data_tst.(align_to{aa})(tt);%atan2d(data.xy_task(tt,2),data.xy_task(tt,1));
                    
                    if ~isnan(rot_by)
                        this_rfTh = chan_centers-rot_by; % rotate basis
                        myb = build_basis_polar_mat(angs,this_rfTh);
                        myb_save = build_basis_polar_mat(angs_sv,this_rfTh);
                    end
                    
                    % myb is length(angs) x n_channels
                    % we want to weight each channel (col) by this trials' channel
                    % activation
                    % (result shoudl be 1 x length(angs))
                    
                    
                    for tpt_idx = 1:length(delay_tpts)
                        % on no-distractor trials, we won't have to rotate by
                        % anything, so those won't have a distractor-aligned
                        % reconstruction
                        if ~isnan(rot_by)
                            recons{aa}(tt,:,tpt_idx) = (myb * chan_resp(tt,:,tpt_idx).').';
                            recons_sv{aa}(tt,:,tpt_idx,shuf_iter) = (myb_save * chan_resp(tt,:,tpt_idx).').';
                        end

                    end
                    
                    clear rot_by myb myb_save;
                end
                
                % compute fidelity
                all_fidelity{aa}(:,:,shuf_iter) = squeeze(mean(cosd(angs) .* recons{aa},2));
                
            end
        end
       
       % TCS 5/9/2019 - we won't do shuffling version of this now...
       
%        % average distractor representation on distractor trials
%        dist_mu = squeeze(mean( recons{2}(data_tst.c_all(:,1)==2,:,:), 1 )); % n_angs x n_tpts
%        
%        % variable for saving distractor-corrected single-item WM
%        % representations
%        recons_nodist = nan(size(recons{1})); 
%        
%        % angular difference - distractor relative to target (+ is CCW, -
%        % CW)
%        ang_diff = atan2d( sind(data_tst.dist_ang_all-data_tst.targ_ang_all), cosd(data_tst.dist_ang_all-data_tst.targ_ang_all) );
%        
%        zero_idx = find(angs==0); % 45
%        
%        % remove that average from each trial, after aligning based on
%        % relative distractor location per trial
%        for tt = 1:size(recons{1},1)
%            
%            % if distractor trial, remove average of all distractors
%            if data_tst.c_all(tt,1)==2
%                
%                % distance between distractor and WM position
%                % ang_diff(tt)
%                
%                % find nearest ang bin - how many units to move?
%                [~,ang_diff_idx] = min(abs(angs-ang_diff(tt)));
%                
%                
%                % circularly shift thismu by -1*(ang_diff_idx-zeroidx)
%                % (if positive, distractor is +ang compared to targ, which
%                % is right, so shift left)
%                recons_nodist(tt,:,:) = recons{1}(tt,:,:) - shiftdim(circshift(dist_mu,-1*(ang_diff_idx-zero_idx),1),-1);
%                
%            end
%        end
%        
       
       % things we want to save
        
        c_all = data_tst.c_all;
        r_all = data_tst.r_all;
        p_all = data_tst.p_all;
        a_all = [data_tst.targ_ang_all data_tst.dist_ang_all];
        sess_all  = data_tst.sess;
        
        
        % save with VE marker when which_vox < 1, otherwise, number of
        % vox
        if which_vox < 1
            fn2s = sprintf('%s/%s_reconstructions/%s_%s_%s_%s_%ichan_VE%03.f_trn%ito%i_recon_thruTime1_shuf%i.mat',root,tst_dir,subj{ss},horzcat(sess{ss}{:}),ROIs{vv},func_suffix,n_chan,100*which_vox,trn_tpts(1),trn_tpts(end),n_shuf_iter);
        else
            fn2s = sprintf('%s/%s_reconstructions/%s_%s_%s_%s_%ichan_%ivox_trn%ito%i_recon_thruTime1_shuf%i.mat'  ,root,tst_dir,subj{ss},horzcat(sess{ss}{:}),ROIs{vv},func_suffix,n_chan,    which_vox,trn_tpts(1),trn_tpts(end),n_shuf_iter);
        end
        fprintf('saving to %s...\n',fn2s);
        
        save(fn2s,'c_all','r_all','p_all','n_chan','delay_tpts','angs','all_fidelity','recons_sv','which_vox','sess_all','these_vox','a_all','fn_trn');
        
        clear data c_all r_all p_all chan_resp w_all recons a_all all_fidelity;
        
        
        
    end
    
    
end
return