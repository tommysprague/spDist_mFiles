% spDist_channelRespAmp_GATdist.m
% adapted from MGSMap_channelRespAmp_catSess_GAT1.m
%
% trains using distractor trials, both on target position and distractor
% position, at each tpt in turn, then reconstructs and computes fidelity at
% each timepoint
%
%
% TCS 8/19/2019
%
function spDist_channelRespAmp_GATdist(subj,sess,ROIs,which_vox,trn_tpts)

tst_dir = 'spDist';


root =  spDist_loadRoot;

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
if nargin < 4
    which_vox = 0.1; % top 1000 vox
end

% if nargin < 5
%     trn_tpts = [7:15]; % what we use to train model!
% end

align_to = {'targ_ang_all','dist_ang_all'};

func_suffix = 'surf';
delay_tpts = -3:26; % 0.8 s TR ---- what we want to reconstruct


% loop over subj, ROIs and load each session, concatenate, and process
for ss = 1:length(subj)
    
    for vv = 1:length(ROIs)
        
        % load TESTING (and training for GAT...) data from each session and concatenate
        data = [];
        
        for sess_idx = 1:length(sess{ss})
            
            fn = sprintf('%s/%s_trialData/%s_%s_%s_%s_trialData.mat',root,tst_dir,subj{ss},sess{ss}{sess_idx},ROIs{vv},func_suffix);
            fprintf('loading TESTING data from %s...\n',fn);
            thisdata = load(fn);
            
            thisdata.sess = sess_idx*ones(size(thisdata.r_all));
            
            data = cat_struct(data,thisdata,{'rf','TR','which_TRs'}); % skip 'rf', these will be the same
            
        end
        
        
        which_TRs = data.which_TRs;
        %which_TRs_tst = data_tst.which_TRs;
        %which_TRs_trn = data_trn.which_TRs;
        
        % write down an index we can use for LORO - 3 digits, first is
        % sessidx, then run index w/in each session (just for convenience)
        data.r_LORO = 100*data.sess+data.r_all;
        
        
        % because which_TRs doesn't necessarily start at 1...
        delay_idx = find(ismember(which_TRs,delay_tpts));
        %IEM_trn_tpt_idx = find(ismbember(which_TRs,trn_tpts)); % tpts to average over when training IEM
        
        % save out recons rotated to align with target and with distractor
        % (nans for no-distractor trials)
        % n_trn_tpts x n_tst_tpts x 2 (train w/ target location, train w/
        % distractor location)
        recons = cell(length(delay_tpts),length(delay_tpts), 2); 
        recons_raw = cell(length(delay_tpts),length(delay_tpts), 2); 
        chan_resp = cell(length(delay_tpts),length(delay_tpts), 2);
        
        % NOTE: to keep things simple, we'll fill in no-distractor trials
        % w/ NaN above
        
        
        
        this_ru = unique(data.r_LORO); % all the runs we'll CV over
        n_folds = length(this_ru);
        
        for trn_tpt_idx = 1:length(delay_tpts)
            for tst_tpt_idx = 1:length(delay_tpts)
                
                % use align_to variable for training/testing?
                for aa = 1:length(align_to)
                    
                    fprintf('Training TPT: %i, Testing TPT: %i\n',trn_tpt_idx,tst_tpt_idx);
                    
                    chan_resp{trn_tpt_idx,tst_tpt_idx,aa}  = nan(size(data.c_all,1),n_chan);
                    recons{trn_tpt_idx,tst_tpt_idx,aa}     = nan(size(data.c_all,1),length(angs));
                    recons_raw{trn_tpt_idx,tst_tpt_idx,aa} = nan(size(data.c_all,1),length(angs));
                    
                    
                    
                    for fold_idx = 1:n_folds
                        
                        % ~~~~~~~ first, estimate IEM ~~~~~~~~~~
                        
                        
                        
                        
                        % pick CV sets
                        trn_runs = ones(length(unique(data.r_LORO)),1);
                        
                        
                        % only use distractor-present trials for
                        % training/testing - selecting here should keep all
                        % indices correct throughout...
                        trn_idx = data.r_LORO~=this_ru(fold_idx) & data.c_all(:,1)==2;
                        tst_idx = data.r_LORO==this_ru(fold_idx) & data.c_all(:,1)==2;
                        
                        
                        % select voxels, data
                        trndata = mean(data.dt_allz(trn_idx,:,delay_idx(trn_tpt_idx)),3);
                        mystd = std(trndata,[],1);
                        
                        tstdata = mean(data.dt_allz(tst_idx,:,delay_idx(tst_tpt_idx)),3);
                        
                        % if which_vox < 1, means we're using VE threshold
                        if which_vox < 1
                            these_vox = mystd~=0 & ~isnan(mystd) & data.rf.ve >= which_vox;
                            
                            % otherwise, rank-order voxels by quadrant-wise
                            % F-score and select top N (NOTE: by dropping NaN
                            % F-scores from sorting to select threshold, we'll
                            % also implicitly exclude those voxels - so don't
                            % need to chop them off first, which creates
                            % indexing complications)
                        else
                            
                            
                            allF = nan(size(trndata,2),1);
                            allp = nan(size(trndata,2),1);
                            thisG = data.c_map(trn_idx,2);
                            for voxidx = 1:size(trndata,2)
                                
                                thisX = trndata(:,voxidx);
                                
                                
                                [p,tab,stats] = anova1(thisX,thisG,'off');
                                allF(voxidx) = tab{2,5};
                                allp(voxidx) = p;
                                clear thisX p tab stats;
                            end
                            
                            f_sorted = sort(allF(~isnan(allF)),'descend');
                            if which_vox <= length(allF) % handle case of small ROI
                                f_thresh = f_sorted(which_vox);
                            else
                                f_thresh = f_sorted(end);
                            end
                            
                            these_vox = allF>=f_thresh;
                            if sum(these_vox)>which_vox
                                allidx = find(these_vox);
                                these_vox(allidx((which_vox+1):end)) = 0;
                            end
                            
                        end
                        
                        trndata = trndata(:,these_vox);
                        tstdata = tstdata(:,these_vox);
                        
                        clear mystd;
                        
                        
                        % n_trials x n_vox (shorter var names...)
                        trn = trndata;
                        tst = tstdata;
                        
                        X = spDist_makeX1(data.(align_to{aa})(trn_idx),chan_centers);
                        X = X./max(X(:)); % normalize design matrix to 1
                        
                        w = X\trn;
                        %w_all{trn_tpt_idx,tst_tpt_idx,fold_idx} = w;
                        
                        % channel responses
                        chan_resp{trn_tpt_idx,tst_tpt_idx,aa}(tst_idx,:) = (inv(w*w.')*w*tst.').';

                        clear trn tst trn_idx tst_idx trn_runs tst_runs X_trn w trndata tstdata allF f_thresh thisG;

                        
                    end
                    
                    
                    for tt = 1:size(chan_resp{trn_tpt_idx,tst_tpt_idx,aa},1)
                        
                        % we have to build a basis set for each trial, each target
                        
                        % remove the polar angle of the aligned target
                        rot_by = data.(align_to{aa})(tt,1);%atan2d(data.xy_task(tt,2),data.xy_task(tt,1));
                        
                        this_rfTh = chan_centers-rot_by; % rotate basis
                        myb = build_basis_polar_mat(angs,this_rfTh);
                        
                        % myb is length(angs) x n_channels
                        % we want to weight each channel (col) by this trials' channel
                        % activation
                        % (result shoudl be 1 x length(angs))
                        
                        myb_orig = build_basis_polar_mat(angs,chan_centers);
                        
                        recons{trn_tpt_idx,tst_tpt_idx,aa}(tt,:) = (myb * chan_resp{trn_tpt_idx,tst_tpt_idx,aa}(tt,:).').';
                        
                        recons_raw{trn_tpt_idx,tst_tpt_idx,aa}(tt,:) = (myb_orig * chan_resp{trn_tpt_idx,tst_tpt_idx,aa}(tt,:).').';
                        
                    end

                    
                    
                    
                end
            end
        end
  
        
        
        
      
       
       
       % for reference: 
       % (saved from spDist_scoreEyeData.m; via eyeData)
       
       % 1: distractor condition (1 = no, 2 = distractor)
       % 2,3: target position X, Y (dva)
       % 4,5: distractor position X,Y (or NaN; dva)
       % 6: distractor bin (-3:3; NaN); Cartesian, so + is CCW
       % 7: distractor direction (1 = CCW, 2 = CW, NaN = no dst)
       % 8: correct? (0 or 1; NaN)
       % 9: RT
       
       % NOTE: we don't do this for GAT (yet)
       % average distractor representation on distractor trials
%        dist_mu = squeeze(mean( recons{2}(data.c_all(:,1)==2,:,:), 1 )); % n_angs x n_tpts
%        
%        % variable for saving distractor-corrected single-item WM
%        % representations
%        recons_nodist = nan(size(recons{1})); 
%        
%        % angular difference - distractor relative to target (+ is CCW, -
%        % CW)
%        ang_diff = atan2d( sind(data.dist_ang_all-data.targ_ang_all), cosd(data.dist_ang_all-data.targ_ang_all) );
%        
%        zero_idx = find(angs==0); % 45
%        
%        % remove that average from each trial, after aligning based on
%        % relative distractor location per trial
%        for tt = 1:size(recons{1},1)
%            
%            % if distractor trial, remove average of all distractors
%            if data.c_all(tt,1)==2
%                
%                % distance between distractor and WM position
%                % ang_diff(tt)
%                
%                % find nearest ang bin - how many units to move?
%                [~,ang_diff_idx] = min(abs(angs-ang_diff(tt)));
%                
%                
%                % circularly shift thismu by (ang_diff_idx-zeroidx)
%                % (if positive, distractor is +ang compared to targ, which
%                % is right, so shift RIGHT) (TCS, updated 5/9/2019)
%                recons_nodist(tt,:,:) = recons{1}(tt,:,:) - shiftdim(circshift(dist_mu,(ang_diff_idx-zero_idx),1),-1);
%                
%            end
%        end
%        
       
       % things we want to save
        
        c_all = data.c_all;
        r_all = data.r_all;
        p_all = data.p_all;
        a_all = [data.targ_ang_all data.dist_ang_all];
        sess_all  = data.sess;
        
        
        % save with VE marker when which_vox < 1, otherwise, number of
        % vox
        if which_vox < 1
            fn2s = sprintf('%s/%s_reconstructions/%s_%s_%s_%s_%ichan_VE%03.f_GATdist.mat',root,tst_dir,subj{ss},horzcat(sess{ss}{:}),ROIs{vv},func_suffix,n_chan,100*which_vox);
        else
            fn2s = sprintf('%s/%s_reconstructions/%s_%s_%s_%s_%ichan_%ivox_GATdist.mat'  ,root,tst_dir,subj{ss},horzcat(sess{ss}{:}),ROIs{vv},func_suffix,n_chan,    which_vox);
        end
        fprintf('saving to %s...\n',fn2s);
        
        save(fn2s,'c_all','r_all','p_all','n_chan','delay_tpts','angs','recons','recons_raw','chan_resp','which_vox','sess_all','these_vox','a_all');
        
        clear data c_all r_all p_all chan_resp recons a_all;
        
        
        
    end
    
    
end
return