function spDist_plotEyeData(subj,sess,WHICH_EXCL)


root = spDist_loadRoot;

if nargin < 1
    subj = {'KD','CC','AY','MR','XL'};
end
%subj = {'CC'};
if nargin < 2
    % TODO: if no defined sessions, use all....
    sess_template = {'spDist1','spDist2'};
    sess = cell(length(subj),1); for ss = 1:length(subj); sess{ss} = sess_template; end 
    clear sess_template
    %sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}};
end

if nargin < 3
    WHICH_EXCL = [13 20 21 22]; % don't exclude trials w/ calibration failures for now...
    %WHICH_EXCL = [20 21 22]; % don't exclude trials w/ calibration failures for now...
end


% for now, let's use cat_struct to load/concatenate all data...
all_subj = nan(1000*length(subj),1);
%u_subj = unique(cellfun(@(s) s(1:2),subj,'uniformoutput',0));

TARG_ECC = 12;

niter = 1000;

all_data = [];

startidx = 1;

for ss = 1:length(subj)
    for sessidx = 1:length(sess{ss})
%     fn = sprintf('%s/data/%s_wmPri_scored.mat',root,subj{ss});
%     fprintf('Loading trial information from %s\n',fn);
%     this_data = load(fn);
    
    fn = sprintf('%s/spDist_behav/%s_%s_scored.mat',root,subj{ss},sess{ss}{sessidx});
    fprintf('Loading scored eye data from %s\n',fn);
    this_scored = load(fn);
    
    this_data.s_all = this_scored.ii_sess;
    this_data.sess_all = sessidx;
    
    this_subj = ss;%find(strcmpi(u_subj,subj{ss}(1:2)));
    
    all_data = cat_struct(all_data,this_data);
    all_subj(startidx:(startidx-1+size(this_scored.ii_sess.trialinfo,1))) = this_subj;
    
    startidx = startidx+size(this_scored.ii_sess.trialinfo,1);
    
    clear this_subj this_data;
    end
end

% let's try this pattern for now
all_subj = all_subj(1:(startidx-1));
all_data.subj_all = all_subj;

% determine which trials to include
% first, narrow based on saccade preprocessing/scoring exclusions
% (wmChoose_extractSaccadeData1.m)
all_data.use_trial = ~cellfun( @any, cellfun( @(a) ismember(a, WHICH_EXCL), all_data.s_all.excl_trial, 'UniformOutput',false));

% drop trials with very short (< 100 ms) or very long RT (> 1 s)
all_data.use_trial(all_data.s_all.i_sacc_rt<0.1 | all_data.s_all.i_sacc_rt>1.0) = 0;
%all_data.use_trial(all_data.s_all.i_sacc_err>5) = 0; % kill 'bad' trials (errors)


%% compare error (std dev) during distractor/no-distractor trials
%  for distractor trials, compare each bin

distractor_bins = unique(all_data.s_all.trialinfo(all_data.s_all.trialinfo(:,1)~=1,6));
distractor_spacing = 360/length(distractor_bins);


% within each:
% - row: bin
% - col: initial, final
% - page1:rad/tang (x/y)
% - page2:subj

all_err = cell(2,1);
all_mu = cell(2,1);

cond_str = {'No distractor','Distractor'};
cond_colors = lines(2);

params_of_interest = {'i_sacc','f_sacc'};
param_str = {'Primary saccade','Final saccade'};

% first, error for no-distractor trials
tmp_err = nan(1,length(params_of_interest),2,length(subj)); % initial, final
tmp_mu  = nan(1,length(params_of_interest),2,length(subj)); % initial, final
for ss = 1:length(subj)
    thisidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==1 & all_data.use_trial==1;
    for pp = 1:length(params_of_interest)
        % distractor bin x param x [radial; tangential] x subj
        tmp_err(1,pp,:,ss) = std( all_data.s_all.(params_of_interest{pp})(thisidx,:), [], 1 );
        tmp_mu(1,pp,:,ss)  = mean( all_data.s_all.(params_of_interest{pp})(thisidx,:), 1 );
    end
end

all_err{1} = tmp_err; 
all_mu{2} = tmp_mu;
clear tmp_err tmp_mu;


tmp_err = nan(length(distractor_bins),length(params_of_interest),2,length(subj));
tmp_mu = nan(length(distractor_bins),length(params_of_interest),2,length(subj));

for ss = 1:length(subj)
    for bb = 1:length(distractor_bins)
        thisidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==2 & all_data.use_trial==1 & all_data.s_all.trialinfo(:,6)==distractor_bins(bb);
        for pp = 1:length(params_of_interest)
            % distractor bin x param x [radial; tangential] x subj
            tmp_err(bb,pp,:,ss) = std( all_data.s_all.(params_of_interest{pp})(thisidx,:), [], 1 );
            tmp_mu(bb,pp,:,ss) = mean( all_data.s_all.(params_of_interest{pp})(thisidx,:),  1 );

        end
    end
end

all_err{2} = tmp_err; 
all_mu{2} = tmp_mu;
clear tmp_err tmp_mu;

%% first figure: just directly compare error for no- and with-distractor
% trials (subplot for each param); averaged over radial/tang...
figure;
for pp = 1:length(params_of_interest)
    
    subplot(1,length(params_of_interest),pp); hold on;
    
    % distractor cond x subj
    thise = nan(length(all_err),length(subj));
    for ii = 1:length(all_err)
        thise(ii,:) = mean(mean(all_err{ii}(:,pp,:,:),3),1); % mean over radial/tangential; distractor bin
    end
    
    %plot(1:size(thise,1),thise,'-','Color',[0.3 0.3 0.3]);
    %plot(1:size(thise,1),thise,'o-','Color',[0.5 0.5 0.5],'MarkerFaceColor','w','MarkerSize',5,'LineWidth',1);
    for ii = 1:length(all_err)
        %plot(ii+[-0.2 0.2],[1 1]*mean(thise(ii,:)),'-','LineWidth',2,'Color',cond_colors(ii,:));
        tmpe = std(thise(ii,:))/sqrt(length(subj));
        plot(ii*[1 1],mean(thise(ii,:))+tmpe*[-1 1],'-','LineWidth',1.5,'Color',cond_colors(ii,:));
        plot(ii,mean(thise(ii,:)),'o','LineWidth',1.5,'Color',cond_colors(ii,:),'MarkerSize',7,'MarkerFaceColor','w');
    end
    
    xlim([0 length(all_err)+1]);
    ylim([0 2.5]);
    set(gca,'XTick',1:length(all_err),'XTickLabel',cond_str','XTickLabelRotation',-45,'TickDir','out');
    xlabel('Condition');
    if pp == 1
        ylabel('Precision (avg std dev, \circ)');
    end
    title(param_str{pp});
    ylim([0 2]);
end



% second figure: error as a function of distractor bin
figure;


for pp = 1:length(params_of_interest)
    
    subplot(1,length(params_of_interest),pp); hold on;
    
    % nbins x nsubj
    thise = squeeze(mean(all_err{2}(:,pp,:,:),3)); % nbins x subj
    
    plot(distractor_bins*distractor_spacing,thise,'-','Color',[0.3 0.3 0.3]);
    plot(distractor_bins*distractor_spacing,mean(thise,2),'o','MarkerSize',5,'Color',cond_colors(2,:),'MarkerFaceColor',cond_colors(2,:));
    
    title(param_str{pp});
    xlabel('Distractor offset (\circ polar angle)');
    if pp == 1
        ylabel('Precision (avg std dev, \circ)');
    end
    
    set(gca,'XTick',distractor_bins*distractor_spacing,'TickDir','out','XTickLabelRotation',-45)
end
match_ylim(get(gcf,'Children'));

% third figure: bias as a function of distractor bin
% NOTE: + bins of distractor correspond to clockwise, which is NEGATIVE in
% y - so we plot -1*that....

figure;


for pp = 1:length(params_of_interest)
    
    subplot(1,length(params_of_interest),pp); hold on;
    
    % nbins x nsubj
    thism = -1*squeeze(all_mu{2}(:,pp,2,:)); % nbins x subj <--- extract 'y'
    
    plot(distractor_bins*distractor_spacing,thism,'-','Color',[0.3 0.3 0.3]);
    plot(distractor_bins*distractor_spacing,mean(thism,2),'o','MarkerSize',5,'Color',cond_colors(2,:),'MarkerFaceColor',cond_colors(2,:));
    
    title(param_str{pp});
    xlabel('Distractor offset (\circ polar angle)');
    if pp == 1
        ylabel('Bias (\circ)');
    end
    ylim([-3 3]);
    set(gca,'XTick',distractor_bins*distractor_spacing,'TickDir','out','XTickLabelRotation',-45)
end
match_ylim(get(gcf,'Children'));