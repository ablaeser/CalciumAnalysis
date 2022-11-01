minIso = 10;
for x = xPrePost
    figure('WindowState','maximized')
    for runs = find([periBout{x}.Nbout] > 0) % 1:expt{x}.Nruns

        goodBouts = find(periBout{x}(runs).iso(:,1) > minIso)';

        %for bout = 1:periBout{x}(runs).Nbout
        %end
        subplot(2,1,1);
        plot( periBout{x}(runs).Tstart(goodBouts)-csdBout{x}(expt{x}.csd).Tstart,  periBout{x}(runs).stat.speed.run(goodBouts,:,1), '.'); hold on;
        ylabel('Speed effect');

        subplot(2,1,2);
        plot( periBout{x}(runs).Tstart(goodBouts)-csdBout{x}(expt{x}.csd).Tstart,  mean(periBout{x}(runs).stat.fluor.effect(goodBouts,boutResponse(x).pre.exc,1), 2, 'omitnan'), '.'); hold on;
        ylabel('dF/F');

    end
    pause;
    clf;
end
close all

%% Histograms of bout properties
bout_dur = cell(1,Nexpt);
for x = xPresent
    bout_dur{x} = [periBout{x}.dur]';
end

%cell2padmat(bout_dur)
ecdfCol(cell2padmat(bout_dur))

%% Recalculate bouts with new criteria

% Get peri-locomotion, stillness, and peri-CSD bout data
onBout = cell(1,Nexpt); onParam = cell(1,Nexpt);
for x = xPrePost
    for runs = flip(1:expt{x}.Nruns)
        [onBout{x}(runs), onParam{x}(runs), loco{x}(runs)] = PeriLoco3D(expt{x}, Tscan{x}{runs}, loco{x}(runs), deform{x}(runs), fluor{x}(runs).F.ROI, defVars, 'base',10, 'run',2, 'iso',[10,0], 'min_vel_on',3); %periParam(x)
        %[onBout{x}(runs), onParam{x}(runs), loco{x}(runs)] = PeriLoco3D(expt{x}, Tscan{x}{runs}, loco{x}(runs), deform{x}(runs), fluor{x}(runs).F.axon, defVars, 'base',10, 'run',2, 'iso',[10,0], 'min_vel_on',3);
        %InspectPeriDeform3D( expt{x}, periBout{x}(run), {'fluor','scaleMag'} ); %defVars
    end
end

%% Summarize peri-onset data by phase
on = repmat( struct('T',[], 'Nscan',NaN, 'ind_base',[], 'ind_loco',[], 'Nroi',NaN, 'pre',[], 'acute',[], 'post',[]), 1, Nexpt );
classPct = nan(Nexpt, 3, 3); % FOV x class x phase
Tbase = -2; Tloco = 0;
minDiff = 0.05;
minBouts = 3;

%figure('WindowState','maximized')
for x = xPrePost
    %clf;
    % Break the onBout structure up by phase
    on_pre = [onBout{x}(expt{x}.preRuns).on];
    on_acute = [onBout{x}(expt{x}.acuteRuns).on];
    on_post = [onBout{x}(expt{x}.postRuns).on];
    % Get peri-onset time vector and find baseline and locomotion time ranges
    %on(x).T = [on_pre.T, on_acute.T, on_post.T];
    on(x).T = on_pre(1).T; %on(x).T(:,1);
    on(x).Nscan = numel(on(x).T);
    on(x).ind_base = find( on(x).T <= Tbase );
    on(x).ind_loco = find( on(x).T >= Tloco );

    % PRE-CSD PHASE
    on(x).pre.velocity = [on_pre.velocity];
    on(x).pre.bouts = find(all(~isnan(on(x).pre.velocity), 1));
    on(x).pre.fluor.raw = cat(2, on_pre.fluor);
    on(x).Nroi = size(on(x).pre.fluor.raw,3);
    on(x).pre.Nbout = numel(on(x).pre.bouts);
    if on(x).pre.Nbout > 0
        on(x).pre.velocity_mean = mean(on(x).pre.velocity(:,on(x).pre.bouts), 2);
        % Gather, normalize and average fluor data
        on(x).pre.fluor.base = median(on(x).pre.fluor.raw(on(x).ind_base,:,:), 1, 'omitnan'); %
        on(x).pre.fluor.loco = max(on(x).pre.fluor.raw(on(x).ind_loco,:,:), [], 1);
        on(x).pre.fluor.norm = (on(x).pre.fluor.raw-on(x).pre.fluor.base)./on(x).pre.fluor.base;
        on(x).pre.fluor.norm_diff = (on(x).pre.fluor.loco - on(x).pre.fluor.base)./on(x).pre.fluor.base;
        % Average over ROIs for each bout
        on(x).pre.fluor.norm_bout = squeeze(mean(on(x).pre.fluor.norm(:,on(x).pre.bouts,:), 3, 'omitnan')); % average of noramlized fluor over all ROIs, for each bout
        on(x).pre.fluor.norm_bout_mean = mean(on(x).pre.fluor.norm_bout, 2); % average over all bouts
        % Average over bouts for each ROI
        on(x).pre.fluor.norm_ROI = squeeze(mean(on(x).pre.fluor.norm(:,on(x).pre.bouts,:), 2, 'omitnan'));
        on(x).pre.fluor.norm_ROI_mean = mean(on(x).pre.fluor.norm_ROI, 2);

        % Gather, normalize and average deformation data
        for v = 1:NdefVars
            on(x).pre.(defVars{v}).raw = cat(2, on_pre.(defVars{v}));
            % Gather, normalize and average fluor data
            on(x).pre.(defVars{v}).base = median(on(x).pre.(defVars{v}).raw(on(x).ind_base,:,:), 1); %, 'omitnan'
            on(x).pre.(defVars{v}).loco = max(on(x).pre.(defVars{v}).raw(on(x).ind_loco,:,:), [], 1);
            on(x).pre.(defVars{v}).norm = (on(x).pre.(defVars{v}).raw-on(x).pre.(defVars{v}).base)./on(x).pre.(defVars{v}).base;
            on(x).pre.(defVars{v}).norm_diff = (on(x).pre.(defVars{v}).loco - on(x).pre.(defVars{v}).base)./on(x).pre.(defVars{v}).base;
            % Average over ROIs for each bout
            on(x).pre.(defVars{v}).norm_bout = squeeze(mean(on(x).pre.(defVars{v}).norm(:,on(x).pre.bouts,:), 3, 'omitnan')); % average of noramlized (defVars{v}) over all ROIs, for each bout
            on(x).pre.(defVars{v}).norm_bout_mean = mean(on(x).pre.(defVars{v}).norm_bout, 2); % average over all bouts
            % Average over bouts for each ROI
            on(x).pre.(defVars{v}).norm_plane = squeeze(mean(on(x).pre.(defVars{v}).norm(:,on(x).pre.bouts,:), 2, 'omitnan'));
            on(x).pre.(defVars{v}).norm_plane_mean = mean(on(x).pre.(defVars{v}).norm_plane, 2);
        end

    else
        on(x).pre.fluor.base = nan(1,0,on(x).Nroi); %
        on(x).pre.fluor.loco = nan(1,0,on(x).Nroi); %
        on(x).pre.fluor.norm = nan(on(x).Nscan,0,on(x).Nroi);
        on(x).pre.fluor.norm_diff = nan(1,0,on(x).Nroi);
        on(x).pre.fluor.norm_bout = nan(on(x).Nscan, 0);
        on(x).pre.fluor.norm_bout_mean = nan(on(x).Nscan, 1);
        on(x).pre.fluor.norm_ROI = nan(on(x).Nscan, on(x).Nroi);
        on(x).pre.fluor.norm_ROI_mean = nan(on(x).Nscan, 1);
    end
    % Classify ROIs based on response to pre-CSD locomotion onset
    if on(x).pre.Nbout >= minBouts %&& on(x).post.Nbout > minBouts
        pre_norm_diff = squeeze(on(x).pre.fluor.norm_diff);
        [~,pre_norm_diff_p] = ttest(pre_norm_diff);
        on(x).pre.fluor.norm_diff_mean = mean(pre_norm_diff, 1, 'omitnan');
        on(x).pre.fluor.roi_act = find(on(x).pre.fluor.norm_diff_mean >= minDiff & pre_norm_diff_p < 0.05);
        on(x).pre.fluor.roi_supp = find(on(x).pre.fluor.norm_diff_mean <= -minDiff & pre_norm_diff_p < 0.05);
        on(x).pre.fluor.roi_neut = setdiff(1:on(x).Nroi, [on(x).pre.fluor.roi_act, on(x).pre.fluor.roi_supp]); % find(abs(on(x).pre.fluor.norm_diff_mean) < minDiff);

        % Sort ROIs by norm_diff_mean, within class
        [~,sortInd] = sort(on(x).pre.fluor.norm_diff_mean(on(x).pre.fluor.roi_act), 'descend');
        on(x).pre.fluor.roi_act = on(x).pre.fluor.roi_act(sortInd);
        [~,sortInd] = sort(on(x).pre.fluor.norm_diff_mean(on(x).pre.fluor.roi_neut), 'descend');
        on(x).pre.fluor.roi_neut = on(x).pre.fluor.roi_neut(sortInd);
        [~,sortInd] = sort(on(x).pre.fluor.norm_diff_mean(on(x).pre.fluor.roi_supp), 'descend');
        on(x).pre.fluor.roi_supp = on(x).pre.fluor.roi_supp(sortInd);
        on(x).pre.fluor.Nroi_act = numel(on(x).pre.fluor.roi_act);
        on(x).pre.fluor.Nroi_neut = numel(on(x).pre.fluor.roi_neut);
        on(x).pre.fluor.Nroi_supp = numel(on(x).pre.fluor.roi_supp);
        %{
        for roi = on(x).pre.fluor.roi_act % on(x).pre.fluor.roi_supp %  find(abs(on(x).pre.fluor.norm_diff_mean) > minDiff) % 1:on(x).Nroi
            cla;
            JitterPlot(cell2padmat({pre_norm_diff(:,roi), acute_norm_diff(:,roi), post_norm_diff(:,roi)}), 'new',false); hold on;
            line([0,4], [0,0], 'color','k');
            set(gca,'Xtick',1:3, 'XtickLabel',{'Pre','Acute','Post'})
            pause;
        end
        %}
    else
        on(x).pre.fluor.norm_diff_mean = nan(1,on(x).Nroi);
        on(x).pre.fluor.roi_act = [];
        on(x).pre.fluor.roi_neut = [];
        on(x).pre.fluor.roi_supp = [];
        on(x).pre.fluor.Nroi_act = NaN;
        on(x).pre.fluor.Nroi_neut = NaN;
        on(x).pre.fluor.Nroi_supp = NaN;
    end

    % Acute peri-CSD phase
    %{
    on(x).acute.velocity = [on_acute.velocity];
    on(x).acute.bouts = find(all(~isnan(on(x).acute.velocity), 1));
    on(x).acute.fluor.raw = cat(2, on_acute.fluor);
    on(x).acute.Nbout = numel(on(x).acute.bouts);
    if on(x).acute.Nbout > 0
        on(x).acute.velocity_mean = mean(on(x).acute.velocity(:,on(x).acute.bouts), 2);
        % Gather, normalize and average fluor data
        on(x).acute.fluor.base = median(on(x).acute.fluor.raw(on(x).ind_base,:,:), 1); %, 'omitnan'
        on(x).acute.fluor.loco = max(on(x).acute.fluor.raw(on(x).ind_loco,:,:), [], 1);
        on(x).acute.fluor.norm = (on(x).acute.fluor.raw-on(x).acute.fluor.base)./on(x).acute.fluor.base;
        on(x).acute.fluor.norm_diff = (on(x).acute.fluor.loco - on(x).acute.fluor.base)./on(x).acute.fluor.base;
        % Average over ROIs for each bout
        on(x).acute.fluor.norm_bout = squeeze(mean(on(x).acute.fluor.norm(:,on(x).acute.bouts,:), 3, 'omitnan')); % average of noramlized fluor over all ROIs, for each bout
        on(x).acute.fluor.norm_bout_mean = mean(on(x).acute.fluor.norm_bout, 2); % average over all bouts
        % Average over bouts for each ROI
        on(x).acute.fluor.norm_ROI = squeeze(mean(on(x).acute.fluor.norm(:,on(x).acute.bouts,:), 2, 'omitnan'));
        on(x).acute.fluor.norm_ROI_mean = mean(on(x).acute.fluor.norm_ROI, 2);

        % Gather, normalize and average deformation data
        for v = 1:NdefVars
            on(x).acute.(defVars{v}).raw = cat(2, on_acute.(defVars{v}));
            % Gather, normalize and average fluor data
            on(x).acute.(defVars{v}).base = median(on(x).acute.(defVars{v}).raw(on(x).ind_base,:,:), 1); %, 'omitnan'
            on(x).acute.(defVars{v}).loco = max(on(x).acute.(defVars{v}).raw(on(x).ind_loco,:,:), [], 1);
            on(x).acute.(defVars{v}).norm = (on(x).acute.(defVars{v}).raw-on(x).acute.(defVars{v}).base)./on(x).acute.(defVars{v}).base;
            on(x).acute.(defVars{v}).norm_diff = (on(x).acute.(defVars{v}).loco - on(x).acute.(defVars{v}).base)./on(x).acute.(defVars{v}).base;
            % Average over ROIs for each bout
            on(x).acute.(defVars{v}).norm_bout = squeeze(mean(on(x).acute.(defVars{v}).norm(:,on(x).acute.bouts,:), 3, 'omitnan')); % average of noramlized (defVars{v}) over all ROIs, for each bout
            on(x).acute.(defVars{v}).norm_bout_mean = mean(on(x).acute.(defVars{v}).norm_bout, 2); % average over all bouts
            % Average over bouts for each ROI
            on(x).acute.(defVars{v}).norm_plane = squeeze(mean(on(x).acute.(defVars{v}).norm(:,on(x).acute.bouts,:), 2, 'omitnan'));
            on(x).acute.(defVars{v}).norm_plane_mean = mean(on(x).acute.(defVars{v}).norm_plane, 2);
        end
    else
        on(x).acute.fluor.base = nan(1,0,on(x).Nroi); %
        on(x).acute.fluor.loco = nan(1,0,on(x).Nroi); %
        on(x).acute.fluor.norm = nan(on(x).Nscan,0,on(x).Nroi);
        on(x).acute.fluor.norm_diff = nan(1,0,on(x).Nroi);
        on(x).acute.fluor.norm_bout = nan(on(x).Nscan, 0);
        on(x).acute.fluor.norm_bout_mean = nan(on(x).Nscan, 1);
        on(x).acute.fluor.norm_ROI = nan(on(x).Nscan, on(x).Nroi);
        on(x).acute.fluor.norm_ROI_mean = nan(on(x).Nscan, 1);
    end
    % Classify ROIs based on response to acute-phase locomotion onset
    if on(x).acute.Nbout >= minBouts %&& on(x).post.Nbout > minBouts
        acute_norm_diff = squeeze(on(x).acute.fluor.norm_diff);
        [~,acute_norm_diff_p] = ttest(acute_norm_diff);
        on(x).acute.fluor.norm_diff_mean = mean(acute_norm_diff, 1, 'omitnan');
        on(x).acute.fluor.roi_act = find(on(x).acute.fluor.norm_diff_mean > minDiff & acute_norm_diff_p < 0.05);
        on(x).acute.fluor.roi_supp = find(on(x).acute.fluor.norm_diff_mean < -minDiff & acute_norm_diff_p < 0.05);
        on(x).acute.fluor.roi_neut = setdiff(1:on(x).Nroi, [on(x).acute.fluor.roi_act, on(x).acute.fluor.roi_supp]); % find(abs(on(x).acute.fluor.norm_diff_mean) < minDiff);

        % Sort ROIs by norm_diff_mean, within class
        [~,sortInd] = sort(on(x).acute.fluor.norm_diff_mean(on(x).acute.fluor.roi_act), 'descend');
        on(x).acute.fluor.roi_act = on(x).acute.fluor.roi_act(sortInd);
        [~,sortInd] = sort(on(x).acute.fluor.norm_diff_mean(on(x).acute.fluor.roi_neut), 'descend');
        on(x).acute.fluor.roi_neut = on(x).acute.fluor.roi_neut(sortInd);
        [~,sortInd] = sort(on(x).acute.fluor.norm_diff_mean(on(x).acute.fluor.roi_supp), 'descend');
        on(x).acute.fluor.roi_supp = on(x).acute.fluor.roi_supp(sortInd);
        on(x).acute.fluor.Nroi_act = numel(on(x).acute.fluor.roi_act);
        on(x).acute.fluor.Nroi_neut = numel(on(x).acute.fluor.roi_neut);
        on(x).acute.fluor.Nroi_supp = numel(on(x).acute.fluor.roi_supp);
        %{
        for roi = on(x).acute.fluor.roi_supp %  find(abs(on(x).acute.fluor.norm_diff_mean) > minDiff) % 1:on(x).Nroi
            cla;
            JitterPlot(cell2padmat({acute_norm_diff(:,roi), acute_norm_diff(:,roi), post_norm_diff(:,roi)}), 'new',false); hold on;
            line([0,4], [0,0], 'color','k');
            set(gca,'Xtick',1:3, 'XtickLabel',{'Pre','Acute','Post'})
            pause;
        end
        %}
    else
        on(x).acute.fluor.norm_diff_mean = nan(1,on(x).Nroi);
        on(x).acute.fluor.roi_act = [];
        on(x).acute.fluor.roi_neut = [];
        on(x).acute.fluor.roi_supp = [];
        on(x).acute.fluor.Nroi_act = NaN;
        on(x).acute.fluor.Nroi_neut = NaN;
        on(x).acute.fluor.Nroi_supp = NaN;
    end
    %}


    % POST CSD PHASE
    on(x).post.velocity = [on_post.velocity];
    on(x).post.bouts = find(all(~isnan(on(x).post.velocity), 1));
    on(x).post.fluor.raw = cat(2, on_post.fluor);
    on(x).post.Nbout = numel(on(x).post.bouts);
    if on(x).post.Nbout > 0
        on(x).post.velocity_mean = mean(on(x).post.velocity(:,on(x).post.bouts), 2);
        % Gather, normalize and average fluor data
        on(x).post.fluor.base = median(on(x).post.fluor.raw(on(x).ind_base,:,:), 1); %, 'omitnan'
        on(x).post.fluor.loco = max(on(x).post.fluor.raw(on(x).ind_loco,:,:), [], 1);
        on(x).post.fluor.norm = (on(x).post.fluor.raw-on(x).post.fluor.base)./on(x).post.fluor.base;
        on(x).post.fluor.norm_diff = (on(x).post.fluor.loco - on(x).post.fluor.base)./on(x).post.fluor.base;
        % Average over ROIs for each bout
        on(x).post.fluor.norm_bout = squeeze(mean(on(x).post.fluor.norm(:,on(x).post.bouts,:), 3, 'omitnan')); % average of noramlized fluor over all ROIs, for each bout
        on(x).post.fluor.norm_bout_mean = mean(on(x).post.fluor.norm_bout, 2); % average over all bouts
        % Average over bouts for each ROI
        on(x).post.fluor.norm_ROI = squeeze(mean(on(x).post.fluor.norm(:,on(x).post.bouts,:), 2, 'omitnan'));
        on(x).post.fluor.norm_ROI_mean = mean(on(x).post.fluor.norm_ROI, 2);

        % Gather, normalize and average deformation data
        for v = 1:NdefVars
            on(x).post.(defVars{v}).raw = cat(2, on_post.(defVars{v}));
            % Gather, normalize and average fluor data
            on(x).post.(defVars{v}).base = median(on(x).post.(defVars{v}).raw(on(x).ind_base,:,:), 1); %, 'omitnan'
            on(x).post.(defVars{v}).loco = max(on(x).post.(defVars{v}).raw(on(x).ind_loco,:,:), [], 1);
            on(x).post.(defVars{v}).norm = (on(x).post.(defVars{v}).raw-on(x).post.(defVars{v}).base)./on(x).post.(defVars{v}).base;
            on(x).post.(defVars{v}).norm_diff = (on(x).post.(defVars{v}).loco - on(x).post.(defVars{v}).base)./on(x).post.(defVars{v}).base;
            % Average over ROIs for each bout
            on(x).post.(defVars{v}).norm_bout = squeeze(mean(on(x).post.(defVars{v}).norm(:,on(x).post.bouts,:), 3, 'omitnan')); % average of noramlized (defVars{v}) over all ROIs, for each bout
            on(x).post.(defVars{v}).norm_bout_mean = mean(on(x).post.(defVars{v}).norm_bout, 2); % average over all bouts
            % Average over bouts for each ROI
            on(x).post.(defVars{v}).norm_plane = squeeze(mean(on(x).post.(defVars{v}).norm(:,on(x).post.bouts,:), 2, 'omitnan'));
            on(x).post.(defVars{v}).norm_plane_mean = mean(on(x).post.(defVars{v}).norm_plane, 2);
        end
    else
        on(x).post.fluor.base = nan(1,0,on(x).Nroi); %
        on(x).post.fluor.loco = nan(1,0,on(x).Nroi); %
        on(x).post.fluor.norm = nan(on(x).Nscan,0,on(x).Nroi);
        on(x).post.fluor.norm_diff = nan(1,0,on(x).Nroi);
        on(x).post.fluor.norm_diff_mean = nan(1, on(x).Nroi);
        on(x).post.fluor.norm_bout = nan(on(x).Nscan, 0);
        on(x).post.fluor.norm_bout_mean = nan(on(x).Nscan, 1);
        on(x).post.fluor.norm_ROI = nan(on(x).Nscan, on(x).Nroi);
        on(x).post.fluor.norm_ROI_mean = nan(on(x).Nscan, 1);
    end
    % Classify ROIs based on response to postCSD phase locomotion onset
    if on(x).post.Nbout >= minBouts %&& on(x).post.Nbout > minBouts
        post_norm_diff = squeeze(on(x).post.fluor.norm_diff);
        [~,post_norm_diff_p] = ttest(post_norm_diff);
        on(x).post.fluor.norm_diff_mean = mean(post_norm_diff, 1, 'omitnan');
        on(x).post.fluor.roi_act = find(on(x).post.fluor.norm_diff_mean > minDiff & post_norm_diff_p < 0.05);
        on(x).post.fluor.roi_supp = find(on(x).post.fluor.norm_diff_mean < -minDiff & post_norm_diff_p < 0.05);
        on(x).post.fluor.roi_neut = setdiff(1:on(x).Nroi, [on(x).post.fluor.roi_act, on(x).post.fluor.roi_supp]); %find(abs(on(x).post.fluor.norm_diff_mean) < minDiff);

        % Sort ROIs by norm_diff_mean, within class
        [~,sortInd] = sort(on(x).post.fluor.norm_diff_mean(on(x).post.fluor.roi_act), 'descend');
        on(x).post.fluor.roi_act = on(x).post.fluor.roi_act(sortInd);
        [~,sortInd] = sort(on(x).post.fluor.norm_diff_mean(on(x).post.fluor.roi_neut), 'descend');
        on(x).post.fluor.roi_neut = on(x).post.fluor.roi_neut(sortInd);
        [~,sortInd] = sort(on(x).post.fluor.norm_diff_mean(on(x).post.fluor.roi_supp), 'descend');
        on(x).post.fluor.roi_supp = on(x).post.fluor.roi_supp(sortInd);
        on(x).post.fluor.Nroi_act = numel(on(x).post.fluor.roi_act);
        on(x).post.fluor.Nroi_neut = numel(on(x).post.fluor.roi_neut);
        on(x).post.fluor.Nroi_supp = numel(on(x).post.fluor.roi_supp);

        %{
        for roi = on(x).post.fluor.roi_supp %  find(abs(on(x).post.fluor.norm_diff_mean) > minDiff) % 1:on(x).Nroi
            cla;
            JitterPlot(cell2padmat({post_norm_diff(:,roi), post_norm_diff(:,roi), post_norm_diff(:,roi)}), 'new',false); hold on;
            line([0,4], [0,0], 'color','k');
            set(gca,'Xtick',1:3, 'XtickLabel',{'Pre','Acute','Post'})
            pause;
        end
        %}
    else
        on(x).post.fluor.norm_diff_mean = nan(1,on(x).Nroi);
        on(x).post.fluor.roi_act = [];
        on(x).post.fluor.roi_neut = [];
        on(x).post.fluor.roi_supp = [];
        on(x).post.fluor.Nroi_act = NaN;
        on(x).post.fluor.Nroi_neut = NaN;
        on(x).post.fluor.Nroi_supp = NaN;
    end

    % How did ROI class change pre vs post CSD?
    classPct(x,1,1) = on(x).pre.fluor.Nroi_act/on(x).Nroi;
    classPct(x,2,1) = on(x).pre.fluor.Nroi_neut/on(x).Nroi;
    classPct(x,3,1) = on(x).pre.fluor.Nroi_supp/on(x).Nroi;
    %classPct(x,1,2) = on(x).acute.fluor.Nroi_act/on(x).Nroi;
    %classPct(x,2,2) = on(x).acute.fluor.Nroi_neut/on(x).Nroi;
    %classPct(x,3,2) = on(x).acute.fluor.Nroi_supp/on(x).Nroi;
    classPct(x,1,3) = on(x).post.fluor.Nroi_act/on(x).Nroi;
    classPct(x,2,3) = on(x).post.fluor.Nroi_neut/on(x).Nroi;
    classPct(x,3,3) = on(x).post.fluor.Nroi_supp/on(x).Nroi;

    % Act -> act
    on(x).trans.roi_a2a = intersect(on(x).pre.fluor.roi_act, on(x).post.fluor.roi_act);
    on(x).trans.Nroi_a2a = numel(on(x).trans.roi_a2a);
    on(x).trans.frac_a2a = on(x).trans.Nroi_a2a/on(x).pre.fluor.Nroi_act;
    % Act -> neut
    on(x).trans.roi_a2n = intersect(on(x).pre.fluor.roi_act, on(x).post.fluor.roi_neut);
    on(x).trans.Nroi_a2n = numel(on(x).trans.roi_a2n);
    on(x).trans.frac_a2n = on(x).trans.Nroi_a2n/on(x).pre.fluor.Nroi_act;
    % Act -> supp
    on(x).trans.roi_a2s = intersect(on(x).pre.fluor.roi_act, on(x).post.fluor.roi_supp);
    on(x).trans.Nroi_a2s = numel(on(x).trans.roi_a2s);
    on(x).trans.frac_a2s = on(x).trans.Nroi_a2s/on(x).pre.fluor.Nroi_act;

    % Neut -> act
    on(x).trans.roi_n2a = intersect(on(x).pre.fluor.roi_neut, on(x).post.fluor.roi_act);
    on(x).trans.Nroi_n2a = numel(on(x).trans.roi_n2a);
    on(x).trans.frac_n2a = on(x).trans.Nroi_n2a/on(x).pre.fluor.Nroi_neut;
    % Neut -> neut
    on(x).trans.roi_n2n = intersect(on(x).pre.fluor.roi_neut, on(x).post.fluor.roi_neut);
    on(x).trans.Nroi_n2n = numel(on(x).trans.roi_n2n);
    on(x).trans.frac_n2n = on(x).trans.Nroi_n2n/on(x).pre.fluor.Nroi_neut;
    % Neut -> supp
    on(x).trans.roi_n2s = intersect(on(x).pre.fluor.roi_neut, on(x).post.fluor.roi_supp);
    on(x).trans.Nroi_n2s = numel(on(x).trans.roi_n2s);
    on(x).trans.frac_n2s = on(x).trans.Nroi_n2s/on(x).pre.fluor.Nroi_neut;

    % Supp -> act
    on(x).trans.roi_s2a = intersect(on(x).pre.fluor.roi_supp, on(x).post.fluor.roi_act);
    on(x).trans.Nroi_s2a = numel(on(x).trans.roi_s2a);
    on(x).trans.frac_s2a = on(x).trans.Nroi_s2a/on(x).pre.fluor.Nroi_supp;
    % Supp -> neut
    on(x).trans.roi_s2n = intersect(on(x).pre.fluor.roi_supp, on(x).post.fluor.roi_neut);
    on(x).trans.Nroi_s2n = numel(on(x).trans.roi_s2n);
    on(x).trans.frac_s2n = on(x).trans.Nroi_s2n/on(x).pre.fluor.Nroi_supp;
    % Supp -> supp
    on(x).trans.roi_s2s = intersect(on(x).pre.fluor.roi_supp, on(x).post.fluor.roi_supp);
    on(x).trans.Nroi_s2s = numel(on(x).trans.roi_s2s);
    on(x).trans.frac_s2s = on(x).trans.Nroi_s2s/on(x).pre.fluor.Nroi_supp;

    %pause;
end


%% Summarize peri-onset velocity and fluor response, by phase, for each FOV
diff_lims = 0.3*[-1,1];
figure('WindowState','maximized')
for x = 11 %xPrePost
    clf
    sortROI = [on(x).pre.fluor.roi_act, on(x).pre.fluor.roi_neut, on(x).pre.fluor.roi_supp];
    %if isempty(sortROI), sortROI = [boutResponse(x).pre.exc, boutResponse(x).pre.neut, boutResponse(x).pre.inh]; end
    if on(x).pre.Nbout > 0

        Nans_pre = [on(x).pre.fluor.Nroi_act, on(x).pre.fluor.Nroi_neut, on(x).pre.fluor.Nroi_supp];
        cum_Nans_pre = cumsum(Nans_pre); 

        subplot(2,2,1);
        imagesc(on(x).pre.fluor.norm_ROI(:,sortROI)'); hold on;
        for t = 1:2
            line(0.5+[0,numel(on(x).T)], cum_Nans_pre(t)*[1,1], 'color','k')
        end

        caxis(diff_lims); colormap bluewhitered;
        CB = colorbar('Location','northoutside');
        CB.Label.String = '\DeltaF/F (averaged over all bouts)';
        set(gca,'Xtick',[on(x).ind_base(end), on(x).ind_loco(1)], 'XtickLabel',{'Base','Loco'}, 'Ytick', cum_Nans_pre(1:2), 'YtickLabel',{'Acivated','Neutral'}, 'TickDir','out')
        ylabel('ROI')
        impixelinfo;
        %plot(on(x).T, on(x).pre.fluor.norm_bout_mean, 'color','k', 'linewidth',1.5 );  hold on;
        %plot(on(x).T, on(x).pre.fluor.norm_bout, 'color',[0,0,0,0.3]);
        %plot(on(x).T, on(x).pre.fluor.norm_ROI, 'color',[0,0,0,0.3]); hold on;
        %plot(on(x).T, on(x).pre.fluor.norm_ROI_mean, 'color','k', 'linewidth',1.5 );
        axis square; axis tight;

        subplot(2,2,3);
        plot(on(x).T, on(x).pre.velocity_mean, 'color','k', 'linewidth',1.5);  hold on;
        plot(on(x).T, on(x).pre.velocity(:,on(x).pre.bouts), 'color',[0,0,0,0.3]);
        line(Tbase*[1,1], [-2, 15], 'color','r', 'linestyle','--');
        line(Tloco*[1,1], [-2, 15], 'color','g', 'linestyle','--');
        legend('Boutwise mean','Bouts', 'Location','northwest');
        axis square; axis tight;
        title(sprintf('%i pre-CSD bouts', on(x).pre.Nbout))
        xlabel('Peri-onset time (s)'); ylabel('Velocity (cm/s)')
    end
    %{
    if on(x).acute.Nbout > 0
        subplot(2,3,2);
        imagesc(on(x).acute.fluor.norm_ROI(:,sortROI)');
        caxis(diff_lims); colormap bluewhitered;
        colorbar('Location','northoutside')
        set(gca,'Xtick',[on(x).ind_base(end), on(x).ind_loco(1)], 'XtickLabel',{'Base','Loco'}, 'TickDir','out')
        impixelinfo;
        %plot(on(x).T, on(x).acute.fluor.norm_bout_mean, 'color','k', 'linewidth',1.5 );  hold on;
        %plot(on(x).T, on(x).acute.fluor.norm_bout, 'color',[0,0,0,0.3]);
        %plot(on(x).T, on(x).acute.fluor.norm_ROI, 'color',[0,0,0,0.3]); hold on;
        %plot(on(x).T, on(x).acute.fluor.norm_ROI_mean, 'color','k', 'linewidth',1.5 );
        axis square; axis tight;

        subplot(2,3,5);
        plot(on(x).T, on(x).acute.velocity_mean, 'color','k', 'linewidth',1.5);  hold on;
        plot(on(x).T, on(x).acute.velocity(:,on(x).acute.bouts), 'color',[0,0,0,0.3]);
        line(Tbase*[1,1], [-2, 15], 'color','r', 'linestyle','--');
        line(Tloco*[1,1], [-2, 15], 'color','g', 'linestyle','--');
        legend('Boutwise mean','Bouts', 'Location','northwest');
        axis square; axis tight;
        title(sprintf('%i acute-CSD bouts', on(x).acute.Nbout))
        xlabel('Peri-onset time (s)');
    end
    %}

    if on(x).post.Nbout > 0

        Nans_post = [on(x).post.fluor.Nroi_act, on(x).post.fluor.Nroi_neut, on(x).post.fluor.Nroi_supp];
        cum_Nans_post = cumsum(Nans_post);

        subplot(2,2,2);
        imagesc(on(x).post.fluor.norm_ROI(:,sortROI)');
        colorbar('Location','northoutside')
        caxis(diff_lims); colormap bluewhitered;
        impixelinfo;
        set(gca,'Xtick',[on(x).ind_base(end), on(x).ind_loco(1)], 'XtickLabel',{'Base','Loco'}, 'Ytick', cum_Nans_post(1:2), 'YtickLabel',{'Acivated','Neutral'}, 'TickDir','out')
        %plot(on(x).T, on(x).post.fluor.norm_bout_mean, 'color','k', 'linewidth',1.5 );  hold on;
        %plot(on(x).T, on(x).post.fluor.norm_bout, 'color',[0,0,0,0.3]);
        %plot(on(x).T, on(x).post.fluor.norm_ROI, 'color',[0,0,0,0.3]); hold on;
        %plot(on(x).T, on(x).post.fluor.norm_ROI_mean, 'color','k', 'linewidth',1.5 );
        axis square; axis tight;

        subplot(2,2,4);
        plot(on(x).T, on(x).post.velocity_mean, 'color','k', 'linewidth',1.5);  hold on;
        plot(on(x).T, on(x).post.velocity(:,on(x).post.bouts), 'color',[0,0,0,0.3]);
        line(Tbase*[1,1], [-2, 15], 'color','r', 'linestyle','--');
        line(Tloco*[1,1], [-2, 15], 'color','g', 'linestyle','--');
        legend('Boutwise mean','Bouts', 'Location','northwest');
        axis square; axis tight;
        title(sprintf('%i post-CSD bouts', on(x).post.Nbout))
        xlabel('Peri-onset time (s)');
    end
    pause;
end

%% Plot fractions of ROI in each class pre/post CSD
className = ["Activated","Neutral","Suppressed"];
p_class = nan(1,3);
for cl = 1:3
    tempMat = squeeze(classPct(xPrePost,cl,[1,3]));
    [~,p_class(cl)] = ttest(tempMat(:,1), tempMat(:,2));
    subplot(1,3,cl)
    JitterPlot( squeeze(classPct(xPrePost,cl,[1,3])), 'new',false, 'paired',true  ); axis square;
    ylabel(sprintf('Fraction %s', className(cl)))
    set(gca,'Xtick',[1,2], 'XtickLabel',["Pre","Post"]);
    ylim([0,1]);
    title( sprintf('p = %2.3f', p_class(cl)) );
end

%% Track transitions between classes pre/post
figure;
for x = xPrePost
    trans_mat(1,1) = on(x).trans.frac_a2a; % a2a
    trans_mat(2,1) = on(x).trans.frac_n2a; % n2a
    trans_mat(3,1) = on(x).trans.frac_s2a; % s2a

    trans_mat(1,2) = on(x).trans.frac_a2n; % a2n
    trans_mat(2,2) = on(x).trans.frac_n2n; % n2n
    trans_mat(3,2) = on(x).trans.frac_s2n; % s2n

    trans_mat(1,3) = on(x).trans.frac_a2s; % a2s
    trans_mat(2,3) = on(x).trans.frac_n2s; % n2s
    trans_mat(3,3) = on(x).trans.frac_s2s; % s2s

    imagesc(trans_mat); caxis([0,1]); axis square;
    CB = colorbar; %CB.Position = 'EastOutside';
    set(gca, 'Xtick',1:3, 'XtickLabel',className, 'Ytick',1:3,'YtickLabel',className);
    ylabel('Pre-CSD class'); xlabel('Post-CSD class');
    title(sprintf('x = %i: %s', x, expt{x}.name), 'Interpreter','none')
    impixelinfo;
    pause;
end

%%
Ntrans = nan(9,Nexpt); Nroi = nan(1,Nexpt);
for x = xPrePost
    if on(x).pre.Nbout > minBouts && on(x).post.Nbout > minBouts % ~isempty(on(x).post)
        Ntrans(:,x) = [on(x).trans.Nroi_a2a, on(x).trans.Nroi_a2n, on(x).trans.Nroi_a2s, on(x).trans.Nroi_n2a, on(x).trans.Nroi_n2n, on(x).trans.Nroi_n2s, on(x).trans.Nroi_s2a, on(x).trans.Nroi_s2n, on(x).trans.Nroi_s2s ]; %
        Nroi(x) = expt{x}.Nroi;
    end
end
xTrans = find(~isnan(sum(Ntrans,1)));
trans_frac = Ntrans(:,xTrans)./Nroi(xTrans);
transExpt = [expt{xTrans}];

figure;
bar(trans_frac', 'stacked')
legend({'a2a','a2n','a2s','n2a','n2n','n2s','s2a','s2n','s2s'})
axis square;
set(gca,'Xtick',1:numel(xTrans), 'XtickLabel',{transExpt.name}, 'TickLabelInterpreter','none')
ylabel('Fraction of ROIs');


%% Plot velocity and fluor
figure('Windowstate','maximized');
sp(1) = subplot(2,2,1);
sp(2) = subplot(2,2,3);
SP(1) = subplot(2,2,2);
SP(2) = subplot(2,2,4);
linkaxes(sp,'x');
linkaxes(SP,'x');
 LS = '-';

for x = 11 %xPrePost % x2Dcsd
    % Get pre-CSD data
    preVelocity = [onBout{x}(expt{x}.preRuns).velocity];
    preVelocity = vertcat(preVelocity{:});
    preLims = cumsum([onBout{x}(expt{x}.preRuns).Nscan]);
    preFluor = [onBout{x}(expt{x}.preRuns).fluor];
    preFluor = cat(1, preFluor{:});
    % Get post-CSD data
    postVelocity = [onBout{x}([expt{x}.acuteRuns, expt{x}.postRuns]).velocity];
    postVelocity = vertcat(postVelocity{:});
    postLims = cumsum([onBout{x}([expt{x}.acuteRuns, expt{x}.postRuns]).Nscan]);   
    postFluor = [onBout{x}([expt{x}.acuteRuns, expt{x}.postRuns]).fluor];
    postFluor = cat(1, postFluor{:});
    velocityLims = [min(vertcat(preVelocity, postVelocity)),  max(vertcat(preVelocity, postVelocity))];

    subplot(sp(2)); % %sp(1) = subplot(2,2,3);
    plot(preVelocity); hold on;
    for lim = 1:numel(preLims)
        line(preLims(lim)*[1,1], velocityLims, 'color','k', 'linestyle',LS)
    end
    set(gca, 'Xtick', preLims, 'XtickLabel',1:numel(preLims), 'TickDir','out')
    ylim(velocityLims)
    xlim([1,length(preVelocity)]);
    title(sprintf('%s: Pre-CSD', expt{x}.name ), 'Interpreter','none')
    ylabel('Velocity (cm/s)');
    MakeScaleBar([30*expt{x}.scanRate, 0], {get(gca,'Xlim'), get(gca,'Ylim')}, [0.1,0.1], [0.5,0], 'Xtext','30 s')

    subplot(SP(2)); % %SP(1) = subplot(2,2,4);
    plot(postVelocity); hold on;
    for lim = 1:numel(postLims)
        line(postLims(lim)*[1,1], velocityLims, 'color','k', 'linestyle',LS)
    end
    ylim(velocityLims)
    set(gca, 'Xtick', postLims, 'XtickLabel',1:numel(postLims), 'TickDir','out')
    xlim([1,length(postVelocity)]);
    title(sprintf('Post-CSD'), 'Interpreter','none')
    MakeScaleBar([30*expt{x}.scanRate, 0], {get(gca,'Xlim'), get(gca,'Ylim')}, [0.1,0.1], [0.5,0], 'Xtext','30 s')

    for roi = on(x).trans.roi_n2n % 1:on(x).Nroi
        subplot(sp(1)); % %sp(2) = subplot(2,2,1);
        cla;
        tempPreFluor = preFluor(:,roi);
        plot(tempPreFluor); hold on;
        tempPreFluorLims = [min(tempPreFluor), max(tempPreFluor)];
        for lim = 1:numel(preLims)
            line(preLims(lim)*[1,1], tempPreFluorLims, 'color','k', 'linestyle',LS)
        end
        ylabel('F');
        ylim(tempPreFluorLims); xlim([1,length(preVelocity)]);
        set(gca, 'Xtick', preLims, 'XtickLabel',1:numel(preLims), 'TickDir','out')
        title(sprintf('ROI %i', roi))

        subplot(SP(1)) %SP(1) = subplot(2,2,2);
        cla;
        tempPostFluor = postFluor(:,roi);
        tempPostFluorLims = [min(tempPostFluor), max(tempPostFluor)];
        plot(tempPostFluor); hold on;
        for lim = 1:numel(postLims)
            line(postLims(lim)*[1,1], tempPostFluorLims, 'color','k', 'linestyle',LS)
        end
        set(gca, 'Xtick', postLims, 'XtickLabel',1:numel(postLims), 'TickDir','out')
        ylim(tempPostFluorLims); xlim([1,length(postVelocity)]);
        pause;
    end
    %pause;
    clf;
end

%% 
figure;
axis square;
for x = xPrePost
    preDiff = squeeze(on(x).pre.fluor.norm_diff);
    postDiff = squeeze(on(x).post.fluor.norm_diff);
    preDiffMean = mean(preDiff,1, 'omitnan');
    postDiffMean = mean(postDiff,1, 'omitnan');
    
    [~, p_sens] = ttest2(preDiff, postDiff);
    on(x).trans.roi_sens = find(p_sens < 0.05 & postDiffMean-preDiffMean > 0);
    on(x).trans.Nroi_sens = numel(on(x).trans.roi_sens);
    on(x).trans.roi_desens = find(p_sens < 0.05 & postDiffMean-preDiffMean <= 0);
    on(x).trans.Nroi_desens = numel(on(x).trans.roi_desens);
    on(x).trans.roi_nochange = setdiff(1:expt{x}.Nroi, [on(x).trans.roi_sens, on(x).trans.roi_desens] );
    on(x).trans.Nroi_nochange = numel(on(x).trans.roi_nochange);

    for roi = on(x).trans.roi_sens %1:expt{x}.Nroi
        cla;
        JitterPlot(cell2padmat({preDiff(:,roi), postDiff(:,roi)}), 'new',false); 
        title(sprintf('p = %2.3f', p_sens(roi)))
        pause;
    end
end
