% Find fluorEvent associated with stillness and locomotion - run BoutResponse first - NEW VERSION
stillEpoch = cell(1,Nexpt); stillBin = cell(1,Nexpt); stillSumm = cell(1,Nexpt);
SAFE = cell(1,Nexpt); Nsafe = cell(1,Nexpt); % LAFE = cell(1,Nexpt); Nlafe = cell(1,Nexpt); CAFE = cell(1,Nexpt); Ncafe = cell(1,Nexpt);
tic; close all;
for x = xPresent
    % find long periods of stillness between locomotive bouts
    [stillEpoch{x}, stillBin{x}, stillSumm{x}] = BinStillEpochs(expt{x}, Tscan{x}, loco{x}, fluor{x}, deform{x}, csdBout{x}, 'show',false, 'criterion','state', 'bin',60);  %

    % Get event times of all events for each ROI
    roiEvent_Tstart = cell(expt{x}.Nroi, 1); roiEvent_magPeak = cell(expt{x}.Nroi, 1); roiEvent_dur = cell(expt{x}.Nroi, 1);
    for roi = 1:expt{x}.Nroi
        for runs = 1:expt{x}.Nruns
            tempRunEvents = [fluorEvent{x}{runs}{roi}];
            % Gather event properties of interest
            if numel(tempRunEvents) > 0
                roiEvent_Tstart{roi} = [roiEvent_Tstart{roi}, [tempRunEvents.Tstart]];
                roiEvent_dur{roi} = [roiEvent_dur{roi}, [tempRunEvents.dur]];
                roiEvent_magPeak{roi} = [roiEvent_magPeak{roi}, [tempRunEvents.magPeak]];
            end
        end
    end

    % Check each epoch for fluorEvents
    for runs = 1:expt{x}.Nruns %expt{x}.preRuns
        stillEpoch{x}(runs).Nevent = nan(expt{x}.Nroi,stillEpoch{x}(runs).Nepoch); stillEpoch{x}(runs).eventRate = nan(expt{x}.Nroi,stillEpoch{x}(runs).Nepoch); stillEpoch{x}(runs).eventPercent = nan(expt{x}.Nroi,stillEpoch{x}(runs).Nepoch);
        SAFE{x}{runs} = cell(expt{x}.Nroi,stillSumm{x}.Nepoch); % still-associated fluorescence events
        for ep = flip(1:stillEpoch{x}(runs).Nepoch)
            epochTimeLims = [stillEpoch{x}(runs).Tstart_trim(ep), stillEpoch{x}(runs).Tstop_trim(ep)];
            for roi = flip(find(NfluorEvent{x}(runs,:)))
                eventStartTime = [fluorEvent{x}{runs}{roi}.Tstart];
                epEvtInd = find(eventStartTime >= epochTimeLims(1) & eventStartTime < epochTimeLims(2));
                epEvent = fluorEvent{x}{runs}{roi}(epEvtInd);
                SAFE{x}{runs}{roi,ep} = epEvent; % save the events
                stillEpoch{x}(runs).Nevent(roi,ep) = numel(epEvtInd);
                stillEpoch{x}(runs).eventRate(roi,ep) = stillEpoch{x}(runs).Nevent(roi,ep)/stillEpoch{x}(runs).dur_trim(ep);
                stillEpoch{x}(runs).eventPercent(roi,ep) = numel(intersect(stillEpoch{x}(runs).scan_trim{ep}, vertcat(epEvent.runScan))')/stillEpoch{x}(runs).Nscan_trim(ep); % percent of epoch duration spent in event state
            end
        end
        stillEpoch{x}(runs).event_tot = sum(stillEpoch{x}(runs).Nevent, 2);
    end
    Nsafe{x} = cellfun(@numel, SAFE{x});
    % Check each bin for associated ROI events
    for bin = 1:numel(stillBin{x})
        binTimeLims = [stillBin{x}(bin).Tstart, stillBin{x}(bin).Tstop];
        for roi = 1:expt{x}.Nroi
            binEventInd = find(roiEvent_Tstart{roi} >= binTimeLims(1) & roiEvent_Tstart{roi} < binTimeLims(2));
            stillBin{x}(bin).Nevent(roi,1) = numel(binEventInd);
            stillBin{x}(bin).eventRate(roi,1) = stillBin{x}(bin).Nevent(roi,1)/stillBin{x}(bin).dur;
            stillBin{x}(bin).medPeakMag(roi,1) = median(roiEvent_magPeak{roi}(binEventInd));
            stillBin{x}(bin).medEventDur(roi,1) = median(roiEvent_dur{roi}(binEventInd));
        end
        %{
        for roi = flip(find(NfluorEvent{x}(runs,:)))
            eventStartTime = [fluorEvent{x}{runs}{roi}.Tstart];
            binEvtInd = find(eventStartTime >= binTimeLims(1) & eventStartTime < binTimeLims(2));
            binEvent = fluorEvent{x}{runs}{roi}(binEvtInd);
            stillBin{x}(bin).Nevent(roi,1) = numel(binEvtInd);
            stillBin{x}(bin).eventRate(roi,1) = stillBin{x}(bin).Nevent(roi)/stillBin{x}(bin).dur;
            stillBin{x}(bin).eventPercent(roi,1) = numel( intersect(stillBin{x}(bin).scan, vertcat(binEvent.runScan)) )/stillBin{x}(bin).Nscan; % percent of bin spent in event state
        end
        %}
    end
end

%% EVent rate ANOVA - compare invidiual ROI level binned event rates baseline (runs 1&2) vs  acute (run 3) vs post CSD (run 4)
minBin = 10;
Ntype = nan(0,3); % # of ROI neutral/suppressed/activated by CSD, with regards to event rate
typeFrac = nan(0,3);  % fraction of ROI neutral/suppressed/activated by CSD, with regards to event rate
k = 0;
roi_CSD_supp = cell(1, Nexpt); roi_CSD_neut = cell(1, Nexpt); roi_CSD_act = cell(1, Nexpt);
for x = [5,8,9,13,28,30,31,36,37,41]
    % Gather bin-level event statistiscs
    binEventRate = [stillBin{x}.eventRate];
    binNevent = [stillBin{x}.Nevent];
    binDur = [stillBin{x}.dur];
    % Group the data:  Which bins are pre/acute/post CSD?
    pre_bins = find([stillBin{x}.run] < expt{x}.csd);  Npre_bin = numel(pre_bins);
    acute_bins = find([stillBin{x}.run] == expt{x}.csd);  Nacute_bin = numel(acute_bins);
    post_bins = find([stillBin{x}.run] == expt{x}.csd+1);  Npost_bin = numel(post_bins);
    groupVec = [ones(1,numel(pre_bins)), 2*ones(1,numel(acute_bins)), 3*ones(1,numel(post_bins))];
    groupDur = [sum(binDur(pre_bins)), sum(binDur(acute_bins)), sum(binDur(post_bins))];
    Ngroup = size(groupDur,2);
    roi_groupRate = ([sum(binNevent(:,pre_bins), 2), sum(binNevent(:,acute_bins), 2), sum(binNevent(:,post_bins), 2)]./groupDur)'; % total events per period for each

    if Npre_bin >= minBin && Nacute_bin >= minBin && Npost_bin >= minBin
        k = k+1;
        % Individual ROI level
        pANOVA = nan(1,expt{x}.Nroi);
        mcResult = false(expt{x}.Nroi,3); % which differences were significant? 1) pre vs acute? , 2) pre vs post? 3) acute vs post?
        roi_group_rate = nan(expt{x}.Nroi, Ngroup);
        %roi_CSD_supp = []; roi_CSD_act = []; roi_CSD_neut = [];
        roi_CSD_supp_recov = []; roi_CSD_act_recov = [];
        for roi = 1:expt{x}.Nroi
            roiRate_pre = binEventRate(roi,pre_bins);
            roiRate_acute = binEventRate(roi,acute_bins);
            roiRate_post = binEventRate(roi,post_bins);
            %JitterPlot({roiRate_pre,roiRate_acute,roiRate_post}, 'Nboot',500)
            [pANOVA(roi),~, statANOVA] = anova1([roiRate_pre, roiRate_acute, roiRate_post], groupVec, 'off');
            if pANOVA(roi) < 0.05
                mcANOVA = multcompare(statANOVA, 'Display','on');
                if mcANOVA(1,6) < 0.05
                    mcResult(roi,:) = (mcANOVA(:,6) < 0.05)';
                    if mcANOVA(1,4) > 0 % is there an acute suppression?
                        roi_CSD_supp{x} = [roi_CSD_supp{x}, roi];
                        if mcANOVA(2,6) >= 0.05 || (mcANOVA(2,6) < 0.05 && mcANOVA(2,4) < 0) % did the suppressed unit return to or exceed, pre-CSD level?
                            roi_CSD_supp_recov = [roi_CSD_supp_recov, roi];
                        end
                    elseif mcANOVA(1,4) < 0 % is there an acute activation?
                        roi_CSD_act{x} = [roi_CSD_act{x}, roi];
                        if mcANOVA(2,6) >= 0.05 || (mcANOVA(2,6) < 0.05 && mcANOVA(2,4) < 0) % did the activated unit return to or below, pre-CSD level?
                            roi_CSD_act_recov = [roi_CSD_act_recov, roi];
                        end
                    end
                else
                    roi_CSD_neut{x} = [roi_CSD_neut{x}, roi];
                end
            else
                roi_CSD_neut{x} = [roi_CSD_neut{x}, roi];
            end
        end

        Ntype(k,3) = numel(roi_CSD_act{x});
        Ntype(k,2) = numel(roi_CSD_supp{x});
        Ntype(k,1) = numel(roi_CSD_neut{x});
        typeFrac(k,:) = Ntype(k,:)/expt{x}.Nroi;

        %numel(roi_CSD_supp_recov)
        %numel(roi_CSD_act_recov)
        %{
        figure;
        %plot(roi_groupRate); hold on;
        plot(roi_groupRate(:,roi_CSD_neut{x}), 'color',0.5*[1,1,1]); hold on;
        plot(roi_groupRate(:,roi_CSD_supp{x}), 'color','r'); hold on;
        plot(roi_groupRate(:,roi_CSD_act{x}), 'color','b'); hold on;
        pause;
        %}
    else
        fprintf('\nx = %i: %s, insufficient bins to run the ANOVA', x, expt{x}.name)
    end
end

StillEventRate_subtypes = figure('Units','inches', 'Position',[6, 5, 3, 2], 'Color','w');
%subplot(1,2,1);
JitterPlot(typeFrac, 0.2, 'paired',true,'new',false, 'Nboot',500);
axis square;
set(gca,'Xtick',1:3, 'XtickLabel',[], 'Ytick',0:0.2:1, 'YtickLabel',[], 'TickDir','out', 'TickLength',[0.01,0], 'box','off', 'Position',[0.02,0.02,0.96,0.96]); % {'Neutral','Suppressed','Activated'}
%ylabel('Fraction of ROIs')
figPath = sprintf('%sOngoing_subtypeBreakdown.pdf', figDir);
exportgraphics(StillEventRate_subtypes, figPath, 'Resolution',300); fprintf('\nSaved %s', figPath)

%% Heatmaps pre/acute/post CSD, sorted by ROI response subtype
close all;
StillEvents_heatmaps = figure('Units','inches', 'Position',[6, 5, 3, 2], 'Color','w');
for x = 5
    % Sort ROIs by spontaneous response to CSD
    roiTypeSort = [roi_CSD_supp{x}, roi_CSD_neut{x}, roi_CSD_act{x}];
    roiTypeTicks = cumsum([numel(roi_CSD_supp{x}), numel(roi_CSD_neut{x}), numel(roi_CSD_act{x})])+0.5;

    % Concatenate run data
    zTemp = [fluor{x}(1:4).z];
    zRunsCat = vertcat(zTemp.ROI); % roiTypeSort
    imagesc(zRunsCat(:,roiTypeSort)'); hold on;
    for runs = 2:4
        line(expt{x}.scanLims(runs)*[1,1], expt{x}.Nroi*[0,1], 'color','k')
    end

    for typ = 1:2
        line( [0,size(zRunsCat, 1)], roiTypeTicks(typ)*[1,1], 'color','k' )
    end

    set(gca,'Xtick',[], 'Ytick',[], 'TickDir','out', 'TickLength',[0.01,0], 'box','off', 'Position',[0.02,0.02,0.96,0.96])
end
caxis([-5,5])
colormap bluewhitered;

%% Event rates pre/post CSD - bins, absolute rate
close all;
figure;
for x = 28% [5,8,9,13,28,30,31,36,37,41]
    bin_pre = find([stillBin{x}.run] < expt{x}.csd);
    bin_rate = 60*mean([stillBin{x}.eventRate], 'omitnan');
    rate_pre = median(bin_rate(bin_pre));
    clf;
    line(0*[1,1], [0,max(bin_rate)], 'color','k'); hold on;
    line([-60,120], rate_pre*[1,1], 'linestyle','--','color','r');
    for bins = 1:numel(stillBin{x})
        %subplot(2,1,1);
        plot(([stillBin{x}.Tmid]-csdBout{x}(expt{x}.csd).Tstart)/60, bin_rate, '.' ); hold on;
        plot(([stillBin{x}.Tmid]-csdBout{x}(expt{x}.csd).Tstart)/60, bin_rate )
        %subplot(2,1,2); %cla;
        %plot([stillBin{x}.Tmid]-csdBout{x}(expt{x}.csd).Tstart, mean([stillBin{x}.medEventDur], 'omitnan'), '.' ); hold on; % [stillBin{x}.medPeakMag]
        %plot(stillEpoch{x}(runs).Tmid_trim, mean(stillEpoch{x}(runs).eventPercent,1 ,'omitnan') ); hold on;
    end
    ylabel('Event rate (events per minute)'); xlabel('Time post-CSD wave (mins)'); title(sprintf('x=%i: %s',x, expt{x}.name), 'Interpreter','none');
    pause;
end

%% Event rates pre/post CSD - bins, rate relative to pre-CSD baseline
discrete_width = 5; % minutes
Tlims = -120:discrete_width:150; 
Tcent = Tlims + discrete_width/2;
norm_rate_cell = cell(1,numel(Tlims));
close all;
figure;
for x = [5,8,9,13,28,30,31,36,37,41]
    bin_pre = find([stillBin{x}.run] < expt{x}.csd);
    bin_rate = 60*mean([stillBin{x}.eventRate], 'omitnan');
    med_rate_pre = median(bin_rate(bin_pre));
    bin_rate_norm = bin_rate/med_rate_pre;
    Tbin = ([stillBin{x}.Tmid]-csdBout{x}(expt{x}.csd).Tstart)/60;
    % Copy each normalized rate sample, into the corresponding (pooled) cell
    disc_cell = discretize(Tbin, Tlims);
    for bin = 1:numel(Tbin)
        norm_rate_cell{disc_cell(bin)} = [norm_rate_cell{disc_cell(bin)}; bin_rate_norm(bin)];
    end
    %{
    %clf;
    line(0*[1,1], [0,max(bin_rate_norm)], 'color','k'); hold on;
    line([-60,120], [1,1], 'linestyle','--','color','r');
    plot(Tbin, bin_rate_norm, '.' ); hold on;
    plot(Tbin, bin_rate_norm )
    title(sprintf('x=%i: %s',x, expt{x}.name), 'Interpreter','none');
    pause;
    %}
end

cellMean = cellfun(@mean, norm_rate_cell);
cellSEM = cellfun(@SEM , norm_rate_cell); % , 'UniformOutput',true

close all;
StillEventRate_pooled = figure('Units','inches', 'Position',[6, 5, 3, 2], 'Color','w');
errorbar( Tcent, cellMean, cellSEM, 'capsize',0 );
yLims = get(gca,'Ylim');
line([0,0], [0,yLims(2)], 'color','k','linestyle','--')
xlim([-60,60]); ylim([-Inf,Inf]);
ylabel('Event rate  during stillness (normalized to baseline)'); xlabel('Time post-CSD wave (mins)');
set(gca,'Xtick',-60:30:120, 'Ytick',0:0.5:1.5, 'TickDir','out', 'TickLength',[0.01,0], 'box','off', 'Position',[0.02,0.02,0.96,0.96])
figPath = sprintf('%sOngoing_PrePost_PooledErrorbar_120mins.pdf', figDir);
exportgraphics(StillEventRate_pooled, figPath, 'Resolution',300); fprintf('\nSaved %s', figPath)


%%
for x = xPresent
    %preScans = []; %[, stillEpoch{x}(2).scan_trim{:}];
    for roi = 1
        subplot(1,2,1)
        for runs = 1:2
            runEpochScans = stillEpoch{x}(runs).scan_trim{:};
            plot( fluor{x}(runs).dFF.ROI(runEpochScans,roi)); hold on; % Tscan{x}{runs}(runEpochScans),
        end
        ylabel('dF/F'); xlabel('Scans (concatenated)')

        subplot(1,2,2)
        for runs = 3:expt{x}.Nruns
            runEpochScans = stillEpoch{x}(runs).scan_trim{:};
            plot(Tscan{x}{runs}(runEpochScans), fluor{x}(runs).dFF.ROI(runEpochScans,roi)); hold on;
        end
        xlabel('Scans (concatenated)'); % ylabel('dF/F');
    end
end

%% Concatenated heatmaps of all still epochs, pre vs post CSD, sorted by ROI response type 
close all;
for x = [5,8,9,13,28,30,31,36,37,41] %xPresent
    dFFcat = [fluor{x}.dFF];
    dFFcat = vertcat(dFFcat.ROI);
    zCat = [fluor{x}.z]; zCat = vertcat(zCat.ROI);
    Nepoch = [stillEpoch{x}.Nepoch];
    roiTypeSort = [roi_CSD_supp{x}, roi_CSD_neut{x}, roi_CSD_act{x}];
    roiTypeTicks = cumsum([numel(roi_CSD_supp{x}), numel(roi_CSD_neut{x}), numel(roi_CSD_act{x})])+0.5;
    % Get the (concatenated) scans associated with still epochs for each run, pre and post csd
    preEpCatScans = [];  preEpCatLims = 0;
    postEpCatScans = []; postEpCatLims = 0;
    for runs = find(Nepoch(1:expt{x}.csd-1)>0)%
        tempEpCatScans = expt{x}.scanLims(runs) + [stillEpoch{x}(runs).scan_trim{:}];
        preEpCatScans = [preEpCatScans, tempEpCatScans];
        preEpCatLims = [preEpCatLims, preEpCatLims(end)+cumsum(cellfun(@length, stillEpoch{x}(runs).scan_trim))];
    end

    for runs = expt{x}.csd-1+find(Nepoch(expt{x}.csd:expt{x}.Nruns)>0) %
        tempEpCatScans = expt{x}.scanLims(runs) + [stillEpoch{x}(runs).scan_trim{:}];
        postEpCatScans = [postEpCatScans, tempEpCatScans];
        postEpCatLims = [postEpCatLims, postEpCatLims(end)+cumsum(cellfun(@length, stillEpoch{x}(runs).scan_trim))];
    end
    preEpCatLims = preEpCatLims+0.5; postEpCatLims = postEpCatLims+0.5;  % prevent line from covering the last sample

    figure('Name',expt{x}.name);
    sp(1) = subplot(1,2,1); cla;
    imagesc(zCat(preEpCatScans,roiTypeSort)'); caxis([-5,5]); hold on;
    for lim = 2:numel(preEpCatLims)
        line(preEpCatLims(lim)*[1,1], [0,expt{x}.Nroi+0.5], 'color','k' )
    end
    for typ = 1:2
        line( [0,numel(preEpCatScans)], roiTypeTicks(typ)*[1,1], 'color','k' )
    end
    title('Pre-CSD');
    sp(2) = subplot(1,2,2); cla;
    imagesc(zCat(postEpCatScans,roiTypeSort)'); caxis([-5,5])
    for lim = 2:numel(postEpCatLims)
        line(postEpCatLims(lim)*[1,1], [0,expt{x}.Nroi+0.5], 'color','k' )
    end
    for typ = 1:2
        line( [0,numel(postEpCatScans)], roiTypeTicks(typ)*[1,1], 'color','k' )
    end
    title('Post-CSD');
    colormap bluewhitered
    pause;
    %{
    % Individual ROI traces
    for roi = 1:expt{x}.Nroi
        sp(1) = subplot(1,2,1); cla;
        plot(dFFcat(preEpCatScans,roi)); hold on;
        axis tight;
        sp(2) = subplot(1,2,2); cla;
        plot(dFFcat(postEpCatScans,roi)); hold on;
        axis tight;
        linkaxes(sp,'y')
        pause;
    end
    %}
end

%% Make movies contrasting spontaneous activity pre/acute/post
movParam.dir = 'D:\MATLAB\Figures\CSD figures\'; mkdir(movParam.dir); % 'D:\MATLAB\LevyLab\Figures\CSD\Movies\';
movParam.fmtSpec = '%2.2f';
movParam.Tperi = [0,0]; % time before/after bout to show
movParam.zProj = [];
movParam.binT = 16;
movParam.displayPct = [5,99.9];
movParam.sbx = [];
movParam.regType =  'affine'; %'raw'; % % _axons
movParam.boutType = 'still';
movParam.edges = [70,190,40,40]; % [80,80,20,20]; %[80,90,60,110]; %
movParam.aviRate = 20; % frames per second
movParam.Toffset = 0;
movParam.level = 'ind'; %'both'; %'none'; %
for x = 28%xPresent 
    %movParam.edges = segParams{x}.edges;
    movParam.zProj = segParams{x}.zProj; %segPlanes;
    if expt{x}.Nplane == 1
        movParam.sourceSbx = 'sbx_affine'; % %'sbxz'; %  
    else
        movParam.sourceSbx = 'sbx_interp';
    end
    %movParam.scalebar = [];
    movParam.scalebar = MakeScaleBar( round(expt{x}.umPerPixel*[50,0]), {[0,expt{x}.Ncol-segParams{x}.edges(1)-segParams{x}.edges(2)]+0.5, [0,expt{x}.Nrow-segParams{x}.edges(3)-segParams{x}.edges(4)]+0.5},...
        [0.1,0.95], [0,0], 'label',false, 'color','w', 'show',false );
     % smuggle stillBins into WriteBoutMovies by using bin # as 'run'
    [boutStack{x}, Tbout{x}, boutSpeed{x}, storeFrames{x}] = WriteBoutMovies(expt{x}, catInfo{x}, Tscan{x}, loco{x}, stillBin{x}, movParam, 'run',[43,44,46]); % [14,30,51]
end

