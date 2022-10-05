minStillDur = 120; 
stillClip = 30; % clip the first few seconds of stillness to guard against residual locomotion-associated activity
stillChunkDur = 60;
maxMissingFrac = 0.1; % exclude epochs of stillness where the median ROI is missing at least this much data
minFluorDur = 1; 
minScaleDur = 0.2;
minZ = 2.5; 
minDFF = 0.05;
minROIdur = 50;
Tpad = [2,2];
NeventMin = 5;
acuteLimit = 15;
postLimit = 60;
Ndim = 1;
show = false; %  true; %   
stillEpoch = cell(1,Nexpt); stillSumm = cell(1,Nexpt); fluorEvent = cell(1,Nexpt); scaleEvent = cell(1,Nexpt); 
if show 
    close all; 
    figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'color','w'); 
    opt = {[0.03,0.02], [0.07,0.04], [0.05,0.05]};  % {[vert, horz], [bottom, top], [left, right] }
    clearvars sp;
    for p = flip(1:3), sp(p) = subtightplot(4,1,p,opt{:});   end
end
for x = xPresent
    % Find sufficiently long inter-bout periods of stillness
    stillClipScan = round(stillClip*expt(x).scanRate);
    stillChunkNscan = round(stillChunkDur*expt(x).scanRate);
    S = 0;
    for run = 1:expt(x).Nruns
        tempStillCC = bwconncomp(loco{x}(run).bout ~= 1); %bwconncomp(loco{x}(run).stateDown == 1);
        tempStillDur = cellfun(@numel, tempStillCC.PixelIdxList)/expt(x).scanRate;
        for s = find(tempStillDur >= minStillDur)
            % Break each epoch of stillness into epochs of predetermined size
            tempStillScanRange = tempStillCC.PixelIdxList{s}([stillClipScan,end-stillClipScan]);
            [tempChunkLims, ~, tempChunkSize] = MakeChunkLims( tempStillScanRange(1), tempStillScanRange(end), tempStillScanRange(end), 'size',stillChunkNscan);
            % Absorb remainder into last chunk
            if tempChunkSize(end) < stillChunkNscan
                tempChunkLims(end-1,2) = tempChunkLims(end,2); 
                tempChunkLims(end,:) = [];
                tempChunkSize = diff(tempChunkLims,1,2)+1;
            end
            for c = find(tempChunkSize >= stillChunkNscan)'
                S = S+1;
                stillEpoch{x}(S).run = run;
                stillEpoch{x}(S).scan = tempChunkLims(c,1):tempChunkLims(c,2); %tempStillCC.PixelIdxList{s}(stillClipScan:end-stillClipScan);
                stillEpoch{x}(S).Nscan = numel(stillEpoch{x}(S).scan);
                stillEpoch{x}(S).T = Tscan{x}{run}(stillEpoch{x}(S).scan); % -Tscan{x}{run}(stillEpoch{x}(S).scan(1))
                stillEpoch{x}(S).Tstart = stillEpoch{x}(S).T(1);
                stillEpoch{x}(S).Tmid = median(stillEpoch{x}(S).T);
                stillEpoch{x}(S).Tstop = stillEpoch{x}(S).T(end);
                stillEpoch{x}(S).dur = stillEpoch{x}(S).T(end) - stillEpoch{x}(S).T(1);
                % get fluor/deformation data during each epoch
                stillEpoch{x}(S).Froi = median( fluor{x}(run).Froi.ROI(stillEpoch{x}(S).scan,:), 1, 'omitnan' );
                stillEpoch{x}(S).Fnp = median( fluor{x}(run).Fnp.ROI(stillEpoch{x}(S).scan,:), 1, 'omitnan' );
                stillEpoch{x}(S).Fo = median( fluor{x}(run).Fo.ROI(stillEpoch{x}(S).scan,:), 1, 'omitnan' );
                stillEpoch{x}(S).dff = fluor{x}(run).dFF.ROI(stillEpoch{x}(S).scan,:);
                if expt(x).Nplane > 1
                    stillEpoch{x}(S).scale = mean(deform{x}(run).scaleMag(stillEpoch{x}(S).scan,4:end-3), 2, 'omitnan');
                else
                    stillEpoch{x}(S).scale = deform{x}(run).scaleMag(stillEpoch{x}(S).scan,:);
                end
                % account for missing (NaN) data
                stillEpoch{x}(S).Nmissing = sum(isnan( stillEpoch{x}(S).dff), 1); 
                stillEpoch{x}(S).missingFrac = stillEpoch{x}(S).Nmissing/stillEpoch{x}(S).Nscan;
                stillEpoch{x}(S).missingFracMed = median(stillEpoch{x}(S).missingFrac);
                stillEpoch{x}(S).roiDur = stillEpoch{x}(S).dur*(1-stillEpoch{x}(S).missingFrac);
                stillEpoch{x}(S).roiDurMed = median(stillEpoch{x}(S).roiDur);
            end
        end
    end
    stillEpoch{x}([stillEpoch{x}.missingFracMed] > maxMissingFrac) = []; % Remove any epochs where fluor data is mostly bad  [stillEpoch{x}.roiDurMed] < minROIdur
    % Divide up pre-CSD, acute and post-CSD phase epochs
    if ~isnan(expt(x).csd) 
        stillEpoch{x}([stillEpoch{x}.run] == expt(x).csd & [stillEpoch{x}.Tstart] < csdBout{x}(expt(x).csd).Tstop) = []; % Remove any epochs that include the CSD wave itself
        stillSumm{x}.Trel = ([stillEpoch{x}.Tmid]-csdBout{x}(expt(x).csd).Tstart)/60; % epoch timing, relative to CSD, in minutes
    else
        stillSumm{x}.Trel = ([stillEpoch{x}.Tmid]-Tscan{x}{end}(end))/60; %Trel = [stillEpoch{x}.Tmid]/60; % epoch timing, relative to CSD, in minutes
    end
    stillSumm{x}.sPre = find(stillSumm{x}.Trel < 0); %stillSumm{x}.sPre = find([stillEpoch{x}.Tmid] < csdBout{x}(expt(x).csd).Tstart);
    stillSumm{x}.sAcute = find(stillSumm{x}.Trel >= 0 & stillSumm{x}.Trel < acuteLimit);
    stillSumm{x}.sPost = find(stillSumm{x}.Trel >= acuteLimit);
    
    % Summarize epoch-level data
    stillSumm{x}.Nepoch = numel(stillEpoch{x});
    stillSumm{x}.dur = [stillEpoch{x}.dur];
    stillSumm{x}.TotDur = sum(stillSumm{x}.dur); 
    stillSumm{x}.PreDur = sum(stillSumm{x}.dur(stillSumm{x}.sPre), 2);
    stillSumm{x}.AcuteDur = sum(stillSumm{x}.dur(stillSumm{x}.sAcute), 2);
    stillSumm{x}.PostDur = sum(stillSumm{x}.dur(stillSumm{x}.sPost), 2);
    stillSumm{x}.roiDur = vertcat(stillEpoch{x}.roiDur);
    stillSumm{x}.roiTotDur = sum(stillSumm{x}.roiDur, 1);
    stillSumm{x}.roiPreDur = sum(stillSumm{x}.roiDur(stillSumm{x}.sPre,:), 1);
    stillSumm{x}.roiAcuteDur = sum(stillSumm{x}.roiDur(stillSumm{x}.sAcute,:), 1);
    stillSumm{x}.roiPostDur = sum(stillSumm{x}.roiDur(stillSumm{x}.sPost,:), 1);
    stillFluorPoolZ = normalize(vertcat(stillEpoch{x}(:).dff));
    stillScalePoolZ = normalize(vertcat(stillEpoch{x}(:).scale));
    % Unpool the normalized data
    for s = 1:stillSumm{x}.Nepoch
        tempScan = 1:stillEpoch{x}(s).Nscan;
        stillEpoch{x}(s).z = stillFluorPoolZ(tempScan,:);
        stillEpoch{x}(s).zScale = stillScalePoolZ(tempScan,:);
        stillFluorPoolZ(tempScan,:) = [];
        stillScalePoolZ(tempScan,:) = [];
    end

    % Get distribution of deformation variables pre/acute/post-CSD
    if ~isempty(stillSumm{x}.sPre)
        tempPreDeform = vertcat( stillEpoch{x}(stillSumm{x}.sPre).scale ); % Pre
        [stillSumm{x}.scale.CDFf{1}, stillSumm{x}.scale.CDFx{1}] = ecdf(tempPreDeform);
    end
    if ~isempty(stillSumm{x}.sAcute)
        tempAcuteDeform = vertcat( stillEpoch{x}(stillSumm{x}.sAcute).scale ); % Acute
        [stillSumm{x}.scale.CDFf{2}, stillSumm{x}.scale.CDFx{2}] = ecdf(tempAcuteDeform);
    end
    if ~isempty(stillSumm{x}.sPost)
        tempPostDeform = vertcat( stillEpoch{x}(stillSumm{x}.sPost).scale ); % Post
        [stillSumm{x}.scale.CDFf{3}, stillSumm{x}.scale.CDFx{3}] = ecdf(tempPostDeform);
    end    
    % Get fluor and deformation events within each epoch of stillness
    padScans = round(expt(x).scanRate*Tpad);
    for s = 1:stillSumm{x}.Nepoch      
        % Look for fluor events
        tempNormFluor = stillEpoch{x}(s).z./max(stillEpoch{x}(s).z, [], 1, 'omitnan');
        tempNormFluor(tempNormFluor <= 0) = 0;
        %{
        if show
            subplot(sp(1)); %sp(1) = subplot(4,1,1); 
            imagesc(stillEpoch{x}(s).dff'); tempPos = get(gca,'Position');
            CB(1) = colorbar; CB(1).Label.String = 'dF/F';
            set(gca,'Xtick',[], 'Position',tempPos, 'Ytick',[1,expt(x).Nroi], 'TickDir','out');
            title( sprintf('%s: [x, s, run] = [%i, %i, %i]', expt(x).name, x, s, stillEpoch{x}(s).run), 'Interpreter','none');
            ylabel('ROI');

            subplot(sp(2)); %sp(2) = subplot(4,1,2); 
            imagesc(stillEpoch{x}(s).z'); tempPos = get(gca,'Position');
            CB(2) = colorbar; CB(2).Label.String = 'Z (still epochs only)';
            set(gca,'Xtick',[], 'Position',tempPos, 'TickDir','out'); 
            ylabel('ROI');

            subplot(sp(3)); %sp(3) = subplot(4,1,3); 
            imagesc(tempNormFluor'); tempPos = get(gca,'Position');
            CB(3) = colorbar; CB(3).Label.String = 'Normalized';
            set(gca,'Xtick',[], 'Position',tempPos, 'TickDir','out');
            ylabel('ROI');
            linkaxes(sp,'xy');
            impixelinfo;
        end
        %}
        for roi = 1:expt(x).Nroi
            %{
            if show
                subplot(sp(3)); %
                set(gca, 'Ytick',roi, 'YtickLabel',num2str(roi));
                
                subtightplot(4,1,4,opt{:}); cla; %subplot(sp(4));
                plot( stillEpoch{x}(s).T-stillEpoch{x}(s).T(1), stillEpoch{x}(s).dff(:,roi) ); hold on;
                line([0,stillEpoch{x}(s).T(end)-stillEpoch{x}(s).T(1)], [0,0], 'lineStyle','--', 'color','k');
                line([0,stillEpoch{x}(s).T(end)-stillEpoch{x}(s).T(1)], minDFF*[1,1], 'lineStyle','-.', 'color','k');
                xlabel('Time still (s)'); ylabel('dF/F_0');
                xlim([-Inf,inf]); ylim([-0.1, Inf]); %
            end
            %}
            tempNormCC = bwconncomp(tempNormFluor(:,roi) > 0);
            putativeEventDur = cellfun(@numel, tempNormCC.PixelIdxList)/expt(x).scanRate;
            E = 0;
            for e = find(putativeEventDur > minFluorDur)
                evtDFF = stillEpoch{x}(s).dff( tempNormCC.PixelIdxList{e}, roi);
                evtZ = stillEpoch{x}(s).z( tempNormCC.PixelIdxList{e}, roi );
                peakDFF = max(evtDFF);
                peakZ = max(evtZ);
                if peakZ >= minZ && peakDFF >= minDFF
                    E = E+1;
                    tempScan = tempNormCC.PixelIdxList{e}';
                    tempScan = tempScan(1)-padScans(1):tempScan(end)+padScans(2);
                    tempScan(tempScan < 1) = []; 
                    tempScan(tempScan > stillEpoch{x}(s).Nscan) = []; % expt(x).Nscan(run)
                    fluorEvent{x}(s,roi).scan{E} = tempScan; % scans in reference to this still epoch, not the full data
                    fluorEvent{x}(s,roi).preScan{E} = find( fluorEvent{x}(s,roi).scan{E} < tempNormCC.PixelIdxList{e}(1));
                    fluorEvent{x}(s,roi).boutScan{E} = find(ismember(fluorEvent{x}(s,roi).scan{E}, tempNormCC.PixelIdxList{e} ));
                    fluorEvent{x}(s,roi).postScan{E} = find( fluorEvent{x}(s,roi).scan{E} > tempNormCC.PixelIdxList{e}(end));
                    fluorEvent{x}(s,roi).T{E} = stillEpoch{x}(s).T(tempScan); % Tscan{x}{run}(roiEvent{x}(s,roi).scan{E});
                    fluorEvent{x}(s,roi).Tstart(E) = stillEpoch{x}(s).T(tempNormCC.PixelIdxList{e}(1));
                    fluorEvent{x}(s,roi).Tstop(E) = stillEpoch{x}(s).T(tempNormCC.PixelIdxList{e}(end));
                    fluorEvent{x}(s,roi).dur(E) = putativeEventDur(e);
                    fluorEvent{x}(s,roi).dffPeak(E) = peakDFF;
                    fluorEvent{x}(s,roi).zpeak(E) = peakZ;
                    %{
                    if show
                        try
                            plot( stillEpoch{x}(s).T(roiEvent{x}(s,roi).scan{E})-stillEpoch{x}(s).T(1), stillEpoch{x}(s).dff(roiEvent{x}(s,roi).scan{E},roi), '.','MarkerSize',7 ); % k
                            hold all;
                        catch
                            fprintf('\n[x, run, s, ROI, event] = [%i, %i, %i, %i, %i]: edge issue', x, run, s, roi, E);
                        end
                        title( sprintf('[ROI, event] = [%i, %i]', roi, E), 'Interpreter','none')
                        %pause;
                    end
                    %}
                end
            end
            if stillEpoch{x}(s).missingFrac(roi) < maxMissingFrac %stillSumm{x}.roiDur(s,roi) >= minROIdur
                fluorEvent{x}(s,roi).Nevent = E;
                fluorEvent{x}(s,roi).eventRate = 60*fluorEvent{x}(s,roi).Nevent/stillEpoch{x}(s).roiDur(roi); % events/min
            else
                fprintf('\n[x,s,roi] = [%i, %i, %i]: missing %2.2f percent of epoch data', x, s, roi, 100*stillEpoch{x}(s).missingFrac(roi)) %  stillSumm{x}.roiDur(s,roi)
                fluorEvent{x}(s,roi).Nevent = NaN; 
                fluorEvent{x}(s,roi).eventRate = NaN;
            end
            %pause;
        end
        %}

        % Look for DEFORMATION events within each epoch of stillness
        tempNormScale = stillEpoch{x}(s).zScale./max(stillEpoch{x}(s).zScale, [], 1, 'omitnan');
        tempNormScale(tempNormScale <= 0) = 0;
        if show
            subplot(sp(1)); %sp(1) = subplot(4,1,1); 
            plot(stillEpoch{x}(s).T-stillEpoch{x}(s).T(1), stillEpoch{x}(s).scale)
            %imagesc(stillEpoch{x}(s).scale'); tempPos = get(gca,'Position');
            %CB(1) = colorbar; CB(1).Label.String = 'Scale';
            set(gca,'XtickLabel',[], 'TickDir','out'); %set(gca,'Xtick',[], 'Position',tempPos, 'Ytick',[1,2], 'TickDir','out');
            ylabel('Scale mag (\mum)'); %ylabel('Axis');
            title( sprintf('%s: [x, s, run] = [%i, %i, %i]', expt(x).name, x, s, stillEpoch{x}(s).run), 'Interpreter','none');
            

            subplot(sp(2)); %sp(2) = subplot(4,1,2); 
            plot(stillEpoch{x}(s).T-stillEpoch{x}(s).T(1), stillEpoch{x}(s).zScale); hold on;
            line([0,stillEpoch{x}(s).T(end)-stillEpoch{x}(s).T(1)], minZ*[1,1], 'color','r', 'lineStyle','--'); hold off;
            %imagesc(stillEpoch{x}(s).zScale'); tempPos = get(gca,'Position');
            %CB(2) = colorbar; CB(2).Label.String = 'Z (still epochs only)';
            set(gca,'XtickLabel',[], 'TickDir','out'); %set(gca,'Xtick',[], 'Position',tempPos, 'TickDir','out', 'Ytick',[1,2]); 
            ylabel('Scale mag (z-score)');%ylabel('Axis');

            subplot(sp(3)); %sp(3) = subplot(4,1,3); 
            %imagesc(tempNormScale'); tempPos = get(gca,'Position');
            %CB(3) = colorbar; CB(3).Label.String = 'Normalized';
            %set(gca,'Xtick',[], 'Position',tempPos, 'TickDir','out');
            plot(stillEpoch{x}(s).T-stillEpoch{x}(s).T(1), tempNormScale);
            set(gca,'XtickLabel',[], 'TickDir','out'); %set(gca,'Xtick',[], 'Position',tempPos, 'TickDir','out');
            ylabel('Normalized Scale'); %ylabel('Axis');
            linkaxes(sp,'x');
            %impixelinfo;
            %pause;
        end
        
        for d = 1:Ndim%:2
            if show
                %subplot(sp(3)); set(gca, 'Ytick',d, 'YtickLabel',num2str(d));
                
                subtightplot(4,1,4,opt{:}); cla; %subplot(sp(4));
                plot( stillEpoch{x}(s).T-stillEpoch{x}(s).T(1), stillEpoch{x}(s).scale(:,d) ); hold on;
                line([0,stillEpoch{x}(s).T(end)-stillEpoch{x}(s).T(1)], [0,0], 'lineStyle','--', 'color','k');
                %line([0,stillEpoch{x}(s).T(end)-stillEpoch{x}(s).T(1)], minDFF*[1,1], 'lineStyle','-.', 'color','k');
                xlabel('Epoch Time (s)'); ylabel('Scale (\mum)');
                xlim([-Inf,inf]); ylim([-0.1, Inf]); %
                %pause;
            end
            tempNormCC = bwconncomp(tempNormScale(:,d) > 0);
            putativeEventDur = cellfun(@numel, tempNormCC.PixelIdxList)/expt(x).scanRate;
            E = 0;
            for e = find(putativeEventDur > minScaleDur)
                evtScale = stillEpoch{x}(s).scale(tempNormCC.PixelIdxList{e}, d);
                peakScale = max(evtScale);
                evtScaleZ = stillEpoch{x}(s).zScale(tempNormCC.PixelIdxList{e}, d );
                peakScaleZ = max(evtScaleZ);
                if peakScaleZ >= minZ % && peakScale >= minDFF
                    E = E+1;
                    tempScan = tempNormCC.PixelIdxList{e}';
                    tempScan = tempScan(1)-padScans(1):tempScan(end)+padScans(2);
                    tempScan(tempScan < 1) = []; 
                    tempScan(tempScan > stillEpoch{x}(s).Nscan) = []; % expt(x).Nscan(run)
                    scaleEvent{x}(s,d).scan{E} = tempScan;
                    scaleEvent{x}(s,d).preScan{E} = find( scaleEvent{x}(s,d).scan{E} < tempNormCC.PixelIdxList{e}(1));
                    scaleEvent{x}(s,d).boutScan{E} = find(ismember(scaleEvent{x}(s,d).scan{E}, tempNormCC.PixelIdxList{e} ));
                    scaleEvent{x}(s,d).postScan{E} = find( scaleEvent{x}(s,d).scan{E} > tempNormCC.PixelIdxList{e}(end));
                    scaleEvent{x}(s,d).T{E} = stillEpoch{x}(s).T(tempScan); % Tscan{x}{run}(scaleEvent{x}(s,d).scan{E});
                    scaleEvent{x}(s,d).Tstart(E) = stillEpoch{x}(s).T(tempNormCC.PixelIdxList{e}(1));
                    scaleEvent{x}(s,d).Tstop(E) = stillEpoch{x}(s).T(tempNormCC.PixelIdxList{e}(end));
                    scaleEvent{x}(s,d).dur(E) = putativeEventDur(e);
                    scaleEvent{x}(s,d).scalePeak(E) = peakScale;
                    scaleEvent{x}(s,d).zpeak(E) = peakScaleZ;
                    if show
                        try
                            plot( stillEpoch{x}(s).T(scaleEvent{x}(s,d).scan{E})-stillEpoch{x}(s).T(1), stillEpoch{x}(s).scale(scaleEvent{x}(s,d).scan{E},d), '.','MarkerSize',7 ); % k
                            hold all;
                        catch
                            fprintf('\n[x, run, s, axis, event] = [%i, %i, %i, %i, %i]: edge issue', x, run, s, d, E);
                        end
                        title( sprintf('[axis, event] = [%i, %i]', d, E), 'Interpreter','none')
                        pause;
                    end
                end
            end
            scaleEvent{x}(s,d).Nevent = E;
            scaleEvent{x}(s,d).eventRate = 60*scaleEvent{x}(s,d).Nevent/stillEpoch{x}(s).dur; % events/min
            %pause;
        end
    end
      
    % Summarize deformation events
    stillSumm{x}.scale.Nevent = reshape( [scaleEvent{x}.Nevent], stillSumm{x}.Nepoch, Ndim ); % # events for each epoch of stillness and axis 
    stillSumm{x}.scale.eventRate = reshape( [scaleEvent{x}.eventRate], stillSumm{x}.Nepoch, Ndim ); % events per time for each epoch and axis
    stillSumm{x}.scale.epochRate = 60*(sum(stillSumm{x}.scale.Nevent, 2)./[stillEpoch{x}.dur]')/Ndim; % total events (across ROI), per epoch duration per axis
    stillSumm{x}.scale.TotEvent = sum(stillSumm{x}.scale.Nevent, 1);
    stillSumm{x}.scale.TotRate = 60*stillSumm{x}.scale.TotEvent./stillSumm{x}.TotDur; % total events (across epochs), per total duration per axis
    stillSumm{x}.scale.PreRate = 60*sum(stillSumm{x}.scale.Nevent(stillSumm{x}.sPre,:), 1)./stillSumm{x}.PreDur;
    stillSumm{x}.scale.AcuteRate = 60*sum(stillSumm{x}.scale.Nevent(stillSumm{x}.sAcute,:), 1)./stillSumm{x}.AcuteDur;
    stillSumm{x}.scale.PostRate = 60*sum(stillSumm{x}.scale.Nevent(stillSumm{x}.sPost,:), 1)./stillSumm{x}.PostDur;
    
    % Summarize ROI base fluor
    stillSumm{x}.fluor.Fo = vertcat(stillEpoch{x}.Fo);
    stillSumm{x}.fluor.Froi = vertcat(stillEpoch{x}.Froi);
    stillSumm{x}.fluor.Fnp = vertcat(stillEpoch{x}.Fnp);
    
    % Summarize ROI fluor events
    stillSumm{x}.fluor.Nevent = reshape( [fluorEvent{x}.Nevent], stillSumm{x}.Nepoch, expt(x).Nroi ); % # events for each epoch of stillness and ROI 
    stillSumm{x}.fluor.eventRate = reshape( [fluorEvent{x}.eventRate], stillSumm{x}.Nepoch, expt(x).Nroi ); % events per time for each epoch and ROI
    stillSumm{x}.fluor.epochRate = 60*(sum(stillSumm{x}.fluor.Nevent, 2, 'omitnan')./[stillEpoch{x}.dur]')/expt(x).Nroi; % total events (across ROI), per epoch duration per ROI % 
    stillSumm{x}.fluor.TotEvent = sum(stillSumm{x}.fluor.Nevent, 1);
    stillSumm{x}.fluor.TotRate = 60*sum(stillSumm{x}.fluor.Nevent, 1)./stillSumm{x}.roiTotDur; % total events (across epochs), per total duration per ROI
    stillSumm{x}.fluor.PreRate = 60*sum(stillSumm{x}.fluor.Nevent(stillSumm{x}.sPre,:), 1)./stillSumm{x}.roiPreDur;
    stillSumm{x}.fluor.AcuteRate = 60*sum(stillSumm{x}.fluor.Nevent(stillSumm{x}.sAcute,:), 1)./stillSumm{x}.roiAcuteDur;
    stillSumm{x}.fluor.PostRate = 60*sum(stillSumm{x}.fluor.Nevent(stillSumm{x}.sPost,:), 1)./stillSumm{x}.roiPostDur;
    stillSumm{x}.fluor.rActive = find(stillSumm{x}.fluor.TotRate > 0);
    stillSumm{x}.fluor.rInactive = find(stillSumm{x}.fluor.TotRate == 0); % Which ROI never had a spontaneous event?
    if ~isempty(stillSumm{x}.sPost)
        [stillSumm{x}.fluor.PostRateCDF(:,2), stillSumm{x}.fluor.PostRateCDF(:,1)] = ecdf(stillSumm{x}.fluor.PostRate);
        [stillSumm{x}.scale.PostRateCDF(:,2), stillSumm{x}.scale.PostRateCDF(:,1)] = ecdf(stillSumm{x}.scale.PostRate);
    else
        stillSumm{x}.fluor.PostRateCDF = nan(1,2);
    end
    if ~isempty(stillSumm{x}.sAcute)
        [stillSumm{x}.fluor.AcuteRateCDF(:,2), stillSumm{x}.fluor.AcuteRateCDF(:,1)] = ecdf(stillSumm{x}.fluor.AcuteRate);
        [stillSumm{x}.scale.AcuteRateCDF(:,2), stillSumm{x}.scale.AcuteRateCDF(:,1)] = ecdf(stillSumm{x}.scale.AcuteRate);
    else
        stillSumm{x}.fluor.AcuteRateCDF = nan(1,2);
        stillSumm{x}.scale.AcuteRateCDF = nan(1,2);
    end
    if ~isempty(stillSumm{x}.sPre)
        [stillSumm{x}.fluor.PreRateCDF(:,2), stillSumm{x}.fluor.PreRateCDF(:,1)] = ecdf(stillSumm{x}.fluor.PreRate);
        [stillSumm{x}.scale.PreRateCDF(:,2), stillSumm{x}.scale.PreRateCDF(:,1)] = ecdf(stillSumm{x}.scale.PreRate);
    else
        stillSumm{x}.fluor.PreRateCDF = nan(1,2);
        stillSumm{x}.scale.PreRateCDF = nan(1,2);
    end
    % Pool event parameters by case: pre-CSD, acute, or post-CSD
    stillSumm{x}.fluor.dffPre = cell(1,expt(x).Nroi); stillSumm{x}.fluor.dffAcute = cell(1,expt(x).Nroi); stillSumm{x}.fluor.dffPost = cell(1,expt(x).Nroi);
    stillSumm{x}.fluor.durPre = cell(1,expt(x).Nroi); stillSumm{x}.fluor.durAcute = cell(1,expt(x).Nroi); stillSumm{x}.fluor.durPost = cell(1,expt(x).Nroi);
    stillSumm{x}.scale.durPre = cell(1,expt(x).Nroi); stillSumm{x}.scale.durAcute = cell(1,expt(x).Nroi); stillSumm{x}.scale.durPost = cell(1,expt(x).Nroi);
    for roi = 1:expt(x).Nroi
        for s = find(stillSumm{x}.fluor.Nevent(:,roi))' % stillSumm{x}.sPre
            if ismember(s,stillSumm{x}.sPre)
                stillSumm{x}.fluor.dffPre{roi} = [stillSumm{x}.fluor.dffPre{roi}, fluorEvent{x}(s,roi).dffPeak];
                stillSumm{x}.fluor.durPre{roi} = [stillSumm{x}.fluor.durPre{roi}, fluorEvent{x}(s,roi).dur];
            elseif ismember(s,stillSumm{x}.sAcute)
                stillSumm{x}.fluor.dffAcute{roi} = [stillSumm{x}.fluor.dffAcute{roi}, fluorEvent{x}(s,roi).dffPeak];
                stillSumm{x}.fluor.durAcute{roi} = [stillSumm{x}.fluor.durAcute{roi}, fluorEvent{x}(s,roi).dur];
            elseif ismember(s,stillSumm{x}.sPost)
                stillSumm{x}.fluor.dffPost{roi} = [stillSumm{x}.fluor.dffPost{roi}, fluorEvent{x}(s,roi).dffPeak];
                stillSumm{x}.fluor.durPost{roi} = [stillSumm{x}.fluor.durPost{roi}, fluorEvent{x}(s,roi).dur];
            end
        end
    end
    Ncase = [cellfun(@numel, stillSumm{x}.fluor.dffPre)', cellfun(@numel, stillSumm{x}.fluor.dffAcute)', cellfun(@numel, stillSumm{x}.fluor.dffPost)']'; % # of events pre/acute/post-CSD
    stillSumm{x}.fluor.dffMed = [cellfun(@median, stillSumm{x}.fluor.dffPre)', cellfun(@median, stillSumm{x}.fluor.dffAcute)', cellfun(@median, stillSumm{x}.fluor.dffPost)']';
    stillSumm{x}.fluor.dffMed(Ncase < NeventMin) = NaN;
    stillSumm{x}.fluor.durMed = [cellfun(@median, stillSumm{x}.fluor.durPre)', cellfun(@median, stillSumm{x}.fluor.durAcute)', cellfun(@median, stillSumm{x}.fluor.durPost)']';
    stillSumm{x}.fluor.durMed(Ncase < NeventMin) = NaN;

    Ncase = [cellfun(@numel, stillSumm{x}.scale.durPre)', cellfun(@numel, stillSumm{x}.scale.durAcute)', cellfun(@numel, stillSumm{x}.scale.durPost)']'; % # of events pre/acute/post-CSD
    for d = 1:Ndim
        for s = find(stillSumm{x}.scale.Nevent(:,d))' % stillSumm{x}.sPre
            if ismember(s,stillSumm{x}.sPre)
                %stillSumm{x}.scale.dffPre{ax} = [stillSumm{x}.scale.dffPre{ax}, scaleEvent{x}(s,ax).dffPeak];
                stillSumm{x}.scale.durPre{d} = [stillSumm{x}.scale.durPre{d}, scaleEvent{x}(s,d).dur];
            elseif ismember(s,stillSumm{x}.sAcute)
                %stillSumm{x}.scale.dffAcute{ax} = [stillSumm{x}.scale.dffAcute{ax}, scaleEvent{x}(s,ax).dffPeak];
                stillSumm{x}.scale.durAcute{d} = [stillSumm{x}.scale.durAcute{d}, scaleEvent{x}(s,d).dur];
            elseif ismember(s,stillSumm{x}.sPost)
                %stillSumm{x}.scale.dffPost{ax} = [stillSumm{x}.scale.dffPost{ax}, scaleEvent{x}(s,ax).dffPeak];
                stillSumm{x}.scale.durPost{d} = [stillSumm{x}.scale.durPost{d}, scaleEvent{x}(s,d).dur];
            end
        end
    end
    stillSumm{x}.scale.durMed = [cellfun(@median, stillSumm{x}.scale.durPre)', cellfun(@median, stillSumm{x}.scale.durAcute)', cellfun(@median, stillSumm{x}.scale.durPost)']';
    stillSumm{x}.scale.durMed(Ncase < NeventMin) = NaN;
end


%% Show concatenated event rasters
epochRaster = cell(1,Nexpt); allStillRaster = cell(1,Nexpt);
LW = 1.25; transparency = 0.7;
figPath = sprintf('%sStillEpochsExptSummary_2D.pdf', figDir );
close all;  figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'color','w'); 
sp(1) = subplot(3,3,[1:2,4:5,7:8]); sp(2) = subplot(3,3,3);  sp(3) = subplot(3,3,6); sp(4) = subplot(3,3,9);
for x = xPresent
    for s = 1:stillSumm{x}.Nepoch
        epochRaster{x}{s} = false(stillEpoch{x}(s).Nscan, expt(x).Nroi+Ndim);
        for roi = find(stillSumm{x}.fluor.Nevent(s,:) > 0)
            for e = 1:fluorEvent{x}(s,roi).Nevent
                epochRaster{x}{s}(fluorEvent{x}(s,roi).scan{e},roi) = true;
            end
        end
        for d = 1:Ndim
            for e = 1:scaleEvent{x}(s,d).Nevent
                epochRaster{x}{s}(scaleEvent{x}(s,d).scan{e},expt(x).Nroi+d) = true;
            end
        end
          %{
        cla;
        imagesc(epochRaster{x}{s}'); hold on;
        axis tight;
        line([0,stillEpoch{x}(s).Nscan+0.5], (expt(x).Nroi+0.5)*[1,1], 'color','k');
        set(gca,'Ytick',[1,expt(x).Nroi+1], 'YtickLabel',{'ROI','Scale'}, 'TickDir','out', 'TickLength',[0.01,0]);
        impixelinfo;
        title( sprintf('[x,s] = [%i, %i]', x, s) );
        pause;
        %}
    end
    
    allStillRaster{x} = vertcat(epochRaster{x}{:});
    allStillLims = [0,cumsum([stillEpoch{x}.Nscan])]+0.5;
    subplot(sp(1));
    cla;
    imagesc(allStillRaster{x}'); hold on;
    line([0,size(allStillRaster{x},1)+0.5], (expt(x).Nroi+0.5)*[1,1], 'color','k', 'LineWidth',LW);
    for s = 1:numel(allStillLims)
        line( allStillLims(s)*[1,1], [0.5,expt(x).Nroi+Ndim+0.5], 'color','c', 'LineWidth',LW);
        if ismember(s+1, stillSumm{x}.sAcute), line( allStillLims(s)*[1,1], [0.5,expt(x).Nroi+Ndim+0.5], 'color','r', 'LineWidth',LW); end
        if ismember(s+1, stillSumm{x}.sPost), line( allStillLims(s)*[1,1], [0.5,expt(x).Nroi+Ndim+0.5], 'color','g', 'LineWidth',LW); end
    end
    axis tight;
    set(gca,'Ytick',[1,expt(x).Nroi+1], 'YtickLabel',{'ROI','Scale'}, 'TickDir','out', 'TickLength',[0.003,0], 'Xtick',allStillLims(1:end-1), 'XtickLabel',sprintfc('%2.1f',stillSumm{x}.Trel)); %{}
    xtickangle(45);
    impixelinfo;
    xlabel('Peri-CSD Time (min)');
    title( sprintf('x = %i: %s', x, expt(x).name), 'Interpreter','none' );
        
    subplot(sp(2));
    cla;
    tempFo = stillSumm{x}.fluor.Fo./median(stillSumm{x}.fluor.Fo(stillSumm{x}.sPre,:),1);
    medFo = median(tempFo, 2);
    tempFroi = vertcat(stillEpoch{x}.Froi);
    tempFroi = tempFroi./median(tempFroi(stillSumm{x}.sPre,:),1);
    medFroi = median(tempFroi, 2);
    tempFnp = vertcat(stillEpoch{x}.Fnp);
    tempFnp = tempFnp./median(tempFnp(stillSumm{x}.sPre,:),1);
    medFnp = median(tempFnp, 2);
    plot( stillSumm{x}.Trel, medFo, 'color','k'); hold on;  % , 'color',exptColor(k,:)
    h(1) = plot( stillSumm{x}.Trel, medFo, 'k.');  % , 'color',exptColor(k,:)
    plot( stillSumm{x}.Trel, medFroi, 'color','b'); hold on;  % , 'color',exptColor(k,:)
    h(2) = plot( stillSumm{x}.Trel, medFroi, 'b.');  % , 'color',exptColor(k,:)
    plot( stillSumm{x}.Trel, medFnp, 'color','r'); hold on;  % , 'color',exptColor(k,:)
    h(3) = plot( stillSumm{x}.Trel, medFnp, 'r.');  % , 'color',exptColor(k,:)
    tempPos = get(gca,'Position');
    legend(h, {'Fo','Froi','Fneuropil'}, 'AutoUpdate',false, 'Location','NorthEastOutside')
    tempYlim = [0.9, 1.1];%get(gca,'Ylim');
    
    line([0,0], tempYlim, 'color','r', 'lineStyle','--');
    line(acuteLimit*[1,1], tempYlim, 'color','g', 'lineStyle','--');
    line(postLimit*[1,1], tempYlim, 'color','g', 'lineStyle','--'); hold off;
    set(gca,'Xtick',[-120:15:120], 'Position',tempPos);
    ylim(tempYlim); xlim([-Inf,120]);
    xlabel('Peri-CSD Time (min)'); ylabel('Median Normalized Fo');
    
    subplot(sp(3));
    cla;
    plot( stillSumm{x}.Trel, stillSumm{x}.fluor.epochRate, 'color','k'); hold on;  % , 'color',exptColor(k,:)
    plot( stillSumm{x}.Trel, stillSumm{x}.fluor.epochRate, 'k.');  % , 'color',exptColor(k,:)
    tempYlim = get(gca,'Ylim');
    line([0,0], tempYlim, 'color','r', 'lineStyle','--');
    line(acuteLimit*[1,1], tempYlim, 'color','g', 'lineStyle','--');
    line(postLimit*[1,1], tempYlim, 'color','g', 'lineStyle','--'); hold off;
    set(gca,'Xtick',[-120:15:120]);
    ylim(tempYlim); xlim([-Inf,120]); % -120
    xlabel('Peri-CSD Time (min)'); ylabel('Events/min/ROI');
    
    subplot(sp(4));
    cla;
    if ~isempty(stillSumm{x}.sPre)
        plot( stillSumm{x}.scale.CDFx{1}, stillSumm{x}.scale.CDFf{1}, 'color', [0,0,1,transparency], 'LineWidth',LW ); hold on;
    end
    if ~isempty(stillSumm{x}.sAcute)
        plot( stillSumm{x}.scale.CDFx{2}, stillSumm{x}.scale.CDFf{2}, 'color', [1,0,0,transparency], 'LineWidth',LW ); hold on;
    end
    if ~isempty(stillSumm{x}.sPost)
        plot( stillSumm{x}.scale.CDFx{3}, stillSumm{x}.scale.CDFf{3}, 'color', [0,1,0,transparency], 'LineWidth',LW ); hold on;
    end
    axis square;
    xlabel('Scale Magnitude (\mum)'); ylabel('CDF'); title('Pre = blue, acute = red, post = green'); %title( sprintf('x = %i', x) );
    xlim([0.1,10]);
    set(gca,'Xscale','log');

    pause(0.2);
    export_fig( figPath, '-pdf', '-painters','-q101', '-append', gcf );
end

%% Do invidivdual ROIs exhibit sustained increase/decrease in event rate post CSD?
k = 0;
barFaceColor = {'k','r','b'}; % distinguishable_colors(3);
for x = xPresent
    sPostGood = find( stillSumm{x}.Trel > acuteLimit & stillSumm{x}.Trel < postLimit ); % restrict analysis to epochs ~15 - 60 mins post-CSD
    stillSumm{x}.fluor.rateEffect = nan(expt(x).Nroi,5); stillSumm{x}.fluor.baseEffect = nan(expt(x).Nroi,5);
    if numel(sPostGood) > 2 && numel(stillSumm{x}.sPre) > 2
        k = k+1;
        for roi = 1:expt(x).Nroi
            preRateTemp = stillSumm{x}.fluor.eventRate(stillSumm{x}.sPre,roi);
            stillSumm{x}.fluor.rateEffect(roi,1) = mean(preRateTemp, 'omitnan'); % preRateMed
            postRateTemp = stillSumm{x}.fluor.eventRate(sPostGood,roi);
            stillSumm{x}.fluor.rateEffect(roi,2) = mean(postRateTemp, 'omitnan'); % postRateMed
            [stillSumm{x}.fluor.rateEffect(roi,5), stillSumm{x}.fluor.rateEffect(roi,4)] = ttest2(preRateTemp, postRateTemp, 'alpha',0.05);
            %{
            JitterPlot([{preRateTemp}, {postRateTemp}] )
            xlim([0.75,2.25]);
            title( sprintf('[x,roi] = [%i, %i]', x, roi));
            pause;
            %}
            preBaseTemp = stillSumm{x}.fluor.Fo(stillSumm{x}.sPre,roi);
            stillSumm{x}.fluor.baseEffect(roi,1) = mean(preBaseTemp, 'omitnan'); % preRateMed
            postBaseTemp = stillSumm{x}.fluor.Fo(sPostGood,roi);
            stillSumm{x}.fluor.baseEffect(roi,2) = mean(postBaseTemp, 'omitnan'); % postBaseMed
            [stillSumm{x}.fluor.baseEffect(roi,5), stillSumm{x}.fluor.baseEffect(roi,4)] = ttest2(preBaseTemp, postBaseTemp, 'alpha',0.05);
            %{
            cla;
            JitterPlot([{preBaseTemp}, {postBaseTemp}] );
            xlim([0.75,2.25]);
            title( sprintf('[x,roi] = [%i, %i]', x, roi));
            pause;
            %}
        end
        % Which units' event rates changed?
        stillSumm{x}.fluor.rateEffect(:,3) = stillSumm{x}.fluor.rateEffect(:,2) - stillSumm{x}.fluor.rateEffect(:,1);
        stillSumm{x}.fluor.rateUp = find(stillSumm{x}.fluor.rateEffect(:,5) == 1 & stillSumm{x}.fluor.rateEffect(:,3) > 0)';
        stillSumm{x}.fluor.NrateUp = numel(stillSumm{x}.fluor.rateUp);
        stillSumm{x}.fluor.rateDown = find(stillSumm{x}.fluor.rateEffect(:,5) == 1 & stillSumm{x}.fluor.rateEffect(:,3) < 0)';
        stillSumm{x}.fluor.NrateDown = numel(stillSumm{x}.fluor.rateDown);
        stillSumm{x}.fluor.rateNon = setdiff(1:expt(x).Nroi, [stillSumm{x}.fluor.rateUp, stillSumm{x}.fluor.rateDown]); %find(~stillSumm{x}.fluor.rateEffect(:,5))';
        stillSumm{x}.fluor.NrateNon = numel(stillSumm{x}.fluor.rateNon);
        % Which units' baseline fluorescence changed?
        stillSumm{x}.fluor.baseEffect(:,3) = (stillSumm{x}.fluor.baseEffect(:,2) - stillSumm{x}.fluor.baseEffect(:,1))./stillSumm{x}.fluor.baseEffect(:,1);
        stillSumm{x}.fluor.baseUp = find(stillSumm{x}.fluor.baseEffect(:,5) == 1 & stillSumm{x}.fluor.baseEffect(:,3) > 0)';
        stillSumm{x}.fluor.NbaseUp = numel(stillSumm{x}.fluor.baseUp);
        stillSumm{x}.fluor.baseDown = find(stillSumm{x}.fluor.baseEffect(:,5) == 1 & stillSumm{x}.fluor.baseEffect(:,3) < 0)';
        stillSumm{x}.fluor.NbaseDown = numel(stillSumm{x}.fluor.baseDown);
        stillSumm{x}.fluor.baseNon = setdiff(1:expt(x).Nroi, [stillSumm{x}.fluor.baseUp, stillSumm{x}.fluor.baseDown]); %find(~stillSumm{x}.fluor.baseEffect(:,5))';
        stillSumm{x}.fluor.NbaseNon = numel(stillSumm{x}.fluor.baseNon);

        subplot(2,3,1); cla;
        if stillSumm{x}.fluor.NrateNon > 0
            plot( [1,2], stillSumm{x}.fluor.baseEffect(stillSumm{x}.fluor.baseNon ,1:2), 'color',[0,0,0,transparency] ); hold on;
            plot( [1,2], stillSumm{x}.fluor.baseEffect(stillSumm{x}.fluor.baseNon,1:2), '.', 'color',[0,0,0,transparency] );
        end
        if stillSumm{x}.fluor.NbaseUp > 0
            plot( [1,2], stillSumm{x}.fluor.baseEffect(stillSumm{x}.fluor.baseUp ,1:2), 'color',[0,0,1,transparency] ); 
            plot( [1,2], stillSumm{x}.fluor.baseEffect(stillSumm{x}.fluor.baseUp,1:2), '.', 'color',[0,0,1,transparency] );
        end
        if stillSumm{x}.fluor.NbaseUp > 0
            plot( [1,2], stillSumm{x}.fluor.baseEffect(stillSumm{x}.fluor.baseDown ,1:2), 'color',[1,0,0,transparency] ); 
            plot( [1,2], stillSumm{x}.fluor.baseEffect(stillSumm{x}.fluor.baseDown,1:2), '.', 'color',[1,0,0,transparency] );
        end
        xlim([0.8,2.2]); set(gca,'Xtick',[1,2], 'XtickLabel',{'Pre-CSD','Post-CSD'});
        ylabel('F_0');
        title( sprintf('x = %i', x) );
        axis square;
        
        subplot(2,3,4);
        h = bar( k, [stillSumm{x}.fluor.NbaseNon, stillSumm{x}.fluor.NbaseDown, stillSumm{x}.fluor.NbaseUp]/expt(x).Nroi, 'stacked' ); hold on;
        for s = 1:numel(h), h(s).FaceColor = barFaceColor{s}; end % barFaceColor(s,:)
        legend(h, {'No change','Down','Up'}, 'AutoUpdate',false, 'Location','northeastoutside');
        %h.FaceColor = barFaceColor;
        axis square; 
        ylabel('Fraction of ROI'); xlabel('Experiment'); title('Base Fluor Effect');
        
        subplot(2,3,2); cla;
        if stillSumm{x}.fluor.NrateNon > 0
            plot( [1,2], stillSumm{x}.fluor.rateEffect(stillSumm{x}.fluor.rateNon ,1:2), 'color',[0,0,0,transparency] ); hold on;
            plot( [1,2], stillSumm{x}.fluor.rateEffect(stillSumm{x}.fluor.rateNon,1:2), '.', 'color',[0,0,0,transparency] );
        end
        if stillSumm{x}.fluor.NrateUp > 0
            plot( [1,2], stillSumm{x}.fluor.rateEffect(stillSumm{x}.fluor.rateUp ,1:2), 'color',[0,0,1,transparency] ); 
            plot( [1,2], stillSumm{x}.fluor.rateEffect(stillSumm{x}.fluor.rateUp,1:2), '.', 'color',[0,0,1,transparency] );
        end
        if stillSumm{x}.fluor.NrateUp > 0
            plot( [1,2], stillSumm{x}.fluor.rateEffect(stillSumm{x}.fluor.rateDown ,1:2), 'color',[1,0,0,transparency] ); 
            plot( [1,2], stillSumm{x}.fluor.rateEffect(stillSumm{x}.fluor.rateDown,1:2), '.', 'color',[1,0,0,transparency] );
        end
        xlim([0.8,2.2]); set(gca,'Xtick',[1,2], 'XtickLabel',{'Pre-CSD','Post-CSD'});
        ylabel('Events per minute');
        title( sprintf('x = %i', x) );
        axis square;
        
        subplot(2,3,5);
        h = bar( k, [stillSumm{x}.fluor.NrateNon, stillSumm{x}.fluor.NrateDown, stillSumm{x}.fluor.NrateUp]/expt(x).Nroi, 'stacked' ); hold on;
        for s = 1:numel(h), h(s).FaceColor = barFaceColor{s}; end % barFaceColor(s,:)
        %legend(h, {'No change','Down','Up'}, 'AutoUpdate',false, 'Location','northwestoutside');
        %h.FaceColor = barFaceColor;
        axis square; 
        ylabel('Fraction of ROI'); xlabel('Experiment'); title('Event Rate Effect');
        
        subplot(2,3,[3,6]); cla;
        plot( stillSumm{x}.fluor.baseEffect(:,3), stillSumm{x}.fluor.rateEffect(:,3), '.' ); hold on;
        line(0.5*[-1,1], [0,0], 'color','k');
        line([0,0], [-2,2], 'color','k');
        axis square;
        xlabel('Relative Change in Baseline Fluor'); ylabel('Change in Event Rate');
        
        %pause;
    end
end

%% Summarize proportions of ROI that went up/down/unchanged in rate or Fo
close all;  clearvars h g sp;
figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'color','w'); 
legendCell = cell(0);
sp(2)= subplot(1,2,2); sp(1)= subplot(1,2,1); 
k = 0;
for x = xPresent
    sPostGood = find( stillSumm{x}.Trel > acuteLimit & stillSumm{x}.Trel < postLimit ); % restrict analysis to epochs ~15 - 60 mins post-CSD
    if numel(sPostGood) > 2 && numel(stillSumm{x}.sPre) > 2
        k = k+1;
        legendCell{k} = expt(x).name;
        subplot(sp(1));
        h = bar( k, [stillSumm{x}.fluor.NbaseNon, stillSumm{x}.fluor.NbaseDown, stillSumm{x}.fluor.NbaseUp]/expt(x).Nroi, 'stacked' ); hold on;
        for s = 1:numel(h), h(s).FaceColor = barFaceColor{s}; end % barFaceColor(s,:)
        axis square; 
        ylabel('Fraction of ROI'); xlabel('Experiment'); title('F_0 Sustained Effect');

        subplot(sp(2));
        g = bar( k, [stillSumm{x}.fluor.NrateNon, stillSumm{x}.fluor.NrateDown, stillSumm{x}.fluor.NrateUp]/expt(x).Nroi, 'stacked' ); hold on;
        for s = 1:numel(g), g(s).FaceColor = barFaceColor{s}; end % barFaceColor(s,:)
        axis square; 
        ylabel('Fraction of ROI'); xlabel('Experiment'); title('Event Rate Sustained Effect');
    end
end
subplot(sp(1)); set(gca,'Position',tempPos, 'Xtick',1:k, 'XtickLabel',legendCell); xtickangle(45);
tempPos = get(gca,'Position');
legend(h, {'No change','Down','Up'}, 'AutoUpdate',false, 'Location','northeastoutside');
set(gca,'Position',tempPos)
subplot(sp(2));
set(gca,'Xtick',1:k, 'XtickLabel',legendCell); xtickangle(45);

%% Compare Fo effect on fibers vs other ROI

%for x = 


%% Does base (pre-CSD) fluor bleach?
figure;
for x = 37 %xPresent([expt(xPresent).Nplane] > 1) % xPresent
    Ttemp = vertcat( Tscan{x}{:} )/60;
    tempFluorF = [fluor{x}.F];
    tempFluorROI = [fluor{x}.Froi];
    tempFluorNP = [fluor{x}.Fnp];
    
    FvolTemp = vertcat(tempFluorF.vol);
    FplaneTemp = vertcat(tempFluorF.plane); 
    FplaneBase = zeros(size(FplaneTemp));
    for z = 1:size(FplaneTemp,2)
        FplaneBase(:,z) = MovingPercentile(FplaneTemp(:,z), 10, round(expt(x).scanRate*300), 'center');
    end
    FnpTemp = vertcat(tempFluorNP.all); % ROI
    FroiTemp = vertcat(tempFluorROI.all);  % ROI
    %baseFvolTemp = MovingPercentile(allFvolTemp, 5, round(expt(x).scanRate*300), 'center');
    %baseFnpTemp = MovingPercentile(allFnpTemp, 5, round(expt(x).scanRate*300), 'center');
    %h(1) = plot( Tbase, FvolTemp, 'k' ); hold on;
    if expt(x).Nplane > 1
        for z = 1:size(FplaneTemp,2)
            plot(Ttemp, FplaneTemp(:,z) ); hold on;
            plot(Ttemp,  FplaneBase(:,z), 'k');
            %pause;
        end
        xlim([-Inf,Inf]);
        title( sprintf('x = %i: %s', x, expt(x).name), 'Interpreter','none');
        ylabel('Mean Planar Fluorescence'); xlabel('Time (min)');
    end
    pause;
    cla;
    %plot( Tbase, FplaneTemp ); 
    %plot( Tbase, FroiTemp );   
    %plot( Tbase, FnpTemp );
    %legend(h, {'Volume','ROI','Neuropil'}, 'AutoUpdate',false);
    %hold off; 
    
    %plot( Ttemp, allFnpTemp ); hold on;
    %plot( Ttemp, baseFnpTemp ); hold off;
end

%% Compare fluor event rates before/after CSD - individual experiments
k = 0;
close all; 
figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'color','w'); 
opt = {[0.05,0.05], [0.08,0.08], [0.05,0.05]};  % {[vert, horz], [bottom, top], [left, right] }
clearvars sp;
for p = flip(1:3), sp(p) = subtightplot(1,3,p,opt{:});   end
exptColor = distinguishable_colors(Npresent);
exptNames = cell(0);
for x = intersect(xPresent, xCSD)
    k = k+1;
    
    subplot(sp(1));
    plot( stillSumm{x}.Trel, stillSumm{x}.fluor.epochRate, 'color',exptColor(k,:)); hold on;
    plot( stillSumm{x}.Trel, stillSumm{x}.fluor.epochRate, '.', 'color',exptColor(k,:)); 
    axis square;
    xlabel('Peri-CSD Time (min)'); ylabel('Events/min/ROI');
    
    subplot(sp(2)); %cla;
    imagesc( stillSumm{x}.fluor.eventRate' ); hold on;
    if ~isempty(stillSumm{x}.sPre)
        line( stillSumm{x}.sPre(end)+0.5*[1,1], [0,expt(x).Nroi] ); hold off;
    end
    colorbar;
    xlabel('Epoch'); ylabel('ROI');
    title( sprintf('%s', expt(x).name), 'Interpreter','none');
    
    subplot(sp(3));
    if ~isempty(stillSumm{x}.sPost)
        h(k) = plot(stillSumm{x}.fluor.PostRateCDF(:,1), stillSumm{x}.fluor.PostRateCDF(:,2), '--', 'color',exptColor(k,:)); hold on;
    end
    if ~isempty(stillSumm{x}.sPre)
        h(k) = plot(stillSumm{x}.fluor.PreRateCDF(:,1), stillSumm{x}.fluor.PreRateCDF(:,2), 'color',exptColor(k,:)); hold on;
    end
    exptNames{k} = expt(x).name;
    xlabel('ROI Event Rate (events/minute stillness)'); ylabel('Fraction of ROI');
    axis square;
    pause;
end
subplot(sp(3));
legend(h, exptNames, 'Location','southeast','Interpreter','none')


%% Compare scale event rates before/after CSD - individual experiments
k = 0;
close all; 
figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'color','w'); 
opt = {[0.05,0.05], [0.08,0.08], [0.05,0.05]};  % {[vert, horz], [bottom, top], [left, right] }
clearvars sp;
for p = flip(1:3), sp(p) = subtightplot(1,3,p,opt{:});   end
exptColor = distinguishable_colors(Npresent);
exptNames = cell(0);
for x = intersect(xPresent, xCSD)
    k = k+1;
    
    subplot(sp(1));
    plot( stillSumm{x}.Trel, stillSumm{x}.scale.epochRate, 'color',exptColor(k,:)); hold on;
    plot( stillSumm{x}.Trel, stillSumm{x}.scale.epochRate, '.', 'color',exptColor(k,:)); 
    axis square;
    xlabel('Peri-CSD Time (min)'); ylabel('Events/min/ROI');
    
    subplot(sp(2)); %cla;
    imagesc( stillSumm{x}.scale.eventRate' ); hold on;
    if ~isempty(stillSumm{x}.sPre)
        line( stillSumm{x}.sPre(end)+0.5*[1,1], [0,expt(x).Nroi] ); hold off;
    end
    colorbar;
    xlabel('Epoch'); ylabel('ROI');
    title( sprintf('%s', expt(x).name), 'Interpreter','none');
    
    subplot(sp(3));
    if ~isempty(stillSumm{x}.sPost)
        h(k) = plot(stillSumm{x}.scale.PostRateCDF(:,1), stillSumm{x}.scale.PostRateCDF(:,2), '--', 'color',exptColor(k,:)); hold on;
    end
    if ~isempty(stillSumm{x}.sPre)
        h(k) = plot(stillSumm{x}.scale.PreRateCDF(:,1), stillSumm{x}.scale.PreRateCDF(:,2), 'color',exptColor(k,:)); hold on;
    end
    exptNames{k} = expt(x).name;
    xlabel('ROI Event Rate (events/minute stillness)'); ylabel('Fraction of ROI');
    axis square;
    pause;
end
subplot(sp(3));
legend(h, exptNames, 'Location','southeast','Interpreter','none')

%% Compare fluor event rates before/after CSD - individual experiments
k = 0;
close all; 
figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'color','w'); 
opt = {[0.05,0.05], [0.08,0.08], [0.05,0.05]};  % {[vert, horz], [bottom, top], [left, right] }
clearvars sp;
for p = flip(1:3), sp(p) = subtightplot(1,3,p,opt{:});   end
exptColor = distinguishable_colors(Npresent);
exptNames = cell(0);
for x = intersect(xPresent, xCSD)
    k = k+1;
    
    subplot(sp(1));
    plot( stillSumm{x}.Trel, stillSumm{x}.fluor.epochRate, 'color',exptColor(k,:)); hold on;
    plot( stillSumm{x}.Trel, stillSumm{x}.fluor.epochRate, '.', 'color',exptColor(k,:)); 
    axis square;
    xlabel('Peri-CSD Time (min)'); ylabel('Events/min/ROI');
    
    subplot(sp(2)); %cla;
    imagesc( stillSumm{x}.scale.eventRate' ); hold on;
    if ~isempty(stillSumm{x}.sPre)
        line( stillSumm{x}.sPre(end)+0.5*[1,1], [0,expt(x).Nroi] ); hold off;
    end
    colorbar;
    xlabel('Epoch'); ylabel('ROI');
    title( sprintf('%s', expt(x).name), 'Interpreter','none');
    
    subplot(sp(3));
    if ~isempty(stillSumm{x}.sPost)
        h(k) = plot(stillSumm{x}.scale.PostRateCDF(:,1), stillSumm{x}.scale.PostRateCDF(:,2), '--', 'color',exptColor(k,:)); hold on;
    end
    if ~isempty(stillSumm{x}.sPre)
        h(k) = plot(stillSumm{x}.scale.PreRateCDF(:,1), stillSumm{x}.scale.PreRateCDF(:,2), 'color',exptColor(k,:)); hold on;
    end
    exptNames{k} = expt(x).name;
    xlabel('Scale Event Rate (events/minute stillness)'); ylabel('Fraction of ROI');
    axis square;
    pause;
end
subplot(sp(3));
legend(h, exptNames, 'Location','southeast','Interpreter','none')

%% Compare events rates before/after CSD - all experiments
k = 0;
close all; 
figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'color','w'); 
opt = {[0.05,0.08], [0.08,0.08], [0.08,0.08]};  % {[vert, horz], [bottom, top], [left, right] }
clearvars sp h;
for p = flip(1:2), sp(p) = subtightplot(1,2,p,opt{:});   end
exptColor = distinguishable_colors(Npresent);
exptNames = cell(0);
for x = intersect(xPresent, xCSD)
    k = k+1;
    
    subplot(sp(1));
    plot( stillSumm{x}.Trel, stillSumm{x}.fluor.epochRate, 'color',exptColor(k,:)); hold on;
    plot( stillSumm{x}.Trel, stillSumm{x}.fluor.epochRate, '.', 'color',exptColor(k,:)); 
    line( [0,0], [0,2.5], 'lineStyle','--', 'color','k');
    line( acuteLimit*[1,1], [0,2.5], 'lineStyle','--', 'color','r');
    xlim([-Inf,Inf]);
    axis square;
    xlabel('Peri-CSD Time (min)'); ylabel('Events/min/ROI');
    
    subplot(sp(2));
    %{
    if ~isempty(stillSumm{x}.sPost)
        h(k) = plot(stillSumm{x}.fluor.PostRateCDF(:,1), stillSumm{x}.fluor.PostRateCDF(:,2), '--', 'color',exptColor(k,:)); hold on;
    end
    %}
    %if ~isempty(stillSumm{x}.sPre)
    h(k) = plot(stillSumm{x}.fluor.PreRateCDF(:,1), stillSumm{x}.fluor.PreRateCDF(:,2), 'color',exptColor(k,:)); hold on;
    %end
    exptNames{k} = expt(x).name;
    xlabel('ROI Event Rate (events/minute stillness)'); ylabel('Fraction of ROI');
    axis square;
    %pause;
end
subplot(sp(2));
legend(h, exptNames, 'Location','southeast','Interpreter','none')

%% Compare properties of events pre/post-CSD, pooling all ROI

tempStill = [stillSumm{xPresent}];
poolRate = [horzcat(tempStill.PreRate); horzcat(tempStill.AcuteRate); horzcat(tempStill.PostRate)];
poolDur = horzcat(tempStill.durMed);
poolDFF = horzcat(tempStill.dffMed);
close all; 
figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'color','w'); 
opt = {[0.05,0.08], [0.08,0.08], [0.08,0.08]};  % {[vert, horz], [bottom, top], [left, right] }
clearvars sp h;
subplot(1,3,1);
JitterPlot( poolRate', 0.5, 'dependent',true, 'Box',true, 'ErrorBar',0, 'monochrome',0.4 ); %plot( poolDur ) xlim([0.8,2.2]); %ylim([0,4]); 
axis square;
set(gca,'Xtick',1:3, 'XtickLabel',{'Pre','Acute','Post'}, 'TickDir','out', 'TickLength',[0.001,0], 'box','off'); % , 'Yscale','log'
ylim([0,Inf]);
ylabel('Event Rate (events/minute)');
%pRate = signrank(poolRate(1,:), poolRate(3,:));
%title( sprintf('p = %1.5f', pDur ) );

subplot(1,3,2);
JitterPlot( poolDur', 0.5, 'dependent',true, 'Yscale','log', 'Box',true, 'ErrorBar',0, 'monochrome',0.4 ); %plot( poolDur ) xlim([0.8,2.2]); %ylim([0,4]);
axis square;
set(gca,'Xtick',1:3, 'XtickLabel',{'Pre','Acute','Post'}, 'TickDir','out', 'TickLength',[0.001,0], 'box','off');
ylabel(sprintf('Median Event Duration (s, minimum %i events per ROI per condition)', NeventMin) );
%pDur = signrank(poolDur(1,:), poolDur(2,:));
%title( sprintf('p = %1.5f', pDur ) );

subplot(1,3,3);
JitterPlot( poolDFF', 0.5, 'dependent',true, 'Yscale','log', 'Box',true, 'ErrorBar',0, 'monochrome',0.4 ); %plot( poolDFF ); xlim([0.8,2.2]); ylim([0,4]);
axis square;
set(gca,'Xtick',1:3, 'XtickLabel',{'Pre','Acute','Post'}, 'TickDir','out', 'TickLength',[0.001,0], 'box','off');
ylabel(sprintf('Median Event dF/F_0 (minimum %i events per ROI per condition)', NeventMin) ); 
ylim([0,3]);
%pDFF = signrank(poolDFF(1,:), poolDFF(2,:));
%title( sprintf('p = %1.5f', pDFF ) );

%% Write tifs concatenating all the events for each ROI
for x = 11 %xPresent
    Tcat = vertcat(Tscan{x}{:});
    roiEventMovies = cell(1,expt(x).Nroi);
    tic;
    for roi = 30 %find(totEvent{x}) %1:expt(x).Nroi
        ROI{x}(roi).projPath = sprintf('%s%s_ROI_%02i.tif', [expt(x).dir, 'ROI\'], expt(x).name, roi );
        fprintf('\nLoading %s... ', ROI{x}(roi).projPath);
        roiMovie = loadtiff(ROI{x}(roi).projPath); toc
        blankFrame = prctile(roiMovie(:),10)*ones(size(roiMovie,1), size(roiMovie,2), 'uint16'); % median(roiMovie(:))
        roiEventMovies{roi} = cell(1, stillSumm{x}.fluor.TotEvent(roi));
        k = 0;
        for s = find(stillSumm{x}.fluor.Nevent(:,roi))' %1:stillEventSumm{x}.Nepoch
            for e = 1:stillSumm{x}.fluor.Nevent(s,roi)
                eventCatScans = find(ismember(Tcat, fluorEvent{x}(s,roi).T{e}'))';
                k = k+1;
                roiEventMovies{roi}{k} = cat(3, roiMovie(:,:,eventCatScans), blankFrame);
            end
        end

        eventTifPath = sprintf('%sROI\\%s_ROI_%02i_spontEvents.tif', expt(x).dir, expt(x).name, roi );
        saveastiff(cat(3,roiEventMovies{roi}{:}) , eventTifPath);
        toc
    end   
end

%% Write movies of each epoch of stillness
GrayOpt = struct('overwrite',true, 'message',false, 'append',false, 'big',true, 'color',false );
RGBOpt = struct('overwrite',true, 'message',false, 'append',false, 'big',true, 'color',true );

movParam.binT = 10;
epochStack = cell(1,Nexpt); stampedStack = cell(1,Nexpt); Tepoch = cell(1,Nexpt); 
for x = 28
    movParam.dir = sprintf('%sZtifs\\', expt(x).dir); %strcat('D:\MATLAB\LevyLab\Figures\Fibers\', expt(x).name, '\' );  mkdir(movParam.dir);
    movParam.edges = segParams{x}.edges;
    movParam.scalebar = MakeScaleBar( round(expt(x).umPerPixel*[50,50]), {[0,expt(x).Ncol-segParams{x}.edges(1)-segParams{x}.edges(2)]+0.5, [0,expt(x).Nrow-segParams{x}.edges(3)-segParams{x}.edges(4)]+0.5},...
        [0.1,0.95], [0,0], 'label',false, 'color','w', 'show',false );
    if expt(x).Nplane > 1
        movParam.sourceSbx = 'sbx_interp';  
    else
        movParam.sourceSbx = 'sbx_affine'; %'sbxz'; %
    end
    sourceSbxPath = [expt(x).dir, expt(x).name, '.', movParam.sourceSbx]; %[sbxInfo.dir, catInfo(x).exptName, '.', movParam.sourceSbx]; % 
    Tcat = vertcat(Tscan{x}{:});
    epochStack{x} = cell(1,stillSumm{x}.Nepoch); stampedStack{x} = cell(1,stillSumm{x}.Nepoch); 
    for s = 1:stillSumm{x}.Nepoch 
        fprintf('\n[x, epoch] = [%i,  %i]', x, s);
        boutFileName = sprintf('%s%s_epoch%02i', movParam.dir, expt(x).name, s);
        stackPath = strcat(boutFileName,'_stack.tif');
        
        % Determine epoch's indices in concatenated movies frame
        [~,firstScan] = min( abs(Tcat - stillEpoch{x}(s).Tstart ) );  
        lastScan = firstScan + stillEpoch{x}(s).Nscan - 1;
        %{
        if expt(x).Nplane == 1
            epochStack{x}{s} = WriteSbxPlaneTif(sourceSbxPath, catInfo(x), movParam.zProj, 'firstScan',firstScan, 'Nscan',stillEpoch{x}(s).Nscan, 'binT',movParam.binT, 'verbose',true);     
            epochStack{x}{s} = uint16(epochStack{x}{s}(movParam.edges(3)+1:end-movParam.edges(4), movParam.edges(1)+1:end-movParam.edges(2),:,:));
        else
            error('Only designed for single plane for now');
        end
        % Determine the relative timing of each (possibly binned) frame/scan
        %}        
        Tepoch{x}{s} = BinDownMean( Tcat(firstScan:lastScan) - csdBout{x}(expt(x).csd).Tstart , movParam.binT ); % - Tcat(firstScan)
        for z = 1:size(epochStack{x}{s}, 3)
            stampedStack{x}{s}(:,:,z) = rgb2gray(insertText( epochStack{x}{s}(:,:,z), [550,20], sprintf('%2.2f s', Tepoch{x}{s}(z)), 'BoxOpacity',0, 'TextColor','w', 'Font','Arial', 'FontSize',16) );
        end
        
        % Write the movie to tiff
        fprintf('\nWriting %s...', stackPath); tic;
        saveastiff(stampedStack{x}{s}, stackPath, GrayOpt); toc % RGBOpt
    end
end

%% Write movies of each ROI's brightest event
for x = 28
    tic
    roiEventMov = cell(1,expt(x).Nroi);
    sMax = nan(1,expt(x).Nroi); eMax = nan(1,expt(x).Nroi);
    Tcat = vertcat(Tscan{x}{:});
    for roi = 1:18 %expt(x).Nroi
        % Find the ROI's brightest event
        tempMaxDFF = 0;
        for s = find(stillSumm{x}.fluor.Nevent(:,roi))'
            [dffCurrMax, eCurrMax] = max(fluorEvent{x}(s,roi).dffPeak);
            if dffCurrMax > tempMaxDFF
                tempMaxDFF = dffCurrMax;
                sMax(roi) = s;
                eMax(roi) = eCurrMax;
            end
        end
        
        % Write the frames to a tiff and AVI
        if expt(x).Nplane == 1
            [~,~,tempEventCatScans] = intersect( fluorEvent{x}(sMax(roi),roi).T{eMax(roi)}, Tcat );
            roiEventMov{roi} = WriteSbxPlaneTif(expt(x).sbx, catInfo(x), 1, 'firstScan',tempEventCatScans(1), 'Nscan',numel(tempEventCatScans), 'binT',10, 'verbose',true);   % 
        end
    end
    toc
    
    % Concatenate movies and draw each frame, indicating ROI
    roiCatMov = cat(3, roiEventMov{:});
    roiCatMov = roiCatMov(segParams{x}.edges(3)+1:end-segParams{x}.edges(4), segParams{x}.edges(1)+1:end-segParams{x}.edges(2),:);
    NcatFrame = size(roiCatMov, 3);
    displayLims = [prctile(roiCatMov(:), movParam.displayPct(1)), prctile(roiCatMov(:), movParam.displayPct(2))];
    zCurr = 1;
    close all; clearvars storeFrames
    tempFrameFile = strcat(expt(x).dir, 'tempFrame.jpg');
    storeFrames = repmat(struct('cdata',[], 'colormap',[]), 1, NcatFrame);
    aviFig = figure('color','k', 'Units','normalized', 'OuterPosition', [0 0 1 1]); % 'WindowState','maximized',
    w = waitbar(0, 'Writing concatenated movie');
    for roi = 1:expt(x).Nroi
        for z = zCurr:zCurr+size(roiEventMov{roi},3)-1 %1:size(roiEventMov{roi},3)
            imshow(roiCatMov(:,:,z), displayLims, 'border','tight'); hold on;
            plot( ROI{x}(roi).footprintEdge(:,2)-segParams{x}.edges(1)+1, ROI{x}(roi).footprintEdge(:,1)-segParams{x}.edges(3)+1, '.', 'color','w', 'MarkerSize',1 );
            text(600, 20, sprintf('ROI %02i',roi), 'HorizontalAlignment','center', 'color','w', 'FontSize',20 );
            %impixelinfo
            
            exportgraphics(gca, tempFrameFile, 'Resolution',100); % save image with set resolution 
            storeFrames(z) = im2frame(imread(tempFrameFile));  % getframe(aviFig);  % convert image to frame
            waitbar(z/NcatFrame, w)
            
        end
        zCurr = z+1;
        %pause;
    end
    delete(w); delete(tempFrameFile);
    toc
    
    % Write AVI
    aviPath = sprintf('%s%s_ROIevents.avi', expt(x).dir, expt(x).name);
    fprintf('\nWriting %s...  ', aviPath); tic;
    writerObj = VideoWriter(aviPath);
    writerObj.FrameRate = movParam.aviRate; % set the seconds per image
    open(writerObj);
    for z = 1:numel(storeFrames)
        frame = storeFrames(z);    
        writeVideo(writerObj, frame);
    end
    close(writerObj);
    toc 
    
end

%% Write movies of each epoch of stillness for each ROI
GrayOpt = struct('overwrite',true, 'message',false, 'append',false, 'big',true, 'color',false );
RGBOpt = struct('overwrite',true, 'message',false, 'append',false, 'big',true, 'color',true );

movParam.binT = 10;
epochStack = cell(1,Nexpt); Tepoch = cell(1,Nexpt); roiStack = cell(1,Nexpt);
for x = 28
    movParam.dir = sprintf('%sZtifs\\', expt(x).dir); %strcat('D:\MATLAB\LevyLab\Figures\Fibers\', expt(x).name, '\' );  mkdir(movParam.dir);
    movParam.edges = segParams{x}.edges;
    movParam.scalebar = MakeScaleBar( round(expt(x).umPerPixel*[50,50]), {[0,expt(x).Ncol-segParams{x}.edges(1)-segParams{x}.edges(2)]+0.5, [0,expt(x).Nrow-segParams{x}.edges(3)-segParams{x}.edges(4)]+0.5},...
        [0.1,0.95], [0,0], 'label',false, 'color','w', 'show',false );
    if expt(x).Nplane > 1
        movParam.sourceSbx = 'sbx_interp';  
    else
        movParam.sourceSbx = 'sbx_affine'; %'sbxz'; %
    end
    sourceSbxPath = [expt(x).dir, expt(x).name, '.', movParam.sourceSbx]; %[sbxInfo.dir, catInfo(x).exptName, '.', movParam.sourceSbx]; % 
    Tcat = vertcat(Tscan{x}{:});
    epochStack{x} = cell(stillSumm{x}.Nepoch, 1);  
    for s = find(cellfun(@isempty, epochStack{x}))' %1:stillSumm{x}.Nepoch 
        fprintf('\n[x, epoch / ] = [%i,  %i / %i]', x, s, stillSumm{x}.Nepoch );
        boutFileName = sprintf('%s%s_epoch%02i', movParam.dir, expt(x).name, s);
        stackPath = strcat(boutFileName,'_stack.tif');
        % Determine epoch's indices in concatenated movies frame
        [~,firstScan] = min( abs(Tcat - stillEpoch{x}(s).Tstart ) );  
        lastScan = firstScan + stillEpoch{x}(s).Nscan - 1;
        Tepoch{x}{s} = BinDownMean( Tcat(firstScan:lastScan) - csdBout{x}(expt(x).csd).Tstart, movParam.binT ); % - Tcat(firstScan)
        if expt(x).Nplane == 1
            movParam.zProj = 1;
            epochStack{x}{s} = uint16(WriteSbxPlaneTif(sourceSbxPath, catInfo(x), movParam.zProj, 'firstScan',firstScan, 'Nscan',stillEpoch{x}(s).Nscan, 'binT',movParam.binT, 'verbose',true));     
            %epochStack{x}{s} = epochStack{x}{s}(movParam.edges(3)+1:end-movParam.edges(4), movParam.edges(1)+1:end-movParam.edges(2),:,:));
        else
            error('Only designed for single plane for now');
        end
        epochStack{x}{s}(:,:,end+1) = zeros( expt(x).Nrow, expt(x).Ncol, 'uint16' );
    end
    tic
    allEpochStack = cat(3, epochStack{x}{:});
    toc
    saveastiff( allEpochStack(movParam.edges(3)+1:end-movParam.edges(4), movParam.edges(1)+1:end-movParam.edges(2),:), [expt(x).dir, 'allEpoch.tif'], GrayOpt); 
    toc
    
    tic
    median(epochStack{x}{1}, 3);
    toc
    
    epochMedProj = cellfun(@median, epochStack{x}, repmat({3}, stillSumm{x}.Nepoch, 1), 'UniformOutput',false);
    epochMedProj = cat(3, epochMedProj{:});
    % {
    roiStack{x} = cell(stillSumm{x}.Nepoch, expt(x).Nroi); %
    tic
    for roi = 1:expt(x).Nroi
        %tempROIpath = sprintf('%sROI\\%s_ROI%03i_allStillEpoch.tif', expt(x).dir, expt(x).name, roi);%[expt(x).dir,'ROI//, 'allEpoch.tif']
        %tempStack = allEpochStack(ROI{x}(roi).crop(3):ROI{x}(roi).crop(4), ROI{x}(roi).crop(1):ROI{x}(roi).crop(2),:);
        %saveastiff( tempStack, tempROIpath, GrayOpt);
        tempProjPath = sprintf('%sROI\\%s_ROI%03i_allStillEpochProj.tif', expt(x).dir, expt(x).name, roi);%[expt(x).dir,'ROI//, 'allEpoch.tif']
        tempStack = epochMedProj(ROI{x}(roi).crop(3):ROI{x}(roi).crop(4), ROI{x}(roi).crop(1):ROI{x}(roi).crop(2),:);
        saveastiff( tempStack, tempProjPath, GrayOpt);
        %roiStack{x}{s,roi} = epochStack{x}{s}(ROI{x}(roi).crop(3):ROI{x}(roi).crop(4), ROI{x}(roi).crop(1):ROI{x}(roi).crop(2),:);
        toc
    end
    %}
end

%% Write z-proj with (color-coded) CSD-affected ROIs overlaid
% Load zproj tif
movParam.displayPct = [1,99.5];
movParam.aviRate = 10; % frames per second
GrayOpt = struct('overwrite',true, 'message',false, 'append',false, 'big',true, 'color',false );
darkParam = 0.85; %1; % 0.75;
darkBlue = [0,0,darkParam]; 
darkRed = [darkParam,0,0];
for x = 28 %xPresent
    movParam.binT = round(expt(x).scanRate*30); %10;
    Tcat = vertcat( Tscan{x}{:} );
    if ~isnan(expt(x).csd), Tcat = Tcat - csdBout{x}(expt(x).csd).Tstart; end
    Tbin = BinDownMean( Tcat, movParam.binT );
    tempFrameFile = strcat(expt(x).dir, 'tempFrame.jpg');
    tic
    % Generate a projection stack
    if expt(x).Nplane == 1   
        ZprojStack = WriteSbxPlaneTif(expt(x).sbx, catInfo(x), 1, 'firstScan',1, 'Nscan',expt(x).totScan, 'binT',movParam.binT, 'verbose',true, ...
            'dir',expt(x).dir, 'name',sprintf('%s_Bin%i',expt(x).name, movParam.binT), 'edges',segParams{x}.edges);   % 
        NprojFrame = size(ZprojStack,3);
    else
        tempTifName = sprintf('%s%s_bin%i.tif', expt(x).dir, expt(x).name, movParam.binT);
        if ~exist(tempTifName, 'file')
            tempCent = round(vertcat(ROI{x}.cent));
            ZprojStack = pipe.zproj(expt(x).sbx, 1, expt(x).totScan-1, 1, unique(tempCent(:,3)), 'mtype','.sbx_interp', 'registration',false);
            ZprojStack = ZprojStack(segParams{x}.edges(3)+1:end-segParams{x}.edges(4),segParams{x}.edges(1)+1:end-segParams{x}.edges(2),:); % crop edges
            ZprojStack = BinDownMean(permute(ZprojStack,[3,1,2]), movParam.binT); 
            ZprojStack = permute(ZprojStack, [2,3,1]); % Downsample zproj 
            saveastiff(uint16(ZprojStack), tempTifName, GrayOpt);
        else
            ZprojStack = loadtiff(tempTifName);
        end
    end
    toc
    NprojFrame = size(ZprojStack, 3);
    displayLims = [prctile(ZprojStack(:), movParam.displayPct(1)), prctile(ZprojStack(:), movParam.displayPct(2))];
    toc
    % Overlay final ROI over projection stack
    NwriteFrame = NprojFrame; 5; % %100;
    close all; clearvars storeFrames
    storeFrames = repmat(struct('cdata',[], 'colormap',[]), 1, NwriteFrame);
    aviFig = figure('color','k', 'Units','normalized', 'OuterPosition', [0 0 1 1]); % 'WindowState','maximized',

    for z = flip(1:NwriteFrame) %
        cla;
        imshow(ZprojStack(:,:,z), displayLims, 'border','tight'); hold on;
        text(40, size(ZprojStack,1)-20, sprintf('%2.1f s', Tbin(z)), 'color','w' ); % 295
        for roi = stillSumm{x}.fluor.baseUp %1:expt(x).Nroi
            plot( ROI{x}(roi).footprintEdge(:,2)-segParams{x}.edges(1)+1, ROI{x}(roi).footprintEdge(:,1)-segParams{x}.edges(3)+1, '.', 'color',darkBlue, 'MarkerSize',1 ); %  0.5*[1,1,1]
        end
        for roi = stillSumm{x}.fluor.baseDown %1:expt(x).Nroi
            plot( ROI{x}(roi).footprintEdge(:,2)-segParams{x}.edges(1)+1, ROI{x}(roi).footprintEdge(:,1)-segParams{x}.edges(3)+1, '.', 'color',darkRed, 'MarkerSize',1 ); %  0.5*[1,1,1]
        end
        %impixelinfo
        exportgraphics(gca, tempFrameFile , 'Resolution',100); % save image with set resolution
        storeFrames(z) = im2frame(imread(tempFrameFile));  % getframe(aviFig);  % convert image to frame
    end
    close(aviFig);
    delete(tempFrameFile); 
    toc
    
    % Write AVI
    aviPath = sprintf('%s%s_baseEffect.avi', expt(x).dir, expt(x).name);
    fprintf('\nWriting %s...  ', aviPath); tic;
    writerObj = VideoWriter(aviPath);
    writerObj.FrameRate = movParam.aviRate; % set the seconds per image
    open(writerObj);
    for z = 1:NwriteFrame %numel(storeFrames)
        try
            writeVideo(writerObj, storeFrames(z));
        catch
            fprintf('\nx = %i: Failed to add frame %i', x, z);
        end
    end
    close(writerObj);
    toc 
    
end
