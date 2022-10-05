%% Show sensitized units before/after CSD
opt = {[0.06,0.04], [0.07,0.05], [0.07,0.07]};  % {[vert, horz], [bottom, top], [left, right] }
figure('WindowState','maximized');
clearvars sp;
for p = flip(1:4), sp(p) = subtightplot(2,2,p, opt{:}); end %#ok<SAGROW>
for roi = rSensitized{x}
    for runs = 1:expt(x).csd-1
        for b = 1:periBout{x}(runs).Nbout
            subplot(sp(1)); 
            plot( periBout{x}(runs).Ton{b}(periBout{x}(runs).onScan{b}), periBout{x}(runs).fluor{b}(periBout{x}(runs).onScan{b},roi)); hold on;
            subplot(sp(3)); cla;
            plot( periBout{x}(runs).Ton{b}(periBout{x}(runs).onScan{b}), periBout{x}(runs).velocity{b}(periBout{x}(runs).onScan{b})); hold on;
        end
    end
    for runs = expt(x).csd:expt(x).Nruns
        for b = 1:periBout{x}(runs).Nbout
            subplot(sp(2)); cla;
            plot( periBout{x}(runs).Ton{b}(periBout{x}(runs).onScan{b}), periBout{x}(runs).fluor{b}(periBout{x}(runs).onScan{b},roi)); hold on;
            subplot(sp(4)); cla;
            plot( periBout{x}(runs).Ton{b}(periBout{x}(runs).onScan{b}), periBout{x}(runs).velocity{b}(periBout{x}(runs).onScan{b})); hold on;
        end
    end
    subplot(sp(1)); title( sprintf('ROI %i: Pre-CSD', roi )); ylabel('dF/F (z)'); hold off;
    subplot(sp(3)); ylabel('Velocity (cm/s)'); xlabel('Peri-Onset Time (s)');hold off;
    subplot(sp(2)); title('Post-CSD'); hold off; % ylabel('dF/F (z)'); hold off;
    subplot(sp(4)); xlabel('Peri-Onset Time (s)'); hold off;
    pause;
    %clf;
end

%% Show activity/deformation/running around onset/offset, averaged over all bouts per experiment
figPath = strcat(figDir, 'BoutOnsetCSD');
if exist(figPath,'file'), delete(figPath); end
jitterWidth = 0.45;
FS = 8; FW = 'bold'; 
close all; clearvars sp SP;
BoutOnsetEffectPrePost = figure('WindowState','maximized', 'color','w');
opt = {[0.03,0.03], [0.06,0.05], [0.04,0.02]};  % {[vert, horz], [bottom, top], [left, right] }
for x = xCSDbout %intersect( xPresent, xCSD )
    for v = 1:NallVars       
        paddedShortData = cell2padmat({onStruct(x).csd.pre.(allVars{v}).short, onStruct(x).csd.post.(allVars{v}).short});
        paddedLongData = cell2padmat({onStruct(x).csd.pre.(allVars{v}).long, onStruct(x).csd.post.(allVars{v}).long});
        
        sp = subtightplot(2, NallVars, v, opt{:}); cla;
        JitterPlot(paddedShortData, jitterWidth);  % , 'ErrorBar',EB, 'ErrorCap',CS
        hold on;
        text(0.5, 0.95, sprintf('p = %1.4f', onStruct(x).csd.post.(allVars{v}).pShort), 'Units','normalized', 'FontSize',FS, 'HorizontalAlignment','center', 'FontWeight',FW );
        line([0, 3], [0,0], 'color','k', 'lineStyle','-');
        xlim([0.25, 2.75]);
        axis square; 
        if v == 1
            title(sprintf('%s: %s', expt(x).name, allVars{v}), 'Interpreter','none');
            ylabel('Short Effect'); 
        else
            title(allVars{v}, 'Interpreter','none');
        end
        set(gca, 'Xtick',[1,2], 'XtickLabel',{'Pre','Post'}); xtickangle(30);
        sp.YRuler.Exponent = 0;

        SP = subtightplot(2, NallVars, NallVars+v, opt{:});  cla;
        JitterPlot(paddedLongData, jitterWidth); % , 'ErrorBar',EB, 'ErrorCap',CS
        hold on;
        text(0.5, 0.95, sprintf('p = %1.4f', onStruct(x).csd.post.(allVars{v}).pLong), 'Units','normalized', 'FontSize',FS, 'HorizontalAlignment','center', 'FontWeight',FW );
        line([0, 3], [0,0], 'color','k', 'lineStyle','-');
        xlim([0.25, 2.75]);
        axis square; 
        if v == 1, ylabel('Long Effect'); end
        set(gca, 'Xtick',[1,2], 'XtickLabel',{'Pre','Post'}); xtickangle(30);
        SP.YRuler.Exponent = 0;
    end
    %{
    figPath = sprintf('%sBoutOnsetEffectPrePost_%s.tif', figDir, expt(x).name);
    if exist(figPath,'file'), delete(figPath); end
    fprintf('\nSaving %s', figPath);
    print(BoutOnsetEffectPrePost, figPath, '-dtiff' );
    clf;
    %}
    pause;
end

%% Heatmaps of bout onsets pre- and post-CSD
figParams.SubplotOpt = {[0.025, 0.03], [0.06, 0.05], [0.03,0.05]};  % {[vert, horz], [bottom, top], [left, right] }
figParams.FontSize = 8;
figParams.TickLength = [0.003,0];
figParams.vars = allVars;
figParams.lims = viewLims;
figParams.csd = true;
figParams.yStamp = -4;
figParams.dir = 'D:\MATLAB\LevyLab\Figures\3D\CSD\';
figParams.path = strcat(figParams.dir,  'BoutOnsetPrePostCSD.pdf');
figParams.save = false;
for x = xPresent
    [BoutOnsetHeatmaps, figParams] = MakeBoutOnsetHeatmaps(figParams, expt(x), periBout{x}, boutSumm(x), onStruct(x), csdBout{x}(expt(x).csd) );
end

%% Compare GLM coefficients for locomotion sensitive, well-fit units (run onPreGLM and onPostGLM scripts first)
close all;
opt = {[0.025, 0.02], [0.06, 0.05], [0.04,0.03]};  % {[vert, horz], [bottom, top], [left, right] }
GLMcoeffPrePost = figure('WindowState','maximized', 'color','w');
Nsensitized = nan(1, onPrePred(x).N); sensFrac = nan(1, onPrePred(x).N);
for v = 1:onPrePred(x).N
    tempCoeffPrePost = zeros(0,2);
    for x = xCSDbout %xCSD
        tempLocoUnits = unique( [rLocoPre{x}, rLocoPost{x}] ); % intersect( rLocoPre{x}, rLocoPost{x} ); %
        tempCoeffPrePost = vertcat( tempCoeffPrePost, abs([onPreSumm(x).peakCoeff(tempLocoUnits, v), onPostSumm(x).peakCoeff(tempLocoUnits, v)]));
    end
    tempCoeffPrePost( all(isnan(tempCoeffPrePost),2), : ) = [];
    tempCoeffPrePost(isnan(tempCoeffPrePost)) = 0;
    Nsensitized(v) = numel(find( tempCoeffPrePost(:,2)./tempCoeffPrePost(:,1) > 1 ));
    sensFrac(v) = Nsensitized(v)/size(tempCoeffPrePost,1);
    
    subtightplot(1, onPrePred(x).N, v, opt{:});
    plot(tempCoeffPrePost(:,1), tempCoeffPrePost(:,2), '.'); hold on;
    xLimits = get(gca,'Xlim'); yLimits = get(gca,'Ylim');
    line([0,xLimits(2)], [0,xLimits(2)], 'color','k', 'lineStyle','--');
    xlabel('Coeff Magnitude (Pre)'); 
    if v == 1, ylabel('Coeff Magnitude (Post)'); end
    title( sprintf('%s (%2.1f pct)',onPrePred(x).name{v}, 100*sensFrac(v)), 'Interpreter','none')
    axis square;
end
% Save figure;
figPath = sprintf('%sGLMcoeffPrePostCSD.tif', figDir);
if exist(figPath,'file'), delete(figPath); end
fprintf('\nSaving %s', figPath);
print(GLMcoeffPrePost, figPath, '-dtiff' );


%% Compare # of locomotion-sensitive ROI pre vs post CSD
alphaVal = 0.2;
offsetVal = 0.35;
tailVal = 'right';
locoFrac = nan(Nexpt, 2); %locoPreFrac = nan(1,Nexpt); locoPostFrac = nan(1,Nexpt);
LocoSensitivePrePost = figure('WindowState','maximized', 'color','w');
for x = xCSDbout %xCSD
    [~, pPre] = ttest( onStruct(x).fluor.shortEffect(onStruct(x).csd.pre.bout,:), 0, 'tail',tailVal );
    rLocoPre{x} = find(pPre < sigThresh);
    locoFrac(x,1) = numel(rLocoPre{x})/expt(x).Nroi; % locoPreFrac(x) = 
    [~, pPost] = ttest( onStruct(x).fluor.shortEffect(onStruct(x).csd.post.bout,:), 0, 'tail',tailVal );
    rLocoPost{x} = find(pPost < sigThresh);
    locoFrac(x,2) = numel(rLocoPost{x})/expt(x).Nroi; % locoPostFrac(x) = 
    subplot(1,3,[1,2]);
    line( [0,expt(x).Nroi+1], [0,0], 'color','k' ); hold on;
    plot(1:expt(x).Nroi, onStruct(x).fluor.shortEffect(onStruct(x).csd.pre.bout,:), '.', 'Color',alphaVal*[1,1,1] ); 
    errorbar( 1:expt(x).Nroi, mean(onStruct(x).fluor.shortEffect(onStruct(x).csd.pre.bout,:),1,'omitnan'), SEM(onStruct(x).fluor.shortEffect(onStruct(x).csd.pre.bout,:)), 'k.', 'MarkerSize',1 )
    plot( (1:expt(x).Nroi)+offsetVal, onStruct(x).fluor.shortEffect(onStruct(x).csd.post.bout,:), '.', 'Color',alphaVal*[1,0,0] );
    errorbar( (1:expt(x).Nroi)+offsetVal, mean(onStruct(x).fluor.shortEffect(onStruct(x).csd.post.bout,:),1,'omitnan'), SEM(onStruct(x).fluor.shortEffect(onStruct(x).csd.post.bout,:)), 'r.', 'MarkerSize',1 )
    yRange = get(gca,'Ylim');
    for r = 1:expt(x).Nroi
        if pPre(r) < sigThresh, text(r, yRange(2), '*', 'color','k', 'HorizontalAlignment','center'); end
        if pPost(r) < sigThresh, text(r+offsetVal, yRange(2), '*', 'color','r', 'HorizontalAlignment','center'); end
    end
    Npre = sum( pPre < sigThresh );
    Npost = sum( pPost < sigThresh );
    title( sprintf('%s (%i bouts pre-CSD, %i post)', expt(x).name, onStruct(x).csd.pre.Nbout, onStruct(x).csd.post.Nbout), 'Interpreter','none');
    xlim([0,expt(x).Nroi+1]);
    xlabel('ROI'); ylabel('Short-Latency Bout Onset Fluor Effect');
    
    subplot(1,3,3);
    bar(100*[locoPreFrac(x), locoPostFrac(x)]); %bar([100*Npre/expt(x).Nroi, 100*Npost/expt(x).Nroi]);
    ylabel('% locomotion-sensitive');
    ylim([0,100]);
    set(gca, 'XtickLabel',{'Pre','Post'});
    box off;
    
    % Save figure;
    %figPath = sprintf('%sLocoSensitivePrePost_%s.tif', figDir, expt(x).name);
    %if exist(figPath,'file'), delete(figPath); end
    %fprintf('\nSaving %s', figPath);
    %export_fig( figPath, '-pdf', '-painters','-q101', '-append', LocoSensitivePrePost ); pause(1);
    %print(LocoSensitivePrePost, figPath, '-dtiff' ); pause(1);
    pause;
    clf;
end

close all;
LocoSensitivePrePostExpt = figure('WindowState','maximized', 'color','w');
bar(100*locoFrac(xCSDbout,:));
axis square;
set(gca,'Xtick', 1:numel(xCSDbout), 'XtickLabel', {expt(xCSDbout).name}, 'TickLabelInterpreter','none', 'box','off' ); % ,
ylabel('% ROI activated by locomotion onset');
legend('Pre-CSD','Post-CSD'); 

figPath = sprintf('%sLocoSensitivePrePostExpt.tif', figDir);
if exist(figPath,'file'), delete(figPath); end
fprintf('\nSaving %s', figPath);
print(LocoSensitivePrePostExpt, figPath, '-dtiff' );

%% Heatmaps of bout onsets pre- and post-CSD
figDir = 'D:\MATLAB\LevyLab\Figures\3D\CSD';
figPath = strcat(figDir,  'BoutOnsetHeatCSD.pdf');
if exist(figPath,'file'), delete(figPath); end
close all; 
BoutOnsetHeatCSD = figure('WindowState','maximized', 'color','w');
opt = {[0.025, 0.03], [0.06, 0.05], [0.03,0.05]};  % {[vert, horz], [bottom, top], [left, right] }
yText = -4;
FS = 8;
for x = xPresent %intersect(xPresent, xCSD)
    [Nt, Nunit, Nbout] = size(onStruct(x).(allVars{1}).data);
    for v = 1:NallVars%-1
        subEffectData = [];
        for b = 1:Nbout
            subEffectData = [subEffectData; onStruct(x).(allVars{v}).data(:,:,b) - onStruct(x).(allVars{v}).base(:,:,b)];
        end
        if ismember(x, xCSD)
            subEffectData = [subEffectData; mean(onStruct(x).(allVars{v}).data(:,:,onStruct(x).csd.pre.bout) - onStruct(x).(allVars{v}).base(:,:,onStruct(x).csd.pre.bout), 3, 'omitnan'); ...
                mean(onStruct(x).(allVars{v}).data(:,:,onStruct(x).csd.post.bout) - onStruct(x).(allVars{v}).base(:,:,onStruct(x).csd.post.bout), 3, 'omitnan')];
        else
            subEffectData = [subEffectData; mean(onStruct(x).(allVars{v}).data-onStruct(x).(allVars{v}).base, 3, 'omitnan')];
        end
        
        subtightplot(NallVars, 1, v, opt{:});  cla;
        if ~strcmpi(allVars{v}, 'speed')
            imagesc(subEffectData'); hold on;
            axPos = get(gca,'Position');
            CB(v) = colorbar('EastOutside'); CB(v).Label.String = allVars{v}; CB(v).Label.FontWeight = 'bold'; CB(v).Label.Interpreter = 'none';
            set(gca, 'Position',axPos, 'TickDir','out', 'TickLength',[0.003,0], 'Xtick',1:Nt:Nt*Nbout, 'XtickLabel',[], 'box','off' ); % 1:Nbout
        else
            plot( subEffectData );
            ylabel('Speed');
            set(gca, 'TickDir','out', 'TickLength',[0.001,0], 'Xtick',1:Nt:Nt*Nbout, 'XtickLabel',1:Nbout );
            xlim([-inf,inf]);
            xlabel( sprintf('Bouts (%s)', expt(x).name), 'Interpreter','none' );
        end
        yRange = get(gca, 'Ylim');
        for b = 1:Nbout 
            line(b*Nt+0.5*[1,1], yRange, 'color','k');  % [0,Nunit+1]
            if v == 1
                if ismember(x, xCSD)  
                    if b > onStruct(x).csd.pre.bout(end), timeColor = 'r'; else, timeColor = 'k'; end
                    text(0.5+(b-0.5)*Nt, yText, sprintf('%2.1f', (boutSumm(x).Tstart(b)-csdBout{x}(expt(x).csd).Tstart)/60), 'HorizontalAlignment','center', 'FontSize',FS, 'color',timeColor );
                else
                    text(0.5+(b-0.5)*Nt, yText, sprintf('%2.1f', boutSumm(x).Tstart(b)/60), 'HorizontalAlignment','center', 'FontSize',FS, 'color','k' );
                end
            end
        end
        line((Nbout+1)*Nt+0.5*[1,1], yRange, 'color','k'); 

        if v == 1
            text(1, yText, 'Min:', 'HorizontalAlignment','right', 'FontSize',FS, 'color','k' ); 
            if ismember(x, xCSD)
                text(0.5+(Nbout+1-0.5)*Nt, yText, 'Pre', 'HorizontalAlignment','center', 'FontSize',FS, 'color','k' );
                text(0.5+(Nbout+2-0.5)*Nt, yText, 'Post', 'HorizontalAlignment','center', 'FontSize',FS, 'color','r' );
            else
                text(0.5+(Nbout+1-0.5)*Nt, yText, 'Mean', 'HorizontalAlignment','center', 'FontSize',FS, 'color','k' );
            end
            ylabel('ROI');
            set(gca,'Xtick', [1,find(periBout{x}(1).on.T == 0), Nt], 'XtickLabel',sprintfc('%2.1f', periBout{x}(1).on.T([1,find(periBout{x}(1).on.T == 0), end])), 'FontSize',8 );
            %text(-0.1, 1.2, 'Sec', 'HorizontalAlignment','right', 'FontSize',FS, 'color','k', 'Units','normalized' );
        elseif ~strcmpi(allVars{v},'speed')
            ylabel('Plane');
            set(gca,'Xtick', [] );
        end
    end
    
    impixelinfo; pause;
    %export_fig( figPath, '-pdf', '-painters','-q101', '-append', BoutOnsetHeatCSD );  pause(1);
    clf;
end

%% Heatmaps of bout onsets pre- and post-CSD AND CSD itself
figDir = 'D:\MATLAB\LevyLab\Figures\3D\CSD';
figPath = strcat(figDir,  'BoutOnsetCSD.pdf');
if exist(figPath,'file'), delete(figPath); end
close all; 
BoutOnsetCSD = figure('WindowState','maximized', 'color','w');
opt = {[0.025, 0.03], [0.06, 0.05], [0.03,0.05]};  % {[vert, horz], [bottom, top], [left, right] }
FS = 8;
for x = xPresent %xPresent %intersect(xPresent, xCSD)
    [Nt, Nunit, Nbout] = size(onStruct(x).(allVars{1}).data);
    for v = 1:NallVars%-1
        subEffectData = [];
        if ismember(x, xCSD)
            for b = onStruct(x).csd.pre.bout
                subEffectData = [subEffectData; onStruct(x).(allVars{v}).data(:,:,b) - onStruct(x).(allVars{v}).base(:,:,b)]; % pre-CSD bouts
            end
            subEffectData = [subEffectData; onStruct(x).csd.(allVars{v}).data - onStruct(x).csd.(allVars{v}).base]; % CSD bout
            for b = onStruct(x).csd.post.bout
                subEffectData = [subEffectData; onStruct(x).(allVars{v}).data(:,:,b) - onStruct(x).(allVars{v}).base(:,:,b)]; % post CSD bouts
            end
            subEffectData = [subEffectData; mean(onStruct(x).(allVars{v}).data(:,:,onStruct(x).csd.pre.bout) - onStruct(x).(allVars{v}).base(:,:,onStruct(x).csd.pre.bout), 3, 'omitnan'); ...
                mean(onStruct(x).(allVars{v}).data(:,:,onStruct(x).csd.post.bout) - onStruct(x).(allVars{v}).base(:,:,onStruct(x).csd.post.bout), 3, 'omitnan')]; % pre and post-CSD means
        else
            for b = 1:Nbout
                subEffectData = [subEffectData; onStruct(x).(allVars{v}).data(:,:,b) - onStruct(x).(allVars{v}).base(:,:,b)];
            end
            subEffectData = [subEffectData; mean(onStruct(x).(allVars{v}).data-onStruct(x).(allVars{v}).base, 3, 'omitnan')];
        end
        
        subtightplot(NallVars, 1, v, opt{:});  cla;
        if ~strcmpi(allVars{v}, 'speed')
            imagesc(subEffectData'); hold on;
            axPos = get(gca,'Position');
            CB(v) = colorbar('EastOutside'); CB(v).Label.String = allVars{v}; CB(v).Label.FontWeight = 'bold'; CB(v).Label.Interpreter = 'none';
            set(gca, 'Position',axPos, 'TickDir','out', 'TickLength',[0.003,0], 'Xtick',1:Nt:Nt*Nbout, 'XtickLabel',[], 'box','off' ); % 1:Nbout
        else
            plot( subEffectData );
            ylabel('Speed');
            set(gca, 'TickDir','out', 'TickLength',[0.001,0], 'Xtick',1:Nt:Nt*Nbout, 'XtickLabel',1:Nbout );
            xlim([-inf,inf]);
            xlabel( sprintf('Bouts (%s)', expt(x).name), 'Interpreter','none' );
        end
        yRange = get(gca, 'Ylim');
        % Add lines between bouts
        for b = 1:size(subEffectData,1)/Nt,  line(b*Nt+0.5*[1,1], yRange, 'color','k');  end
        % Add time stamps
        if v == 1
            ylabel('ROI');
            set(gca,'Xtick', [1,find(periBout{x}(1).on.T == 0), Nt], 'XtickLabel',sprintfc('%2.1f', periBout{x}(1).on.T([1,find(periBout{x}(1).on.T == 0), end])), 'FontSize',8 );
            text(1, yText, 'Min:', 'HorizontalAlignment','right', 'FontSize',FS, 'color','k' ); 
            if ismember(x, xCSD)  
                for b = onStruct(x).csd.pre.bout
                    text(0.5+(b-0.5)*Nt, yText, sprintf('%2.1f', (boutSumm(x).Tstart(b)-csdBout{x}(expt(x).csd).Tstart)/60), 'HorizontalAlignment','center', 'FontSize',FS, 'color','k' );
                end
                text(0.5+(numel(onStruct(x).csd.pre.bout)+0.5)*Nt, yText, 'CSD', 'HorizontalAlignment','center', 'FontSize',FS, 'color','k' );
                for b = onStruct(x).csd.post.bout
                    text(0.5+(b+0.5)*Nt, yText, sprintf('%2.1f', (boutSumm(x).Tstart(b)-csdBout{x}(expt(x).csd).Tstart)/60), 'HorizontalAlignment','center', 'FontSize',FS, 'color','r' );
                end
                text(0.5+(Nbout+1.5)*Nt, yText, 'Pre', 'HorizontalAlignment','center', 'FontSize',FS, 'color','k' );
                text(0.5+(Nbout+2.5)*Nt, yText, 'Post', 'HorizontalAlignment','center', 'FontSize',FS, 'color','r' );
            else
                for b = 1:Nbout, text(0.5+(b-0.5)*Nt, yText, sprintf('%2.1f', boutSumm(x).Tstart(b)/60), 'HorizontalAlignment','center', 'FontSize',FS, 'color','k' ); end
                text(0.5+(Nbout+0.5)*Nt, yText, 'Mean', 'HorizontalAlignment','center', 'FontSize',FS, 'color','k' );
            end
        elseif ~strcmpi(allVars{v},'speed')
            ylabel('Plane');
            set(gca,'Xtick', [] );
        end
    end
    
    %impixelinfo; pause;
    export_fig( figPath, '-pdf', '-painters','-q101', '-append', BoutOnsetCSD );  pause(1);
    clf;
end

%% Heatmaps of bout onsets pre- and post-CSD WITHOUT CSD itself
figDir = 'D:\MATLAB\LevyLab\Figures\3D\';

close all; 
BoutOnsetPrePost = figure('WindowState','maximized', 'color','w');
opt = {[0.025, 0.03], [0.06, 0.05], [0.03,0.05]};  % {[vert, horz], [bottom, top], [left, right] }
FS = 8;
for x = xPresent %intersect(xPresent, xCSD)
    [Nt, Nunit, Nbout] = size(onStruct(x).(allVars{1}).data);
    for v = 1:NallVars%-1
        subEffectData = [];
        if ismember(x, xCSD)
            for b = onStruct(x).csd.pre.bout
                subEffectData = [subEffectData; onStruct(x).(allVars{v}).data(:,:,b) - onStruct(x).(allVars{v}).base(:,:,b)]; % pre-CSD bouts
            end
            %subEffectData = [subEffectData; onStruct(x).csd.(allVars{v}).data - onStruct(x).csd.(allVars{v}).base]; % CSD bout
            for b = onStruct(x).csd.post.bout
                subEffectData = [subEffectData; onStruct(x).(allVars{v}).data(:,:,b) - onStruct(x).(allVars{v}).base(:,:,b)]; % post CSD bouts
            end
            subEffectData = [subEffectData; mean(onStruct(x).(allVars{v}).data(:,:,onStruct(x).csd.pre.bout) - onStruct(x).(allVars{v}).base(:,:,onStruct(x).csd.pre.bout), 3, 'omitnan'); ...
                mean(onStruct(x).(allVars{v}).data(:,:,onStruct(x).csd.post.bout) - onStruct(x).(allVars{v}).base(:,:,onStruct(x).csd.post.bout), 3, 'omitnan')]; % pre and post-CSD means
        else
            for b = 1:Nbout
                subEffectData = [subEffectData; onStruct(x).(allVars{v}).data(:,:,b) - onStruct(x).(allVars{v}).base(:,:,b)];
            end
            subEffectData = [subEffectData; mean(onStruct(x).(allVars{v}).data-onStruct(x).(allVars{v}).base, 3, 'omitnan')];
        end
        
        subtightplot(NallVars, 1, v, opt{:});  cla;
        if ~strcmpi(allVars{v}, 'speed')
            imagesc(subEffectData'); hold on;
            axPos = get(gca,'Position');
            CB(v) = colorbar('EastOutside'); CB(v).Label.String = allVars{v}; CB(v).Label.FontWeight = 'bold'; CB(v).Label.Interpreter = 'none';
            set(gca, 'Position',axPos, 'TickDir','out', 'TickLength',[0.003,0], 'Xtick',1:Nt:Nt*Nbout, 'XtickLabel',[], 'box','off' ); % 1:Nbout
        else
            plot( subEffectData );
            ylabel('Speed');
            set(gca, 'TickDir','out', 'TickLength',[0.001,0], 'Xtick',1:Nt:Nt*Nbout, 'XtickLabel',1:Nbout );
            xlim([-inf,inf]);
            xlabel( sprintf('Bouts (%s)', expt(x).name), 'Interpreter','none' );
        end
        yRange = get(gca, 'Ylim');
        % Add lines between bouts
        for b = 1:size(subEffectData,1)/Nt,  line(b*Nt+0.5*[1,1], yRange, 'color','k');  end
        % Add time stamps
        if v == 1
            ylabel('ROI');
            set(gca,'Xtick', [1,find(periBout{x}(1).on.T == 0), Nt], 'XtickLabel',sprintfc('%2.1f', periBout{x}(1).on.T([1,find(periBout{x}(1).on.T == 0), end])), 'FontSize',8 );
            text(1, yText, 'Min:', 'HorizontalAlignment','right', 'FontSize',FS, 'color','k' ); 
            if ismember(x, xCSD)  
                for b = onStruct(x).csd.pre.bout
                    text(0.5+(b-0.5)*Nt, yText, sprintf('%2.1f', (boutSumm(x).Tstart(b)-csdBout{x}(expt(x).csd).Tstart)/60), 'HorizontalAlignment','center', 'FontSize',FS, 'color','k' );
                end
                %text(0.5+(numel(onStruct(x).csd.pre.bout)+0.5)*Nt, yText, 'CSD', 'HorizontalAlignment','center', 'FontSize',FS, 'color','k' );
                for b = onStruct(x).csd.post.bout
                    text(0.5+(b-0.5)*Nt, yText, sprintf('%2.1f', (boutSumm(x).Tstart(b)-csdBout{x}(expt(x).csd).Tstart)/60), 'HorizontalAlignment','center', 'FontSize',FS, 'color','r' );
                end
                text(0.5+(Nbout+0.5)*Nt, yText, 'Pre', 'HorizontalAlignment','center', 'FontSize',FS, 'color','k' );
                text(0.5+(Nbout+1.5)*Nt, yText, 'Post', 'HorizontalAlignment','center', 'FontSize',FS, 'color','r' );
            else
                for b = 1:Nbout, text(0.5+(b-0.5)*Nt, yText, sprintf('%2.1f', boutSumm(x).Tstart(b)/60), 'HorizontalAlignment','center', 'FontSize',FS, 'color','k' ); end
                text(0.5+(Nbout+0.5)*Nt, yText, 'Mean', 'HorizontalAlignment','center', 'FontSize',FS, 'color','k' );
            end
        elseif ~strcmpi(allVars{v},'speed')
            ylabel('Plane');
            set(gca,'Xtick', [] );
        end
    end
    
    impixelinfo; pause;
    
    %{
    %export_fig( figPath, '-pdf', '-painters','-q101', '-append', BoutOnsetPrePost );  pause(1);
    figPath = sprintf('%sBoutOnsetPrePost_%s.tif', figDir, expt(x).name);
    if exist(figPath,'file'), delete(figPath); end
    fprintf('\nSaving %s', figPath);
    print(BoutOnsetPrePost, figPath, '-dtiff' );
    %}
    clf;
end

%% How does Fo change post CSD?
preDur = 15; % use the last preDur minute of the run before the CSD run as baseline
binWidth = 15*60; %10; % bin data by binWidth seconds
alphaLevel = 0.05;
close all; clearvars h g;
figure('Units','normalized', 'OuterPosition', [0,0,1,1], 'Color','w', 'PaperOrientation','landscape');
FS = 18;
FpreBin = cell(1,Nexpt); TpreBin = cell(1,Nexpt); FpostBin = cell(1,Nexpt); TpostBin = cell(1,Nexpt);  Tbin = cell(1,Nexpt); FnormBin = cell(1,Nexpt); dFnorm = cell(1,Nexpt);
dimFrac = cell(1,Nexpt); pDim = cell(1,Nexpt); testROI = repmat(struct('non',[], 'Nnon',NaN, 'dim',[], 'Ndim',NaN, 'Nbright',NaN, 'bright',[]), 1, Nexpt);
clustMat = zeros(0, 4); %[2-17 mins, 47-62 mins, x, roi, T1, T2]
for x = setdiff(xPresent, 36)
    % Concatenate variables across runs
    Tcat = vertcat(Tscan{x}{:}) - csdBout{x}(expt(x).csd).Tstart;
    Tcat = Tcat(1:end-1);
    fluorCatROI = [fluor{x}(:).Froi];
    tempFroi = vertcat(fluorCatROI.ROI);
    tempFroi = tempFroi(1:end-1,:);
    fluorCatNP = [fluor{x}(:).Fnp];
    tempFnp = vertcat(fluorCatNP.ROI);
    tempFnp = tempFnp(1:end-1,:);
    
    % Subtract and filter signals 
    tempFsub = tempFroi - tempFnp + median(tempFnp, 1, 'omitnan');
    lpFreq = 2;
    if lpFreq > expt(x).scanRate, lpFreq = expt(x).scanRate/(2*pi); end
    gaussSigma = 1/(2*pi*lpFreq);
    gaussFilt = MakeGaussFilt( 4, 0, gaussSigma, expt(x).scanRate, false ); %  gaussFilt = MakeGaussFilt( 4, 0, 1/expt.scanRate, expt.scanRate, true );
    tempFsub = filtfilt(gaussFilt, 1, tempFsub);

    % Pull out the fluorescence signals, pre-pinprick and post-CSD, bin and normalize
    binSize = round(binWidth*expt(x).scanRate);
    preScans = 1:expt(x).scanLims(expt(x).csd);
    [preScanChunks,~, preChunkLength] = MakeChunkLims(preScans(1), preScans(end), preScans(end), 'size',binSize);
    preScanChunks(preChunkLength < binSize,:) = [];
    for c = 1:size(preScanChunks,1)
        FpreBin{x}(c,:) = median(tempFsub(preScanChunks(c,1):preScanChunks(c,2),:), 1, 'omitnan');
        TpreBin{x}(c,1) = median(Tcat(preScanChunks(c,1):preScanChunks(c,2)), 1, 'omitnan');
    end
    postScans = find(Tcat/60 > 2);
    [postScanChunks,~, postChunkLength] = MakeChunkLims(postScans(1), postScans(end), postScans(end), 'size',binSize);
    postScanChunks(postChunkLength < binSize/2,:) = [];
    for c = 1:size(postScanChunks,1)
        FpostBin{x}(c,:) = median(tempFsub(postScanChunks(c,1):postScanChunks(c,2),:), 1, 'omitnan');
        TpostBin{x}(c,1) = median(Tcat(postScanChunks(c,1):postScanChunks(c,2)), 1, 'omitnan');
    end

    Tbin{x} = vertcat(TpreBin{x}, TpostBin{x})/60; % convert from seconds to minutes
    FnormBin{x} = vertcat(FpreBin{x}, FpostBin{x})./FpreBin{x}(end,:);
    postFirst = find(Tbin{x} > 0, 1, 'first');
    [~,post60] = min(abs(Tbin{x}-55)); % find(Tbin{x} > 0 & Tbin{x} < 62, 1, 'last');
    dFnorm{x} = (FnormBin{x}([postFirst,post60],:)-1)';  % relative changes at ~15 mins and ~60 mins post-CSD
       
    % Calcuate the fraction of ROI that are dimmer than just before CSD
    dimFrac{x} = sum(FnormBin{x} < 1, 2)/expt(x).Nroi;
    dimFrac{x}( find(Tbin{x} < 0, 1, 'last') ) = NaN;
    for roi = flip(1:expt(x).Nroi)
        [~,pDim{x}(roi)] = ttest2( FnormBin{x}(1:postFirst-1,roi), FnormBin{x}(postFirst:post60,roi) );
        if pDim{x}(roi) > alphaLevel
            testROI(x).non = [testROI(x).non, roi];
        elseif mean(FnormBin{x}(postFirst:post60,roi)) < mean(FnormBin{x}(1:postFirst-1,roi))
            testROI(x).dim = [testROI(x).dim, roi];
        else
            testROI(x).bright = [testROI(x).bright, roi];
        end
    end
    testROI(x).Nnon = numel(testROI(x).non);
    testROI(x).Ndim = numel(testROI(x).dim);
    testROI(x).Nbright = numel(testROI(x).bright);
    
    clustMat = vertcat(clustMat, [dFnorm{x}, repmat(x,expt(x).Nroi,1), (1:expt(x).Nroi)']); %Tbin{x}(postFirst)
    
    %{
    %roiRand = datasample( 1:expt(x).Nroi, 5, 'Replace',false);
    for roi = 1:expt(x).Nroi  %57 % roiRand% 170 % 
        subplot(2,1,1); cla;
        h(1) = plot( Tcat/60, tempFroi(:,roi) ); hold all;
        h(2) = plot( Tcat/60, tempFnp(:,roi) );
        ylim([-Inf,Inf]); xlim([-Inf,Inf]);
        legend(h, {'ROI','Neuropil'}, 'AutoUpdate','off');
        ylabel('Raw Fluor');
        title(sprintf('[x, run, roi] = [%i, %i]: Raw Data', x, roi));
        
        subplot(2,1,2); cla;
        g(1) = plot( Tcat/60, Fnorm{x}(:,roi) ); hold on; 
        g(2) = plot( Tbin/60, tempFbin(:,roi), 'k', 'LineWidth',1.5 ); hold on; % , 'color',[0,0,0,0.01]
        legend(g, {'Full','Binned'}, 'AutoUpdate','off');
        line(Tcat(tempPreScans(1))*[1,1]/60, [min(Fnorm{x}(:,roi)), max(Fnorm{x}(:,roi))], 'color','r', 'lineStyle','--');
        line(Tcat(tempPreScans(end))*[1,1]/60, [min(Fnorm{x}(:,roi)), max(Fnorm{x}(:,roi))], 'color','r', 'lineStyle','--');
        ylim([-Inf,Inf]); xlim([-Inf,Inf]);
        title('Processed Signal');
        ylabel('Normalized, Neuropil-subtracted Fluorescence');
        xlabel('Peri-CSD Time (min)');
        pause;
    end
    %}
    %{
    figure;
    plot( Tbin/60, tempFbin, 'color',[0,0,0,0.04], 'LineWidth',1.5 ); hold on;
    ylim([-Inf,Inf]); xlim([-Inf,Inf]);
    title('Processed Signal');
    ylabel('Normalized, Neuropil-subtracted Fluorescence');
    xlabel('Peri-CSD Time (min)');
    cla;
    %}
    %{
    imagesc(tempFbin(Tbin/60>=-10 & Tbin/60 < 30,:)')
    caxis([0.5, 1.5]); 
    CB = colorbar;
    CB.Label.String = 'Norm. Fluor';
    ylabel('ROI');
    tempXlim = get(gca, 'Xlim');
    set(gca,'Xtick',tempXlim, 'XtickLabel',[-10,30], 'TickDir','out', 'FontSize',FS);
    xlabel('Peri-CSD Time (min)');
    title(sprintf('x = %i: %s', x, expt(x).name), 'Interpreter','none');
    impixelinfo
    pause;
    %}
    %{
    subplot(2,1,1); cla;
    line([-60,120], [1,1], 'color','r'); hold on;
    plot( Tbin{x}, FnormBin{x}(:,testROI(x).non), 'color',[0,0,0,0.1] ); hold on;
    plot( Tbin{x}, FnormBin{x}(:,testROI(x).dim), 'color',[1,0,0,0.1] );
    plot( Tbin{x}, FnormBin{x}(:,testROI(x).bright), 'color',[0,0,1,0.1] )
    plot( Tbin{x}, median(FnormBin{x}, 2), 'k.', 'markersize',10)
    set(gca,'Xtick',[-60:10:60]);
    xlabel('Time post-CSD (min)'); ylabel('Baseline-Normalized Median Fluorescence');
    title(sprintf('%s', expt(x).name), 'Interpreter','none');
    %hold off;
    
    subplot(2,1,2); cla;
    line([-60,60], 0.5*[1,1], 'color','r'); hold on;% hold all;
    plot( Tbin{x}, dimFrac{x} ); hold on;
    set(gca,'Xtick',[-60:10:120]);
    ylim([0,1]);
    xlabel('Time post-CSD (min)'); ylabel('Fraction of units dimmed');
    pause;
    %}
    %{
    subplot(2,2,[2,4]);  cla;
    %{
    bar([testROI(x).Nnon, testROI(x).Ndim, testROI(x).Nbright]/expt(x).Nroi ); 
    ylabel('Fraction of ROI'); 
    set(gca, 'Xtick',1:3, 'XtickLabel',{'No effect','Dimmed','Brightened'});
    xtickangle(30);
    axis square;
    title('Effect up to 1 hour post-CSD');
    %}
    % {
    plot( FnormBin{x}(postFirst,:)-1, FnormBin{x}(post60,:)-1, '.' ); hold on;
    line([-0.3,0.3],[0,0], 'color','k');
    line([0,0],[-0.3,0.3], 'color','k');
    xlabel('\DeltaFnorm 0-15 mins post-CSD'); ylabel('\DeltaFnorm 45-60 mins post-CSD');
    axis square;
    %set(gca,'XaxisLocation','origin','YaxisLocation','origin');
    %}

end

% Cluster pooled data
[clustMat(:,5), clustCent] = kmeans(clustMat(:,1:2), 3, 'Replicates',10);
cDim = 1; cNon = 2; cBright = 3;
clustColor = distinguishable_colors(3); %flip(,1);
%plot( clustMat(:,1), clustMat(:,2), '.' )
%% Unpool clustering results and compare to t-testing
clustROI = repmat(struct('non',[], 'Nnon',NaN, 'dim',[], 'Ndim',NaN, 'Nbright',NaN, 'bright',[]), 1, Nexpt);
clustFrac = nan(3,Nexpt); testFrac = nan(3,Nexpt);
for x = setdiff(xPresent, 36)
    tempMat = clustMat(clustMat(:,3) == x,:);
    clustROI(x).non = tempMat(tempMat(:,5) == cNon,4)';
    clustROI(x).Nnon = numel(clustROI(x).non);
    clustROI(x).dim = tempMat(tempMat(:,5) == cDim,4)';
    clustROI(x).Ndim = numel(clustROI(x).dim);
    clustROI(x).bright = tempMat(tempMat(:,5) == cBright,4)';
    clustROI(x).Nbright = numel(clustROI(x).bright);
    clustFrac(:,x) = [clustROI(x).Nnon, clustROI(x).Ndim, clustROI(x).Nbright]/expt(x).Nroi;
    testFrac(:,x) = [testROI(x).Nnon, testROI(x).Ndim, testROI(x).Nbright]/expt(x).Nroi;
end
overallClustFrac = [sum([clustROI(setdiff(xPresent, 36)).Nnon]), sum([clustROI(setdiff(xPresent, 36)).Ndim]), sum([clustROI(setdiff(xPresent, 36)).Nbright])]/sum([expt(setdiff(xPresent, 36)).Nroi]);
overallTestFrac = [sum([testROI(setdiff(xPresent, 36)).Nnon]), sum([testROI(setdiff(xPresent, 36)).Ndim]), sum([testROI(setdiff(xPresent, 36)).Nbright])]/sum([expt(setdiff(xPresent, 36)).Nroi]);

close all; clearvars h; k =0;
figure;
subplot(2,2,1);
line([-0.5,1],[0,0], 'color','k'); hold on;
line([0,0],[-0.5,1], 'color','k');
for c = [cNon, cDim, cBright]
    k = k+1;
    tempInd = find(clustMat(:,5) == c);
    plot( clustMat(tempInd,1), clustMat(tempInd,2), '.', 'color',clustColor(c,:) ); hold on;
    %plot( clustCent(c,1), clustCent(c,2), 'k.', 'MarkerSize',10 )
    h(k) = plot(1,NaN, '.', 'color',clustColor(c,:));
    %pause;
end
xlabel('\DeltaFnorm 0-15 mins post-CSD'); ylabel('\DeltaFnorm 45-60 mins post-CSD'); title('3-means clustering');
axis square;
legend(h, {'No effect','Dimmed','Brightened'}, 'AutoUpdate','off');

subplot(2,2,3);
for x = setdiff(xPresent, 36)
    plot(1:3, clustFrac(:,x)); hold on;
    plot(1:3, clustFrac(:,x), '.');
end
plot(1:3, overallClustFrac, 'k', 'LineWidth',2);
plot(1:3, overallClustFrac, 'k.');
axis square;
xlim([0.5,3.5]); ylim([0,1]);
set(gca,'Xtick',1:3, 'XtickLabel',{'No effect','Dimmed','Brightened'}); xtickangle(30);
ylabel('Fraction of ROI');

subplot(2,2,2);
line([-0.5,1],[0,0], 'color','k'); hold on;
line([0,0],[-0.5,1], 'color','k');
for x = setdiff(xPresent, 36)
    plot( dFnorm{x}(testROI(x).non,1), dFnorm{x}(testROI(x).non,2), '.', 'color',clustColor(cNon,:) );
    plot( dFnorm{x}(testROI(x).dim,1), dFnorm{x}(testROI(x).dim,2), '.', 'color',clustColor(cDim,:) );
    plot( dFnorm{x}(testROI(x).bright,1), dFnorm{x}(testROI(x).bright,2), '.', 'color',clustColor(cBright,:) );
end
xlabel('\DeltaFnorm 0-15 mins post-CSD'); ylabel('\DeltaFnorm 45-60 mins post-CSD'); title('T-testing');
axis square;

subplot(2,2,4);
for x = setdiff(xPresent, 36)
    plot(1:3, testFrac(:,x)); hold on;
    plot(1:3, testFrac(:,x), '.');
end
plot(1:3, overallTestFrac, 'k', 'LineWidth',2);
plot(1:3, overallTestFrac, 'k.');
axis square;
xlim([0.5,3.5]); ylim([0,1]);
set(gca,'Xtick',1:3, 'XtickLabel',{'No effect','Dimmed','Brightened'}); xtickangle(30);
ylabel('Fraction of ROI');

%%
%%  View ROIs
for x = 37 %xPresent
    ViewROI3D(expt(x), ROI{x}, fluor{x}, loco{x}, 'save',false, 'setROI',testROI(x).bright);
end