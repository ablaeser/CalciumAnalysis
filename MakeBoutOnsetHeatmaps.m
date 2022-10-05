function [BoutOnsetHeatmaps, figParams] = MakeBoutOnsetHeatmaps(figParams, expt, periBout, boutSumm, onStruct, csdBout )
%figParams, expt(x), periBout{x}, boutSumm(x), onStruct(x), csdBout{x}(expt(x).csd)

% Heatmaps of bout onsets pre- and post-CSD, optionally including CSD itself
%if exist(figParams.path,'file') && figParams.save, delete(figParams.path); end
close all;  clearvars sp;
BoutOnsetHeatmaps = figure('WindowState','maximized', 'color','w');
[Nscan, ~, Nbout] = size(onStruct.(figParams.vars{1}).data);
boutTicks = 1:Nscan:Nscan*Nbout;
limFields = fieldnames(figParams.lims);
if expt.Nplane == 30
    planeTicks = [1,5:5:30];
elseif expt.Nplane == 15
    planeTicks = [1:3:15,15];
elseif expt.Nplane == 1
    planeTicks = 1;
    figParams.vars( strcmpi(figParams.vars,'shiftZ') ) = [];
    figParams.vars( strcmpi(figParams.vars,'dShiftZ') ) = [];
end
Nvars = numel(figParams.vars);
for v = 1:Nvars
    % Get variable's data
    varEffectData = [];
    if figParams.csd  %ismember(x, xCSD)
        for b = onStruct.csd.pre.bout
            varEffectData = [varEffectData; onStruct.(figParams.vars{v}).data(:,:,b) - onStruct.(figParams.vars{v}).base(:,:,b)]; % pre-CSD bouts
        end
        varEffectData = [varEffectData; onStruct.csd.(figParams.vars{v}).data - onStruct.csd.(figParams.vars{v}).base]; % CSD bout
        for b = onStruct.csd.post.bout
            varEffectData = [varEffectData; onStruct.(figParams.vars{v}).data(:,:,b) - onStruct.(figParams.vars{v}).base(:,:,b)]; % post CSD bouts
        end
        varEffectData = [varEffectData; mean(onStruct.(figParams.vars{v}).data(:,:,onStruct.csd.pre.bout) - onStruct.(figParams.vars{v}).base(:,:,onStruct.csd.pre.bout), 3, 'omitnan'); ...
            mean(onStruct.(figParams.vars{v}).data(:,:,onStruct.csd.post.bout) - onStruct.(figParams.vars{v}).base(:,:,onStruct.csd.post.bout), 3, 'omitnan')]; % pre and post-CSD means
    else
        for b = 1:Nbout
            varEffectData = [varEffectData; onStruct.(figParams.vars{v}).data(:,:,b) - onStruct.(figParams.vars{v}).base(:,:,b)];
        end
        varEffectData = [varEffectData; mean(onStruct.(figParams.vars{v}).data-onStruct.(figParams.vars{v}).base, 3, 'omitnan')];
    end
    
    % Make the plot
    sp(v) = subtightplot(Nvars, 1, v, figParams.SubplotOpt{:});  cla;
    if strcmpi(figParams.vars{v}, 'speed')
        plot( varEffectData );
        ylabel('Speed');
        set(gca, 'TickDir','out', 'TickLength',[0.001,0], 'Xtick',boutTicks, 'XtickLabel',1:Nbout );
        xlim([-inf,inf]);
        xlabel( sprintf('Bouts (%s)', expt.name), 'Interpreter','none' );
    elseif strcmpi(figParams.vars{v}, 'fluor')
        MakeDeformPlot( varEffectData, figParams.vars{v}, figParams.TickLength, boutTicks, unique([1:round(expt.Nroi/5):expt.Nroi, expt.Nroi]), figParams.lims.fluor ) % planeTicks unique([1:10:expt.Nroi,expt.Nroi])
        ylabel('ROI');
    else
        for f = 1:numel(limFields)
            if contains( figParams.vars{v}, limFields{f}, 'ignoreCase',true )
                MakeDeformPlot( varEffectData, figParams.vars{v}, figParams.TickLength, boutTicks, planeTicks, figParams.lims.(limFields{f}) ) 
            end
        end       
    end
    
    %{
    if ~strcmpi(figParams.vars{v}, 'speed')
        imagesc(varEffectData'); hold on;
        axPos = get(gca,'Position');
        CB(v) = colorbar('EastOutside'); CB(v).Label.String = figParams.vars{v}; CB(v).Label.FontWeight = 'bold'; CB(v).Label.Interpreter = 'none';
        set(gca, 'Position',axPos, 'TickDir','out', 'TickLength',[0.003,0], 'Xtick',boutTicks, 'XtickLabel',[], 'box','off' ); % 1:Nbout
        ylabel('Plane');
        set(gca,'Xtick', [] );
        for f = 1:numel(limFields)
            if contains( figParams.vars{v}, limFields{f}, 'ignoreCase',true )
                caxis(figParams.lims.(limFields{f}))
            end
        end
    else
        plot( varEffectData );
        ylabel('Speed');
        set(gca, 'TickDir','out', 'TickLength',[0.001,0], 'Xtick',boutTicks, 'XtickLabel',1:Nbout );
        xlim([-inf,inf]);
        xlabel( sprintf('Bouts (%s)', expt.name), 'Interpreter','none' );
    end
    %}
    % Add lines between bouts
    yRange = get(gca, 'Ylim');
    if isinf(yRange(1)), yRange(1) = min(varEffectData); end
    if isinf(yRange(2)), yRange(2) = max(varEffectData); end
    
    for b = 1:size(varEffectData,1)/Nscan,  line(b*Nscan+0.5*[1,1], yRange, 'color','k');  end
    if v == 1
        set(gca,'Xtick', [1,find(periBout(1).on.T == 0), Nscan], 'XtickLabel',sprintfc('%2.1f', periBout(1).on.T([1,find(periBout(1).on.T == 0), end])), 'FontSize',8 );
        % Add time stamps
        text(1, figParams.yStamp, 'Min:', 'HorizontalAlignment','right', 'FontSize',figParams.FontSize, 'color','k' ); 
        if figParams.csd   
            for b = onStruct.csd.pre.bout
                text(0.5+(b-0.5)*Nscan, figParams.yStamp, sprintf('%2.1f', (boutSumm.Tstart(b)-csdBout.Tstart)/60), 'HorizontalAlignment','center', 'FontSize',figParams.FontSize, 'color','k' );
            end
            text(0.5+(numel(onStruct.csd.pre.bout)+0.5)*Nscan, figParams.yStamp, 'CSD', 'HorizontalAlignment','center', 'FontSize',figParams.FontSize, 'color','r' );
            for b = onStruct.csd.post.bout
                text(0.5+(b+0.5)*Nscan, figParams.yStamp, sprintf('%2.1f', (boutSumm.Tstart(b)-csdBout.Tstart)/60), 'HorizontalAlignment','center', 'FontSize',figParams.FontSize, 'color','r' );
            end
            text(0.5+(Nbout+1.5)*Nscan, figParams.yStamp, 'Pre', 'HorizontalAlignment','center', 'FontSize',figParams.FontSize, 'color','k' );
            text(0.5+(Nbout+2.5)*Nscan, figParams.yStamp, 'Post', 'HorizontalAlignment','center', 'FontSize',figParams.FontSize, 'color','r' );
        else
            for b = 1:Nbout, text(0.5+(b-0.5)*Nscan, figParams.yStamp, sprintf('%2.1f', boutSumm.Tstart(b)/60), 'HorizontalAlignment','center', 'FontSize',figParams.FontSize, 'color','k' ); end
            text(0.5+(Nbout+0.5)*Nscan, figParams.yStamp, 'Mean', 'HorizontalAlignment','center', 'FontSize',figParams.FontSize, 'color','k' );
        end
    end
end
linkaxes(sp,'x');

% Save the figure
if figParams.save
    fprintf('\nSaving %s..', figParams.path );
    export_fig( figParams.path, '-pdf', '-painters','-q101', '-append', BoutOnsetHeatmaps );  pause(1);
else
    impixelinfo; pause;
end
end

