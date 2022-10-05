function InspectPeriOnOff3D(expt, periBout, periStat, showVars, varargin)
% Visualize deformation variables associated with locomotion bout onset/offset
IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'periBout', @isstruct )
addRequired( IP, 'periStat', @isstruct )
addRequired( IP, 'defVars', @iscell )
addOptional( IP, 'limits', [], @isstruct)
%addParameter(IP, 'bSet',1:periBout.Nbout, @isnumeric )
addParameter(IP, 'type', 'both', @ischar )
addParameter(IP, 'show', 'mean', @ischar )
addParameter(IP, 'tick', 1, @isnumeric )
addParameter(IP, 'path', '', @ischar )
parse( IP, expt, periBout, periStat, showVars, varargin{:} ); 
limStruct = IP.Results.limits;
%bSet =  IP.Results.bSet;
figType = IP.Results.type;
showType = IP.Results.show;
tickInt = IP.Results.tick; %2;
figPath = IP.Results.path;
if isempty(limStruct)
    limStruct.z = [-1,3];
    limStruct.trans = [-10, 10]; 
    limStruct.scale = [0.90, 1.1];
    limStruct.shear = [-0.05, 0.05];
    limStruct.shift = [-3, 3];
end

% Onset xticks
onZeroTick = find(periBout(1).on.T == 0);
onTicks = unique([flip(onZeroTick:-tickInt:1), onZeroTick:tickInt:numel(periBout(1).on.ind)]);
TonTicks = arrayfun(  @(x) sprintf('%3.1f',x), periBout(1).on.T(onTicks),  'uniformoutput', false); 
% Offset xticks
Toff = periBout(1).off.T(periBout(1).off.ind);
offZeroTick = find(Toff==0);
offTicks = unique([flip(offZeroTick:-tickInt:1), offZeroTick:tickInt:numel(periBout(1).off.ind)]);
ToffTicks = arrayfun(  @(x) sprintf('%3.1f',x), Toff(offTicks),  'uniformoutput', false);

Nsp = numel(showVars)+1; %10;
opt = {[0.01,0.07], [0.07,0.04], [0.2,0.2]};  % {[vert, horz], [bottom, top], [left, right] }
close all; clearvars sp SP;
if strcmpi(figType, 'on') 
    onFig = figure('WindowState','maximized', 'Color','w', 'PaperOrientation','landscape');
elseif strcmpi(figType, 'off') 
    offFig = figure('WindowState','maximized', 'Color','w', 'PaperOrientation','landscape');
else
    onFig = figure('units','normalized', 'Color','w', 'PaperOrientation','portrait', 'OuterPosition',[0,0,0.5,1]); % outerposition  = [left, bottom, width, height]
    offFig = figure('units','normalized', 'Color','w', 'PaperOrientation','portrait', 'OuterPosition',[0.5,0,0.5,1]);
end

subInd.on.pre = find(periBout(1).on.T < -1);
subInd.on.short = find(periBout(1).on.T > 0 & periBout(1).on.T < 2);
subInd.on.long = find(periBout(1).on.T >= 2);
% INDIVIDUAL BOUTS
if strcmpi(showType, 'ind')
    for runs = 1:numel(periBout)
        if periBout(runs).Nbout > 0 %&& ~isempty( bSet )
            for b = 1:periBout(runs).Nbout % bSet
                % Make the onset figure
                if strcmpi(figType, 'on') || strcmpi(figType, 'both')
                    figure(onFig); clf;

                    % Locomotion Velocity
                    sp(Nsp) = subtightplot(Nsp,1,Nsp,opt{:}); %#ok<AGROW>
                    plot( periBout(runs).on.velocity(:,b), 'k' ); hold on; 
                    set(gca,'TickLength',[0.005,0], 'Xtick', onTicks, 'XtickLabel',TonTicks); % 
                    ylabel('Velocity (cm/s)'); xlabel('Peri-ON Time (s)');  % Time Running (s)
                    ylim( limStruct.velocity );

                    % Deformation variables
                    for s = 1:Nsp-1
                        sp(s) = MakeOnOffSubplot( periBout(runs), periStat(runs), showVars{s}, b, onTicks, [Nsp,1,s], opt, 'on' );
                    end
                    impixelinfo;
                    subplot(sp(1));  title( sprintf('Onset:  bout = %i', b) );
                    linkaxes(sp,'x');
                    axis tight;
                end
                % Make the offset figure
                if strcmpi(figType, 'off') || strcmpi(figType, 'both')
                    figure(offFig); clf;

                    % Locomotion Velocity
                    SP(Nsp) = subtightplot(Nsp,1,Nsp,opt{:}); %#ok<AGROW>
                    plot( periBout(runs).off.velocity(:,b), 'k' ); hold on; 
                    set(gca,'TickLength',[0.005,0], 'Xtick', offTicks, 'XtickLabel',ToffTicks); % 
                    ylabel('Velocity (cm/s)'); xlabel('Peri-OFF Time (s)');  % Time Running (s)
                    ylim( limStruct.velocity );

                    % Deformation variables
                    for s = 1:Nsp-1
                        SP(s) = MakeOnOffSubplot( periBout(runs), periStat(runs), showVars{s}, b, offTicks, [Nsp,1,s], opt, 'off' );
                    end
                    impixelinfo;
                    subplot(SP(1));  title( sprintf('Offset:  bout = %i', b) );
                    linkaxes(SP,'x');
                    axis tight;
                end
                pause;
            end
        end
    end
end

% BOUTWISE MEANS
if strcmpi(showType, 'mean')
    if strcmpi(figType, 'on') || strcmpi(figType, 'both')
        figure(onFig); clf;  
        % Locomotion Velocity
        sp(Nsp) = subtightplot(Nsp,1,Nsp,opt{:}); 
        onTemp = [periBout.on];
        onVelocity = [onTemp.velocity];
        plot( onVelocity, 'color',[0,0,0,0.2] ); hold on; 
        plot( mean(onVelocity, 2), 'color','k', 'LineWidth',1.5 );
        set(gca,'TickLength',[0.005,0], 'Xtick', onTicks, 'XtickLabel',TonTicks); % 
        ylabel('Velocity (cm/s)'); xlabel('Peri-ON Time (s)');  % Time Running (s)
        ylim( limStruct.velocity )
        % Deformation variables
        for s = 1:Nsp-1,  sp(s) = MakeOnOffMeanSubplot( periBout, showVars{s}, 'on', subInd, [Nsp,1,s], opt );  end   % , onTicks
        impixelinfo;
        subplot(sp(1));  title(sprintf('Onset:  mean of %i bouts (%s)', size(onVelocity,2), expt.name), 'Interpreter','none');
        linkaxes(sp,'x');
        axis tight;
        if ~isempty(figPath)
            export_fig( figPath, '-pdf', '-painters','-q101', '-append', onFig );
        else
            pause;
        end
    end
    if strcmpi(figType, 'off') || strcmpi(figType, 'both')
        figure(offFig); clf; 
        % Locomotion Velocity
        SP(Nsp) = subtightplot(Nsp,1,Nsp,opt{:}); 
        offTemp = [periBout.off];
        offVelocity = [offTemp.velocity];
        plot( offVelocity, 'color',[0,0,0,0.2] ); hold on; 
        plot( mean(offVelocity, 2), 'color','k', 'LineWidth',1.5 );
        set(gca,'TickLength',[0.005,0], 'Xtick', offTicks, 'XtickLabel',ToffTicks); % 
        ylabel('Velocity (cm/s)'); xlabel('Peri-OFF Time (s)');  % Time Running (s)
        ylim( limStruct.velocity );
        % Deformation variables
        for s = 1:Nsp-1,  SP(s) = MakeOnOffMeanSubplot( periBout, showVars{s}, 'off', subInd, [Nsp,1,s], opt );  end  % , offTicks
        impixelinfo;
        subplot(SP(1));  title( sprintf('Offset:  mean of %i bouts', size(offVelocity,2) ) );
        linkaxes(SP,'x');
        axis tight;
    end
end

end

function spOut = MakeOnOffSubplot( periBout, periStat, fieldName, b, scanTicks, spInd, opt, spType )
    if strcmpi(spType, 'on' )
        plotMat = periBout.on.(fieldName);
    elseif strcmpi(spType, 'off' )
        plotMat = periBout.off.(fieldName);
    end
    if ~isempty(b)
        plotMat = plotMat(:,:,b);
        tempStat = periStat.(fieldName); 
        if strcmpi(fieldName, 'fluor')
            sigEffectInd = find( tempStat.pEffect(b,:,3) < 0.05 & tempStat.effect(b,:,3) > 0); % which ROI showed significant activation?
        else
            sigEffectInd = find( tempStat.pEffect(b,:,3) < 0.05 ); % which planes showed sustained deformation?
        end
    else
        plotMat = mean(plotMat, 3, 'omitnan');
        sigEffectInd = [];
    end
    spOut = subtightplot(spInd(1), spInd(2), spInd(3), opt{:});
    imagesc( plotMat' ); 
    axPos = get(gca,'Position');
    CB(spInd(3)) = colorbar('EastOutside'); CB(spInd(3)).Label.String = fieldName; CB(spInd(3)).Label.FontWeight = 'bold'; CB(spInd(3)).Label.Interpreter = 'none';
    set(gca, 'Position',axPos, 'YTick',sigEffectInd, 'TickLength',[0.005,0], 'Xtick', scanTicks, 'XtickLabel',[], 'TickDir','out', 'YtickLabel',[] );
end

function spOut = MakeOnOffMeanSubplot( periBout, varName, setType, subInd, spInd, opt ) % , scanTicks
    if strcmpi(setType, 'on' )
        onTemp = [periBout.on];
        varMat = cat(3, onTemp.(varName));
    elseif strcmpi(setType, 'off' )
        offTemp = [periBout.off];
        varMat = cat(3, offTemp.(varName));
    end
    varMat = varMat - mean(varMat(subInd.(setType).pre,:,:), 1); % subtract each bout's pre-onset mean
    varMean = mean(varMat, 3, 'omitnan'); % average over bouts

    spOut = subtightplot(spInd(1), spInd(2), spInd(3), opt{:});
    imagesc( varMean' ); 
    axPos = get(gca,'Position');
    CB(spInd(3)) = colorbar('EastOutside'); CB(spInd(3)).Label.String = varName; CB(spInd(3)).Label.FontWeight = 'bold'; CB(spInd(3)).Label.Interpreter = 'none';
    set(gca, 'Position',axPos, 'TickLength',[0.005,0], 'TickDir','out', 'Xtick',[subInd.(setType).pre', subInd.(setType).short', subInd.(setType).long'],...
        'XtickLabel',[repmat({'Pre'}, 1, numel(subInd.(setType).pre) ), repmat({'Short'}, 1, numel(subInd.(setType).short) ), repmat({'Long'}, 1, numel(subInd.(setType).long) )]); 
end