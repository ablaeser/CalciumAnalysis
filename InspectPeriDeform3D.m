function InspectPeriDeform3D(expt, periBout, defVars, varargin) % periStat, 
% Visualize deformation variables associated with locomotion bouts
IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'periBout', @isstruct )
%addRequired( IP, 'periStat', @isstruct )
addRequired( IP, 'defVars', @iscell )
addOptional( IP, 'limits', [], @isstruct)
addParameter(IP, 'bSet',1:periBout.Nbout, @isnumeric )
addParameter(IP, 'spType', 'all', @ischar )
parse( IP, expt, periBout, defVars, varargin{:} );  % periStat, 
limStruct = IP.Results.limits;
bSet =  IP.Results.bSet;
spType = IP.Results.spType;
if isempty(limStruct)
    limStruct.z = [-1,3];
    limStruct.trans = [-10, 10]; 
    limStruct.scale = [0.90, 1.1];
    limStruct.shear = [-0.05, 0.05];
    limStruct.shift = [-3, 3];
    limStruct.velocity = [-5,20];
end
limFields = fieldnames(limStruct);
Nvars = numel(defVars);
Nsp = Nvars+1; %10;
opt = {[0.01,0.07], [0.07,0.04], [0.2,0.2]};  % {[vert, horz], [bottom, top], [left, right] }
tickInt = round(2*expt.scanRate); %4;
if expt.Nplane == 1
    planeTicks = 1;
elseif expt.Nplane == 15
    planeTicks = [1:3:15, 15];
elseif expt.Nplane == 30
    planeTicks = [1:5:30, 30];
end
TL = [0.005,0];
close all; clearvars sp
figure('WindowState','maximized', 'Color','w', 'PaperOrientation','landscape');
for b = bSet 
    % Set xticks
    Ton = periBout.T{b} - periBout.Tstart(b); % periBout.T{b}(periBout.boutScan{b}(1));
    zeroTick = periBout.boutScan{b}(1); %find(periBout.T{b} - periBout.T{b}(periBout.boutScan{b}(1)) == 0);
    scanTicks = unique([flip(zeroTick:-tickInt:1), zeroTick:tickInt:numel(Ton)]);
    Tticks = arrayfun(@(x) sprintf('%3.1f',x), Ton(scanTicks),  'uniformoutput', false); 
   
    % Deformation variables
    for v = flip(1:Nvars)
        sp(v) = subtightplot(Nsp, 1 ,v, opt{:}); cla;
        [~,f] = intersect(limFields, defVars{v});
        if isempty(f)
            MakeDeformPlot( periBout.(defVars{v}){b}, defVars{v}, TL, scanTicks, planeTicks ); 
        else
            MakeDeformPlot( periBout.(defVars{v}){b}, defVars{v}, TL, scanTicks, planeTicks, limStruct.(limFields{f}) );
        end
        hold on;
        line(zeroTick*[1,1], [0, size(periBout.(defVars{v}){b},2)], 'color','k', 'lineStyle','--');
        line(periBout.boutScan{b}(end)*[1,1], [0, size(periBout.(defVars{v}){b},2)], 'color','k', 'lineStyle','--');
    end
    
    % Locomotion Velocity
    sp(Nsp) = subtightplot(Nsp,1,Nsp,opt{:}); %#ok<AGROW>
    plot( periBout.velocity{b}, 'k' ); hold on;
    %line([1,periBout.Nscan(b)], periStat.speed.pre(b)*[1,1], 'color','r', 'LineStyle','--' )
    %line([1,periBout.Nscan(b)], periStat.speed.run(b)*[1,1], 'color','b', 'LineStyle','--' )
    %line([1,periBout.Nscan(b)], periStat.speed.post(b)*[1,1], 'color','g', 'LineStyle','--' )
    h(1) = plot( periBout.preScan{b}, periBout.velocity{b}(periBout.preScan{b}), 'r.' );
    h(2) = plot( periBout.boutScan{b}, periBout.velocity{b}(periBout.boutScan{b}), 'b.' );
    if ~isempty(periBout.postScan{b}) 
        h(3) = plot( periBout.postScan{b}, periBout.velocity{b}(periBout.postScan{b}), 'g.' );  
        legend(h, 'Pre','Run','Post', 'Location','NorthEast');
    end
    hold off
    %text( 1.02, 0.5, sprintf('effect = %2.1f', periStat.speed.effect(b,3)), 'units','normalized', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'FontSize',9)
    set(gca,'TickLength',[0.005,0], 'Xtick', scanTicks, 'XtickLabel',Tticks); %
    ylabel('Velocity (cm/s)'); xlabel('Peri-bout Time (s)');  % Time Running (s)
    ylim( limStruct.velocity );
      
    if strcmpi(spType, 'all') 
        try impixelinfo; end
    end
    subplot(sp(1));  title( sprintf('bout = %i', b) );
    linkaxes(sp,'x');
    axis tight;
    pause;
end

end
%{
function spOut = MakeDeformSubplot( periBout, periStat, fieldName, b, scanTicks, spInd, opt, type )
    tempBout = periBout.(fieldName); 
    tempStat = periStat.(fieldName); 
    if strcmpi(fieldName, 'fluor')
        sigEffectInd = find( tempStat.pEffect(b,:,3) < 0.05 & tempStat.effect(b,:,3) > 0); % which ROI showed significant activation?
    else
        sigEffectInd = find( tempStat.pEffect(b,:,3) < 0.05 ); % which planes showed sustained deformation?
    end
    spOut = subtightplot(spInd(1), spInd(2), spInd(3), opt{:});
    if strcmpi(type, 'mean')
        tempMean = mean(tempBout{b}(:,sigEffectInd), 2, 'omitnan');
        plot( tempMean, 'k' ); hold on; 
        plot( periBout.preScan{b}, tempMean(periBout.preScan{b}), 'r.' )
        plot( periBout.boutScan{b}, tempMean(periBout.boutScan{b}), 'b.' )
        plot( periBout.postScan{b}, tempMean(periBout.postScan{b}), 'g.' )
        line([1,periBout.Nscan(b)], mean(tempStat.pre(b,sigEffectInd),'omitnan')*[1,1], 'color','r', 'LineStyle','--' ) 
        line([1,periBout.Nscan(b)], mean(tempStat.run(b,sigEffectInd),'omitnan')*[1,1], 'color','b', 'LineStyle','--' ) 
        line([1,periBout.Nscan(b)], mean(tempStat.post(b,sigEffectInd),'omitnan')*[1,1], 'color','g', 'LineStyle','--' ) 
        effectStr = sprintf('mean onset effect = %2.4f,  mean p = %2.4f\nmean offset effect = %2.4f,  mean p = %2.4f\nmean overall effect = %2.4f,  mean p = %2.4f', ...
            mean(tempStat.effect(b,sigEffectInd,1),'omitnan'), mean(tempStat.pEffect(b,sigEffectInd,1),'omitnan'), mean(tempStat.effect(b,sigEffectInd,2),'omitnan'),...
            mean(tempStat.pEffect(b,sigEffectInd,2),'omitnan'), mean(tempStat.effect(b,sigEffectInd,3),'omitnan'), mean(tempStat.pEffect(b,sigEffectInd,3),'omitnan'));
        text( 1.02, 0.5, effectStr, 'units','normalized', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'FontSize',9)
        ylabel(fieldName, 'Interpreter','none');
    elseif strcmpi(type, 'all')
        imagesc( tempBout{b}' ); hold on;
        axPos = get(gca,'Position');
        CB(spInd(3)) = colorbar('EastOutside'); CB(spInd(3)).Label.String = fieldName; CB(spInd(3)).Label.FontWeight = 'bold'; CB(spInd(3)).Label.Interpreter = 'none';
        set(gca, 'Position',axPos, 'YTick',sigEffectInd );
        if strcmpi(fieldName, 'fluor')
            ylabel('ROI');
        else
            ylabel('Plane');
        end
    else
        error('Invalid subplot type');
    end
    set(gca,'TickLength',[0.005,0], 'Xtick', scanTicks, 'XtickLabel',[], 'TickDir','out', 'YtickLabel',[]); %
    hold off;
end
%}