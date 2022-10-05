function ViewResults3D(expt, T, deform, loco, fluor, viewVars, varargin)%ROI,  projRGB = 
%
IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'T', @iscell )
addRequired( IP, 'deform', @isstruct )
addRequired( IP, 'loco', @isstruct )
addRequired( IP, 'fluor', @isstruct ) 
addRequired( IP, 'viewVars', @iscell )
addOptional( IP, 'ROI', [], @isstruct ) 
addOptional( IP, 'axon', [], @isstruct ) 
addParameter( IP, 'proj', '', @ischar )
addParameter( IP, 'limits', struct(), @isstruct)
addParameter( IP, 'fluorType', 'z', @ischar )
addParameter( IP, 'cat', true, @islogical )
addParameter( IP, 'sortROI', [], @isnumeric )
addParameter( IP, 'save', '', @ischar )
addParameter( IP, 'overwrite', false, @islogical )
parse( IP, expt, T, deform, loco, fluor, viewVars, varargin{:} );  % , ROI
ROI = IP.Results.ROI;
axon = IP.Results.axon;
projType = IP.Results.proj;
limStruct = IP.Results.limits;
figDir = IP.Results.save;
sortROI = IP.Results.sortROI;
%fluorType = IP.Results.fluorType;
catRuns = IP.Results.cat;
% Determine figure save path (optional)
if ~isempty(figDir) 
    if catRuns
        tifPath = strcat(figDir, expt.name, '_results_cat.tif'); 
        figPath = strcat(figDir, expt.name, '_results_cat.fig'); 
    else
        tifPath = strcat(figDir, expt.name, '_results_ind.tif'); 
        figPath = strcat(figDir, expt.name, '_results_ind.fig'); 
    end
    if exist(tifPath, 'file'), delete(tifPath); end
    if exist(figPath, 'file'), delete(figPath); end
else
    tifPath = ''; 
    figPath = '';
end

if isempty( fieldnames(limStruct) )
    limStruct = struct('fluor',[-Inf,Inf], 'trans',[-Inf,Inf], 'scale',[-Inf,Inf], 'shear',[-Inf,Inf], 'shift',[-3, 3], 'stretch',[-inf, inf]); 
    %limStruct = struct('dFF',[-3, 3], 'trans',[-10, 10], 'scale',[0.90, 1.1], 'shear',[-0.05, 0.05], 'shift',[-3, 3]); 
end
limFields = fieldnames(limStruct);
leftOpt = {[0.03,0.001], [0.07,0.04], [0.04,0.05]};  % {[vert, horz], [bottom, top], [left, right] }
rightOpt = {[0.01,0.001], [0.07,0.02], [0.04,0.05]};  % {[vert, horz], [bottom, top], [left, right] }
Nscan = cellfun(@numel, T);
runTicks = expt.scanLims; % cumsum(infoStruct.Nscan);
TL = [0.005,0];
if expt.Nplane == 30
    planeTicks = [1,5:5:30];
elseif expt.Nplane == 15
    planeTicks = [1:3:15,15];
elseif expt.Nplane == 1
    planeTicks = 1;
    viewVars( strcmpi(viewVars, 'shiftZ') ) = [];
    viewVars( strcmpi(viewVars,'dShiftZ') ) = [];
else
    planeTicks = round(linspace(1, expt.Nplane, 5));
end
Nview = numel(viewVars);
close all; clearvars sp SP;
ResultsSummFig = figure('Units','normalized', 'OuterPosition', [0,0,1,1], 'Color','w', 'PaperOrientation','landscape'); 
% ROI/Axon masks
if ~isempty(ROI)
    Nroi = numel(ROI);
    labelMat = zeros(expt.Nrow, expt.Ncol, expt.Nplane, 'uint8' );
    if ~isempty(axon)
        Naxon = numel(axon);
        segColor = distinguishable_colors(Naxon);
        for a = 1:Naxon
            for r = axon(a).ROI %1:expt.Nroi
                labelMat( ROI(r).ind ) = a;
            end
            labelTitle = sprintf('%s: Axons', expt.name);
        end
    else
        segColor = distinguishable_colors(Nroi);
        for r = 1:Nroi
            labelMat( ROI(r).ind ) = r;
            %text( )
        end
        labelTitle = sprintf('%s: ROI', expt.name);
    end
    
    projPath = sprintf('%s%s_%sProj.tif', expt.dir, expt.name, projType);
    SP(1) = subtightplot(2,2,1, leftOpt{:});
    if ~isempty(projType) && exist(projPath, 'file')
        projImage = loadtiff(projPath);
        if size(projImage,3) > 1, projImage = max(projImage,[],3); end
        displayPct = [10,99.7];
        displayLims = [prctile(projImage(:),displayPct(1)), prctile(projImage(:),displayPct(2))];
        imshow(projImage, displayLims); hold on;
        title( labelTitle, 'Interpreter','none' ); 
        if ~isempty(axon)
            for a = 1:expt.Naxon
                for r = axon(a).ROI
                    text( ROI(r).cent(1), ROI(r).cent(2), sprintf('%i-%i', r, a), 'HorizontalAlignment','center', 'FontSize',8 );
                end
            end
        else
            for r = 1:Nroi
                plot( ROI(r).footprintEdge(:,2), ROI(r).footprintEdge(:,1), '.', 'color', segColor(r,:), 'MarkerSize',1 )
                text( ROI(r).cent(1), ROI(r).cent(2), sprintf('%i', r), 'HorizontalAlignment','center', 'FontSize',8 );
            end
        end
    else
        projRGB = label2rgb( max( labelMat,[] ,3) );
        imshow( projRGB ); hold on;
        title( labelTitle, 'Interpreter','none' ); 
        if ~isempty(axon)
            for a = 1:expt.Naxon
                for r = axon(a).ROI
                    text( ROI(r).cent(1), ROI(r).cent(2), sprintf('%i-%i', r, a), 'HorizontalAlignment','center', 'FontSize',8 );
                end
            end
        % {
        else
            for r = 1:Nroi
                %plot( ROI(r).footprintEdge(:,2), ROI(r).footprintEdge(:,1), '.', 'color', segColor(r,:), 'MarkerSize',1 )
                text( ROI(r).cent(1), ROI(r).cent(2), sprintf('%i', r), 'HorizontalAlignment','center', 'FontSize',8 );
            end
        %}
        end
    end
end

SP(3) = subtightplot(2,2,3, leftOpt{:}); %subtightplot(9,2,11:2:18, opt{:});
imshow( expt.maxProj, [] ); hold on;
colormap(gca, gray)
linkaxes(SP,'xy');
impixelinfo;
if catRuns 
    tempFluor = [fluor.z]; % [fluor.dFF]; %
    if isfield(tempFluor, 'ROI')
        tempFluorCat = vertcat( tempFluor.ROI ); 
        if ~isempty(sortROI), tempFluorCat = tempFluorCat(:,sortROI); end
        fluorUnit = 'ROI';
    else
        tempFluorCat = vertcat( tempFluor.plane ); 
        fluorUnit = 'Plane';
    end
    %{
    % Correlation
    tempCorr = corr( vertcat( tempFluorCat ), 'Rows','complete' );
    tempCorr(tempCorr==1) = NaN;
    subtightplot(2,2,3, leftOpt{:}); %subtightplot(9,2,11:2:18, opt{:});
    imagesc(tempCorr); axis square; 
    title('Correlation (dF/F)'); xlabel(fluorUnit); ylabel(fluorUnit);
    colorbar;
    %}
    % Plot concatenated variables
    Tint = 180;
    [~,XtickInt] = min( abs(T{1} - Tint) );
    Ttick = sprintfc('%2.1f', T{1}(1:XtickInt:Nscan(1))/60 );
    for v = 1:Nview
        sp(v) = subtightplot(Nview,2,2*v,rightOpt{:}); cla;
        if strcmpi(viewVars{v}, 'speed') || strcmpi(viewVars{v}, 'velocity')
            MakeDeformPlot( vertcat(loco.Vdown), viewVars{v}, TL, runTicks, planeTicks, limStruct.velocity ); 
            hold on;
            if isfield(loco,'bout')
                plot( -1*vertcat(loco.bout), 'color',[0.6,0.6,0.6], 'LineWidth',2 );  % [0.5,0.5,0.5, 0.5]  %'k'
                tempBout = vertcat(loco.bout); tempBout(isnan(tempBout)) = 0;
                tempCC = bwconncomp(vertcat(tempBout));
                for b = 1:tempCC.NumObjects
                    plot( tempCC.PixelIdxList{b}(1), -1, 'ko');
                    plot( tempCC.PixelIdxList{b}(end), -1, 'kx');
                end
            end
            % {
            if isfield(loco,'stateDown')
                stateVec = vertcat(loco.stateDown);
                stateVec(stateVec == 1) = NaN;
                plot( -1.2*stateVec, 'r*', 'markerSize',1 ); 
                plot( -1.2*stateVec, 'color','r', 'LineWidth',2 ); 
            end
            %}
            ylabel('Velocity (cm/s)'); xlabel('Time (min)'); 
            set(gca, 'Xtick',1:XtickInt:Nscan(1), 'XTickLabel',Ttick, 'TickDir','out', 'TickLength',TL);
            xtickangle(30);
        elseif strcmpi(viewVars{v}, 'fluor')
            MakeDeformPlot( tempFluorCat, viewVars{v}, TL, runTicks, unique([1:round(Nroi/5):Nroi, Nroi]), limStruct.fluor ) % planeTicks unique([1:10:Nroi,Nroi])
            ylabel('ROI'); title( expt.name, 'Interpreter','none' )
        else
            for f = 1:numel(limFields)
                if contains( viewVars{v}, limFields{f}, 'ignoreCase',true )
                    MakeDeformPlot( vertcat( deform.(viewVars{v}) ), viewVars{v}, TL, runTicks, planeTicks, limStruct.(limFields{f}) ) 
                end
            end       
        end
        %{
        if strcmpi(viewVars{v}, 'speed')
            % Wheel velocity
            Vcat = vertcat(loco.Vdown);
            boutCat = vertcat(loco.bout);
            plot( Vcat ); hold on; % 'Color',[0,0,0,AlphaVal] periSpeed
            plot( -1.5*boutCat, 'Color','k', 'LineWidth',1.5 ); % [0,0,0,0.5]
            ylim(limStruct.velocity);
            ylabel('Velocity (cm/s)'); xlabel('Time (min)'); 
            set(gca, 'Xtick',1:XtickInt:Nscan(1), 'XTickLabel',Ttick, 'TickDir','out', 'TickLength',[0.005,0]);
            xtickangle(30);
        elseif strcmpi(viewVars{v}, 'fluor')
            imagesc( tempFluorCat' );
            axPos = get(gca,'Position');
            CB(v) = colorbar('EastOutside'); CB(v).Label.String = sprintf('%s (%s)',viewVars{v}, fluorType); CB(v).Label.FontWeight = 'bold'; CB(v).Label.Interpreter = 'none'; 
            set(gca,'Position',axPos, 'TickLength',[0.005,0],'Ytick',planeTicks, 'YTickLabel',[], 'TickDir','out', 'Xtick',[1,cumScan(1:end-1)], 'XTickLabel',1:expt.Nruns); % , 'Xtick',Xticks, 'XTickLabel',[]
            ylabel('ROI');
            caxis(limStruct.fluor); 
        else
            defCat = vertcat( deform.(viewVars{v}) );
            imagesc( defCat' );
            axPos = get(gca,'Position');
            CB(v) = colorbar('EastOutside'); CB(v).Label.String = viewVars{v}; CB(v).Label.FontWeight = 'bold'; CB(v).Label.Interpreter = 'none'; 
            set(gca,'Position',axPos, 'TickLength',[0.005,0],'Ytick',planeTicks, 'YTickLabel',[], 'TickDir','out', 'Xtick',[1,cumScan(1:end-1)], 'XTickLabel',[]); % , 'Xtick',Xticks, 'XTickLabel',[]
            ylabel('Plane'); 
            % Is there a prescribed color range?
            for f = 1:numel(limFields)
                if contains( viewVars{v}, limFields{f}, 'ignoreCase',true )
                    caxis(limStruct.(limFields{f})); 
                    break
                end
            end
        end
        %}
        %yRange = get(gca,'Ylim');
        %for r = 1:expt.Nruns-1, line(cumScan(r)*[1,1], yRange, 'color','k', 'LineWidth',1); end
    end
    linkaxes(sp,'x');
    xlim([-Inf, Inf]);  %axis tight;
    % Save the results to a tiff
    if ~isempty(tifPath)
        tic
        fprintf('\nSaving %s...', tifPath);
        print(ResultsSummFig, tifPath, '-dtiff' );  %export_fig( figPath, '-pdf', '-painters','-q101', '-append', gcf );
        toc
        pause(1);
        fprintf('\nSaving %s...', figPath);
        savefig(figPath);
        toc
    else
        impixelinfo;
        pause;
    end
    %close all;
else
    for r = expt.runs % find(~cellfun(@isempty, T)) %1:numel(infoStruct)
        Xticks = unique([1:tickInt:Nscan(r),Nscan(r)]); 
        Tticks =  num2str(T{r}(Xticks), '%2.1f'); % convert scan ticks to time units
        % Intrarun correlation between fluor signals
        subtightplot(2,2,3, rightOpt{:}); %subtightplot(9,2,11:2:18, opt{:});
        tempCorr = fluor(r).dFF.corr; tempCorr(tempCorr==1) = NaN;
        imagesc(tempCorr); axis square; 
        title('Correlation (dF/F)'); xlabel('ROI'); ylabel('ROI');
        colorbar;
        % Display specified variables
        for v = 1:Nview
            sp(v) = subtightplot(Nview,2,2*v,rightOpt{:}); cla;
            MakeDeformPlot( deformData, deformName, TL, Xticks, Yticks, deformationLims )
            %{
            if strcmpi(viewVars{v}, 'speed')
                % Wheel velocity
                plot( loco(r).Vdown ); hold on; % 'Color',[0,0,0,AlphaVal] periSpeed
                plot( -1*loco(r).bout, 'Color',[0,0,0,0.5], 'LineWidth',1.5 );
                ylim(limStruct.velocity);
                set(gca, 'Xtick',Xticks, 'XTickLabel',Tticks, 'TickDir','out', 'TickLength',[0.005,0]);
                xtickangle(30);
                ylabel('Velocity (cm/s)'); xlabel('Time (s)'); 
            elseif strcmpi(viewVars{v}, 'fluor')
                % fluor signal
                imagesc( fluor(r).(fluorType).ROI' );
                axPos = get(gca,'Position');
                CB(v) = colorbar('EastOutside'); CB(v).Label.String = sprintf('%s (%s)',viewVars{v}, fluorType); CB(v).Label.FontWeight = 'bold'; CB(v).Label.Interpreter = 'none'; 
                set(gca,'Position',axPos, 'TickLength',[0.005,0], 'Xtick',Xticks, 'XTickLabel',[],'Ytick',planeTicks, 'YTickLabel',[], 'TickDir','out')
                ylabel('Plane');
                caxis(limStruct.fluor); 
            else
                imagesc( deform(r).(viewVars{v})' );
                axPos = get(gca,'Position');
                CB(v) = colorbar('EastOutside'); CB(v).Label.String = viewVars{v}; CB(v).Label.FontWeight = 'bold'; CB(v).Label.Interpreter = 'none'; 
                set(gca,'Position',axPos, 'TickLength',[0.005,0], 'Xtick',Xticks, 'XTickLabel',[],'Ytick',planeTicks, 'YTickLabel',[], 'TickDir','out')
                ylabel('Plane'); 
                % Is there a prescribed color range?
                for f = 1:numel(limFields)
                    if contains( viewVars{v}, limFields{f}, 'ignoreCase',true )
                        caxis(limStruct.(limFields{f})); 
                        break
                    end
                end
            end
            %}
        end
        linkaxes(sp,'x');
        axis tight;

        % Save the results to a tiff
        if ~isempty(tifPath)
            tic
            fprintf('\nSaving %s...', tifPath);
            print(ResultsSummFig, tifPath, '-dtiff' );  %export_fig( figPath, '-pdf', '-painters','-q101', '-append', gcf );
            toc
            pause(1);
            fprintf('\nSaving %s...', figPath);
            savefig(figPath, ResultsSummFig);
            toc
        else
            impixelinfo;
            %pause;
        end
    end
    %close all;
end
end