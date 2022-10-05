function ViewROI3D(expt, ROI, fluor, loco, varargin)
%
IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'ROI', @isstruct )
addRequired( IP, 'fluor', @isstruct )
addRequired( IP, 'loco', @isstruct )
addParameter( IP, 'minInt', 2000, @isnumeric ) %
addParameter( IP, 'setROI', 1:expt.Nroi, @isnumeric)
addParameter( IP, 'save', false, @islogical)
parse( IP, expt, ROI, fluor, loco, varargin{:} );  %  exptDir, exptName
setROI = IP.Results.setROI;
%minInt = IP.Results.minInt;
saveToggle = IP.Results.save;
figPath = sprintf('%s%s_ROIsummary.pdf', expt.dir, expt.name );

% Reconcatenate fluor signals
Froi = [fluor.Froi]; Froi = vertcat(Froi.ROI);
Fnp = [fluor.Fnp]; Fnp = vertcat(Fnp.ROI);
F = [fluor.F]; F = vertcat(F.ROI);
Fo = [fluor.Fo]; Fo = vertcat(Fo.ROI);
dFF = [fluor.dFF]; dFF = vertcat(dFF.ROI);
activity = [fluor.act]; activity = vertcat(activity.ROI);
Fsub = cell(1,expt.Nroi); FsubNP = cell(1,expt.Nroi);
tempFroi = [fluor.Froi]; tempFnp = [fluor.Fnp];
%{
tempFsub = vertcat(tempFroi.sub); tempFsubNP = vertcat(tempFnp.sub);
for r = 1:expt.Nroi
    Fsub{r} = vertcat(tempFsub{:,r});
    FsubNP{r} = vertcat(tempFsubNP{:,r});
end
%}
speed = vertcat(loco.speedDown);

if saveToggle && exist(figPath,'file')
    delete(figPath); 
    fprintf('Writing %s', figPath ); 
end
opt = {[0.05,0.09], [0.06,0.04], [0.02,0.03]};  % {[vert, horz], [bottom, top], [left, right] }
close all; clearvars h sp;
figure('WindowState','maximized');  
sp(5) = subtightplot(6,2,12,opt{:});  plot(speed); title('Locomotion'); ylabel('Velocity (cm/s)'); xlabel('Scans per Run');  axis tight; set(gca,'Xtick',expt.scanLims, 'TickDir','out');
sp(4) = subtightplot(6,2,10,opt{:});  sp(3) = subtightplot(6,2,8,opt{:}); sp(2) = subtightplot(6,2,6,opt{:}); sp(1) = subtightplot(6,2,[2,4],opt{:}); % sp(5) = subtightplot(7,2,12,opt{:});
roiColor = [0,0,1];
for r = setROI 
    % Get labeled projections of subROI and convert to RGB 
    %{
    tempLabel = zeros(expt.Nrow, expt.Ny, expt.Nplane, 'uint8');
    tempLabel(ROI(r).ind) = 1; 
    tempRGBstack = zeros(expt.Nrow, expt.Ny, 3, expt.Nplane, 'uint8'); %clearvars tempRGBstack
    for z = flip(1:expt.Nplane)
        tempRGBstack(:,:,:,z) = label2rgb( tempLabel(:,:,z), roiColor);
    end
    tempRGBstack = permute( tempRGBstack, [1,2,4,3] );
    tempRGBproj = 255*ones(expt.Nrow+expt.Nplane, expt.Ny+expt.Nplane, 3, 'uint8' );
    tempRGBproj(1:expt.Nrow, 1:expt.Ny, :) = squeeze( min( tempRGBstack, [], 3) ); % z projection    
    if expt.Nplane > 1
        tempRGBproj(expt.Nrow+1:expt.Nrow+expt.Nplane, 1:expt.Ny, :) = permute( squeeze(min( tempRGBstack, [], 1)), [2,1,3]); % y projection
        tempRGBproj(1:expt.Nrow, expt.Ny+1:end, :) = squeeze(min( tempRGBstack, [], 2)); % x projection
    end
    %}
    % Get edges of 2d preliminary ROI
    %tempEdges = edge(ROI(r).mask_2d);
    %[tempEdgeRow, tempEdgeCol] = ind2sub( size(ROI(r).mask_2d), find( tempEdges ));
    % Show 2D preliminary ROI, and final ROI and subROI, including side projections
    subtightplot(6,2,[1,3,5], opt{:}); cla;
    %imshow( tempRGBproj, [] ); hold on;
    imshow( max(ROI(r).labelVol,[],3), []); hold on;
    plot( ROI(r).footprintEdge(:,2), ROI(r).footprintEdge(:,1), 'k.', 'MarkerSize',2 )
    %plot(tempEdgeCol, tempEdgeRow, 'k.', 'MarkerSize',2 ); 
    line([expt.Ncol, expt.Ncol],[1, expt.Nrow+expt.Nplane],'color','k');
    line([1, expt.Ncol+expt.Nplane],[expt.Nrow, expt.Nrow],'color','k'); 
    title( sprintf('%s:  ROI %i', expt.name, r ), 'Interpreter','none');  %  (vol = %i voxels,  footprint = %i pixels) , ROI(r).Nvox, ROI(r).footprint
    
    % Projection with ROI overlaid
    subtightplot(6,2,[7,9,11], opt{:}); cla;
    imshow( ROI(r).cropMax, [] ); hold on; axis image;
    %plot( ROI(r).footprintEdge(:,2), ROI(r).footprintEdge(:,1), 'b.', 'MarkerSize',2 );
    %xlim([ROI(r).crop(1), expt.Ncol-ROI(r).crop(2)]);
    %ylim([ROI(r).crop(3), expt.Nrow-ROI(r).crop(4)]);
    title('Max Projection');
    if ~saveToggle, impixelinfo; end

    % Raw signal - roi and np
    subplot(sp(1));   cla;
    plot( [Froi(:,r), Fnp(:,r)] ); hold on; ylabel('F'); title('Raw Signals'); legend('ROI','Neuropil'); axis tight; set(gca,'Xtick',expt.scanLims,'XtickLabel',[],'TickDir','out');
    % NP subtracted
    subplot(sp(2)); cla;
    plot(F(:,r)); hold on; plot( Fo(:,r), 'k' )  %plot([F(:,r), Fo(:,r)] ); 
    ylabel('F'); legend('F','Fo'); axis tight; set(gca,'Xtick',expt.scanLims,'XtickLabel',[],'TickDir','out'); title('Low-pass filtered, NP-Subtracted, Censored'); 
    % rescaled
    subplot(sp(3));   cla;
    plot(dFF(:,r)); title('High-Pass Filtered, Rescaled'); ylabel('dF/Fo'); axis tight; set(gca,'Xtick',expt.scanLims,'XtickLabel',[],'TickDir','out');
    % deconvolved
    subplot(sp(4));   cla;
    plot(activity(:,r)); title('Deconvolved'); ylabel('Activity');  % ,expt.scanLims, 'TickDir','out'
    linkaxes(sp,'x');
    set(gca,'Xtick',expt.scanLims,'XtickLabel',[],'TickDir','out');
    xlim([-Inf,Inf]); %axis tight; 
    % Save the figure (optional)
    if saveToggle 
        fprintf('\nROI %i / %i', r, numel(setROI) );
        pause(0.5);
        export_fig( figPath, '-pdf', '-painters','-q101', '-append', gcf );  
    else
        impixelinfo
        pause;
    end
end
close all;

end