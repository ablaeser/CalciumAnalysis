%% View the big picture data
viewLims = struct('trans',[-Inf,Inf], 'scale',[-Inf,Inf], 'stretch',[-Inf,Inf], 'shear',[-0.03,0.03], 'shift',[-3, 3], 'velocity',[-3, 15], 'fluor',[-5, 5]); 
for x = [9,13] % xPrePost %xPresent % xPresent % x2D %xPresent %find(locoFrac(:,2) < 0.02)' %xWave %xPresent
    %sortROI = [boutResponse(x).pre.exc, boutResponse(x).pre.neut, boutResponse(x).pre.inh];
    %sortROI = [preCSD_summary{X}.rDeform, preCSD_summary{X}.rMixed, preCSD_summary{X}.rLoco]; 
    ViewResults3D( expt{x}, Tscan{x}, deform{x}, loco{x}, fluor{x}, allVars, ROI{x}, fiber{x}, 'cat',true, 'limits',viewLims); % , 'sortROI',sortROI

end

%% Inspect well isolated locomotion bouts
for x = 97 %xPresent % xWave %
    for run = expt{x}.preRuns %expt{x}.runs %expt{x}.csd %% 1:numel(periBoutCat{x})
        tempIso = vertcat(periBout{x}(run).iso);
        tempIsoPre = tempIso(:,1);
        isoBouts = find(tempIsoPre >= minIso)';
        if ~isempty(isoBouts)
            InspectPeriDeform3D( expt{x}, periBout{x}(run), {'fluor'}, 'spType','all', 'bSet',isoBouts ); % allVars , viewLims
        end
        %InspectPeriDeform3D( csdBout(x), periStat{x}(run), allVars, viewLims, 'spType','all' );
        %InspectPeriOnOff3D(expt{x}, periBout{x}, periStat{x}, allVars(1:end-1), viewLims, 'type','on', 'show','mean'); % pause(1);  % , 'path',strcat(figDir, 'BoutOnsetMean')
        %InspectPeri3D( csdBout(x), periParam{x}(expt{x}.csd), viewLims ); pause;
    end
end

%% Individual traces from inhibited ROI

preResp = [boutResponse(xPresent).pre];
[~,sortx] = sort([preResp.inhFrac], 'descend');
figure('WindowState','maximized');
for x = 43 %xPresent(sortx) %27 %xPresent
    Tpre = vertcat(Tscan{x}{1:2});
    Vpre = vertcat(loco{x}(1:2).Vdown);
    FroiPre = [fluor{x}(1:2).Froi]; FroiPre = vertcat(FroiPre.ROI);
    FnpPre = [fluor{x}(1:2).Fnp]; FnpPre = vertcat(FnpPre.ROI);
    FPre = [fluor{x}(1:2).F]; FPre = vertcat(FPre.ROI);
    FoPre = [fluor{x}(1:2).Fo]; FoPre = vertcat(FoPre.ROI);
    dFFpre = [fluor{x}(1:2).dFF]; dFFpre = vertcat(dFFpre.ROI);
    zpre = [fluor{x}(1:2).z]; zpre = vertcat(zpre.ROI);
    sp(1) = subplot(5,1,1); sp(2) = subplot(5,1,2); sp(3) = subplot(5,1,3); sp(4) = subplot(5,1,4); sp(5) = subplot(5,1,5);
    linkaxes(sp,'x');
    for roi = flip(boutResponse(x).pre.inh)
        subplot(sp(1)); cla;
        plot(Tpre, FroiPre(:,roi) ); hold on;
        plot(Tpre, FnpPre(:,roi) );
        title('Raw');
               
        subplot(sp(2)); cla;
        plot(Tpre, FPre(:,roi) ); hold on;
        plot(Tpre, FoPre(:,roi) );
        title('F, F0');
        
        subplot(sp(3)); cla;
        plot(Tpre, dFFpre(:,roi) ); hold on;
        title('dF/F0');
        
        subplot(sp(4)); cla;
        plot(Tpre, zpre(:,roi)); 
        title('z-score');
        
        subplot(sp(5)); cla;
        plot(Tpre, Vpre); 
        title('Velocity');
        xlim([-Inf,Inf]);
        pause;
    end
end


%%  View ROIs
for x = 31 %xPresent
    ROI{x} = WriteROIproj(expt{x}, catInfo(x), ROI{x}, 'edges',segParams{x}.edges, 'overwrite',true);
    
    [~,rSNR] = sort(fluor{x}(1).Froi.SNR, 'ascend');
    ViewROI3D(expt{x}, ROI{x}, fluor{x}, loco{x}, 'save',false, 'setROI',rSNR); % ,'minInt',2500 , 'setROI',1:5 stillSumm{x}.fluor.baseDown
    %{
    for f = 1:Nfiber(x)
        ViewROI3D(expt{x}, ROI{x}, fluor{x}, loco{x}, 'save',false, 'setROI',fiber{x}(f).ROI);
    end
    %}
end

%%
for x = 36
    MakeSbxPartial(expt(X), catInfo(X));  % , varargin
end

%% Max projections for 3D

for x = x3D
    zprojMovie = loadtiff( sprintf('%s%s_zproj.tif', expt{x}.dir, expt{x}.name ) );
    zprojMax = max(zprojMovie, [], 3);
    maxProjPath = sprintf('%s%s_maxProj.tif', expt{x}.dir, expt{x}.name );
    saveastiff(zprojMax, maxProjPath);
    fprintf('\nSaved %s', maxProjPath);
end


%% Segmentation and max proj side by side
SegmentationExample = figure('Units','normalized', 'OuterPosition', [0,0,1,1], 'Color','w', 'PaperOrientation','landscape');
opt = {[0.03,0.001], [0.07,0.04], [0.04,0.05]};  % {[vert, horz], [bottom, top], [left, right] }
sp(3) = subtightplot(1,3,3,opt{:}); sp(2) = subtightplot(1,3,2,opt{:}); sp(1) = subtightplot(1,3,1,opt{:}); 
for x = intersect( xPresent, x2D )   
    maxProj = expt{x}.maxProj(segParams{x}.edges(3)+1:end-segParams{x}.edges(4), segParams{x}.edges(1)+1:end-segParams{x}.edges(2));
    meanProj = expt{x}.meanProj(segParams{x}.edges(3)+1:end-segParams{x}.edges(4), segParams{x}.edges(1)+1:end-segParams{x}.edges(2));
    subplot(sp(1)); cla;
    imshow( maxProj, prctile(maxProj(:), [10,99.8]) );  
    hold on;
    %{
    if Nfiber(x) > 0
        fiberColor = distinguishable_colors(Nfiber(x));
        for f = 1:Nfiber(x)
            for roi = fiber{x}(f).ROI
                plot(ROI{x}(roi).footprintEdge(:,2) - segParams{x}.edges(1), ROI{x}(roi).footprintEdge(:,1) - segParams{x}.edges(3), '.', 'color',fiberColor(f,:), 'MarkerSize',1 )
            end
            %pause;
        end
    end
    %}
    %MakeScaleBar( round(expt{x}.umPerPixel*[40,40]), {[0,size(expt{x}.meanProj,2)]+0.5, [0,size(expt{x}.meanProj,1)]+0.5}, [0.05,0.82], [0,0], 'label',false, 'color','w', 'show',true );
    title( sprintf('x = %i: %s  Max Proj', x, expt{x}.name ), 'Interpreter','none');

    subplot(sp(2)); cla;
    imshow( meanProj, prctile(meanProj(:), [10,99.8]) );  
    hold on;
    title('Mean');
    
    subplot(sp(3)); cla;
    imshow( expt{x}.roiProj(segParams{x}.edges(3)+1:end-segParams{x}.edges(4), segParams{x}.edges(1)+1:end-segParams{x}.edges(2),:), [] ); hold on;
    for roi = 1:expt{x}.Nroi
        %text( ROI{x}(roi).cent(1) - segParams{x}.edges(1), ROI{x}(roi).cent(2) - segParams{x}.edges(3), num2str(roi), 'color','k', 'FontSize',8, 'HorizontalAlignment','center');
    end
    box on;
    linkaxes(sp,'xy');
    
    impixelinfo; pause;  
end

%% Write z-proj with ROIs overlaid
% Load zproj tif
movParam.displayPct = [1,99.5];
movParam.aviRate = 10; % frames per second
GrayOpt = struct('overwrite',true, 'message',false, 'append',false, 'big',true, 'color',false );
for x = 22 %xPresent
    movParam.binT = round(expt{x}.scanRate*30); %10;
    Tcat = vertcat( Tscan{x}{:} );
    if ~isnan(expt{x}.csd), Tcat = Tcat - csdBout{x}(expt{x}.csd).Tstart; end
    Tbin = BinDownMean( Tcat, movParam.binT );
    tempFrameFile = strcat(expt{x}.dir, 'tempFrame.jpg');
    tic
    % Generate a projection stack
    if expt{x}.Nplane == 1   
        ZprojStack = WriteSbxPlaneTif(expt{x}.sbx, catInfo(x), 1, 'firstScan',1, 'Nscan',expt{x}.totScan, 'binT',movParam.binT, 'verbose',true, ...
            'dir',expt{x}.dir, 'name',sprintf('%s_Bin%i',expt{x}.name, movParam.binT), 'edges',segParams{x}.edges);   % 
        NprojFrame = size(ZprojStack,3);
        %{
        NwriteFrame = NprojFrame; %
        % Overlay final ROI on downsampled data
        close all; clearvars storeFrames
        storeFrames = repmat(struct('cdata',[], 'colormap',[]), 1, NwriteFrame);
        aviFig = figure('color','k', 'Units','normalized', 'OuterPosition', [0 0 1 1]); % 'WindowState','maximized',
        preColor = distinguishable_colors(numel(preROI{x}));
        %w = waitbar(0, 'Writing concatenated movie');
        for z = flip(1:NwriteFrame) %
            cla;
            imshow(ZprojStack(:,:,z), displayLims, 'border','tight'); hold on;
            text(40, size(ZprojStack,1)-20, sprintf('%2.1f s', Tbin(z)), 'color','w' );
            for roi = 1:numel(preROI{x})
                plot( preROI{x}(roi).footprintEdge(:,2)-segParams{x}.edges(1)+1, preROI{x}(roi).footprintEdge(:,1)-segParams{x}.edges(3)+1, '.', 'color',preColor(roi,:), 'MarkerSize',1 );
                %text(600, 20, sprintf('ROI %02i',roi), 'HorizontalAlignment','center', 'color','w', 'FontSize',20 );
            end
            %impixelinfo
            exportgraphics(gca, tempFrameFile , 'Resolution',100); % save image with set resolution
            storeFrames(z) = im2frame(imread(tempFrameFile));  % getframe(aviFig);  % convert image to frame
            %waitbar(z/NwriteFrame, w)
            %pause;
        end
        close(aviFig);
        delete(tempFrameFile); %delete(w); 
        toc

        % Write AVI
        aviPath = sprintf('%s%s_preROImovie.avi', expt{x}.dir, expt{x}.name);
        fprintf('\nWriting %s...  ', aviPath); tic;
        writerObj = VideoWriter(aviPath);
        writerObj.FrameRate = movParam.aviRate; % set the seconds per image
        open(writerObj);
        for z = 1:NwriteFrame %numel(storeFrames)
            frame = storeFrames(z);    
            writeVideo(writerObj, frame);
        end
        close(writerObj);
        toc 

        % Overlay final ROI on downsampled data
        NwriteFrame = NprojFrame;
        close all; clearvars storeFrames
        storeFrames = repmat(struct('cdata',[], 'colormap',[]), 1, NwriteFrame);
        aviFig = figure('color','k', 'Units','normalized', 'OuterPosition', [0 0 1 1]); % 'WindowState','maximized',
        roiColor = distinguishable_colors(expt{x}.Nroi);
        %w = waitbar(0, 'Writing concatenated movie');
        for z = flip(1:NwriteFrame) %
            cla;
            imshow(ZprojStack(:,:,z), displayLims, 'border','tight'); hold on;
            text(40, size(ZprojStack,1)-20, sprintf('%2.1f s', Tbin(z)), 'color','w' );
            for roi = 1:expt{x}.Nroi
                plot( ROI{x}(roi).footprintEdge(:,2)-segParams{x}.edges(1)+1, ROI{x}(roi).footprintEdge(:,1)-segParams{x}.edges(3)+1, '.', 'color',roiColor(roi,:), 'MarkerSize',1 );
                %text(600, 20, sprintf('ROI %02i',roi), 'HorizontalAlignment','center', 'color','w', 'FontSize',20 );
            end
            %impixelinfo
            exportgraphics(gca, tempFrameFile , 'Resolution',100); % save image with set resolution
            storeFrames(z) = im2frame(imread(tempFrameFile));  % getframe(aviFig);  % convert image to frame
            %waitbar(z/NwriteFrame, w)
            %pause;
        end
        close(aviFig);
        delete(tempFrameFile); %delete(w); 
        toc

        % Write AVI
        aviPath = sprintf('%s%s_ROImovie.avi', expt{x}.dir, expt{x}.name);
        fprintf('\nWriting %s...  ', aviPath); tic;
        writerObj = VideoWriter(aviPath);
        writerObj.FrameRate = movParam.aviRate; % set the seconds per image
        open(writerObj);
        for z = 1:NwriteFrame %numel(storeFrames)
            frame = storeFrames(z);    
            writeVideo(writerObj, frame);
        end
        close(writerObj);
        toc 
        %}
    else
        tempTifName = sprintf('%s%s_bin%i.tif', expt{x}.dir, expt{x}.name, movParam.binT);
        if ~exist(tempTifName, 'file')
            tempCent = round(vertcat(ROI{x}.cent));
            ZprojStack = pipe.zproj(expt{x}.sbx, 1, expt{x}.totScan-1, 1, unique(tempCent(:,3)), 'mtype','.sbx_interp', 'registration',false);
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
    NwriteFrame = NprojFrame; %100;
    close all; clearvars storeFrames
    storeFrames = repmat(struct('cdata',[], 'colormap',[]), 1, NwriteFrame);
    aviFig = figure('color','k', 'Units','normalized', 'OuterPosition', [0 0 1 1]); % 'WindowState','maximized',
    roiColor = distinguishable_colors(expt{x}.Nroi);
    for z = flip(1:NwriteFrame) %
        cla;
        imshow(ZprojStack(:,:,z), displayLims, 'border','tight'); hold on;
        text(40, size(ZprojStack,1)-20, sprintf('%2.1f s', Tbin(z)), 'color','w' ); % 295
        for roi = 1:expt{x}.Nroi
            plot( ROI{x}(roi).footprintEdge(:,2)-segParams{x}.edges(1)+1, ROI{x}(roi).footprintEdge(:,1)-segParams{x}.edges(3)+1, '.', 'color',roiColor(roi,:), 'MarkerSize',1 ); %  0.5*[1,1,1]
            %text(600, 20, sprintf('ROI %02i',roi), 'HorizontalAlignment','center', 'color','w', 'FontSize',20 );
        end
        %impixelinfo
        exportgraphics(gca, tempFrameFile , 'Resolution',100); % save image with set resolution
        storeFrames(z) = im2frame(imread(tempFrameFile));  % getframe(aviFig);  % convert image to frame
    end
    close(aviFig);
    delete(tempFrameFile); 
    toc
    
    % Write AVI
    aviPath = sprintf('%s%s_ROImovie.avi', expt{x}.dir, expt{x}.name);
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



%% Generate movies of registered CSD bouts
movParam.dir = 'D:\MATLAB\Figures\CSD figures\'; mkdir(movParam.dir); % 'D:\MATLAB\LevyLab\Figures\CSD\Movies\';
movParam.fmtSpec = '%2.2f';
movParam.Tperi = [10,10]; % time before/after bout to show
movParam.zProj = [];
movParam.binT = 8;
movParam.displayPct = [5,99.9];
movParam.sbx = [];
movParam.regType = 'raw'; %'affine'; % % _axons
movParam.boutType = 'csd';
movParam.edges = []; % [80,80,20,20]; %[80,90,60,110]; %
movParam.aviRate = 10; % frames per second
movParam.Toffset = 0;
movParam.level = 'cat';  %'ind'; %'none'; %
for x = xCSD %xPresent 
    movParam.edges = segParams{x}.edges;
    movParam.zProj = segParams{x}.zProj; %segPlanes;
    if expt{x}.Nplane == 1
        movParam.sourceSbx = 'sbxz'; %  'sbx_affine'; %
    else
        movParam.sourceSbx = 'sbx_interp';
    end
    movParam.scalebar = [];
    %movParam.scalebar = MakeScaleBar( round(expt{x}.umPerPixel*[50,50]), {[0,expt{x}.Ncol-segParams{x}.edges(1)-segParams{x}.edges(2)]+0.5, [0,expt{x}.Nrow-segParams{x}.edges(3)-segParams{x}.edges(4)]+0.5},...
     %   [0.1,0.95], [0,0], 'label',false, 'color','w', 'show',false );
    [boutStack{x}, Tbout{x}, boutSpeed{x}, storeFrames{x}] = WriteBoutMovies(expt{x}, catInfo{x}, Tscan{x}, loco{x}, csdBout{x}, movParam); % , ROI{x} , 'overwrite',true (1)    , ROI{x}, axon{x}
end


%% Generate movies of fibers during registered CSD bouts
movParam.dir = 'D:\MATLAB\LevyLab\Figures\BoutMovies\CSD\Fibers\'; % 'D:\MATLAB\LevyLab\Figures\CSD\Movies\';
mkdir(movParam.dir)
movParam.fmtSpec = '%2.1f';
movParam.Tperi = [5,5]; %[15,10]; % time before/after bout to show
movParam.zProj = [];
movParam.binT = 1;
movParam.displayPct = [1,99.9];
movParam.sbx = [];
movParam.boutType = 'csd';
movParam.aviRate = 10; % frames per second
movParam.level = 'ind';
mkdir(movParam.dir)
for x = 9 %xPresent 
    movParam.edges = segParams{x}.edges;
    if expt{x}.Nplane > 1
        movParam.sourceSbx = 'sbx_interp';  
    else
        movParam.sourceSbx = 'sbx_affine'; %'sbxz'; %
    end
    for f = 1:Nfiber(x)
        movParam.zProj = unique(round(fiber{x}(f).cent(:,3)))'; %1; % %segParams{x}.zProj; %segPlanes;
        movParam.regType = sprintf('_aff_fiber%02.0f_', f); %'affine_axons'; 
        WriteBoutMovies(expt{x}, catInfo(x), Tscan{x}, loco{x}, csdBout{x}, movParam, ROI{x}, fiber{x}(f)); % , 'overwrite',true (1)  , ROI{x}  
    end
end

%% Generate movies of RAW CSD bouts
movParam.dir = 'D:\MATLAB\LevyLab\Figures\BoutMovies\CSD\'; % 'D:\MATLAB\LevyLab\Figures\CSD\Movies\';
movParam.fmtSpec = '%2.1f';
movParam.Tperi = [10,5]; % time before/after bout to show
movParam.zProj = [];
movParam.binT = 1;
movParam.displayPct = [1,99.5];
movParam.sourceSbx = []; 
movParam.regType = 'raw'; 
movParam.boutType = 'csd';
movParam.edges = []; % [80,80,20,20]; %[80,90,60,110]; %
movParam.aviRate = 10; % frames per second
movParam.level = 'ind';
mkdir(movParam.dir)
for x = 13 %xPresent 
    segParams(x) = GetSegParams( catInfo(x) );
    %[segEdges, segPlanes] = SegmentCat3D(catInfo(x));
    movParam.sbx = runInfo{x}(expt{x}.csd).path;
    movParam.edges = segParams(x).edges;
    movParam.zProj = segParams(x).zProj; %segPlanes;
    WriteBoutMovies(expt{x}, runInfo{x}(expt{x}.csd), Tscan{x}(expt{x}.csd), loco{x}(expt{x}.csd), csdBout{x}(expt{x}.csd), movParam, 'run',1); %  catInfo(x) , 'overwrite',true 
end

%% Generate movies of locomotion bouts

movParam.fmtSpec = '%2.2f';
movParam.Tperi = [10,10]; % time before/after bout to show
movParam.zProj = [];
movParam.binT = 1; % 15;
movParam.displayPct = [5,99.95];
movParam.sourceSbx = 'sbx_affine'; % 'sbx_partial'; %'sbxz'; % 
movParam.regType = 'affine'; % 'raw'; % 'partial'; % %
movParam.boutType = 'loco';
movParam.aviRate = 15; % frames per second
movParam.level = 'both'; % 'cat'; % 'ind'; % 
for x = 5 %xPresent
    movParam.dir = ['D:\MATLAB\LevyLab\Figures\BoutMovies\Locomotion\', expt{x}.name,'\']; % 'D:\MATLAB\LevyLab\Figures\CSD\Movies\';
    movParam.edges = segParams{x}.edges;
    movParam.zProj = 26:29; %segParams{x}.zProj;
    
    %movParam.scalebar = MakeScaleBar( round(expt{x}.umPerPixel*[50,50]), {[0,expt{x}.Ncol-segParams{x}.edges(1)-segParams{x}.edges(2)]+0.5, [0,expt{x}.Nrow-segParams{x}.edges(3)-segParams{x}.edges(4)]+0.5},...
    %    [0.1,0.95], [0,0], 'label',false, 'color','w', 'show',false );
    
    WriteBoutMovies(expt{x}, catInfo(x), Tscan{x}, loco{x}, periBout{x}, movParam); % , 'run',2  , ROI{x} , 'run',1, 'bout',5 , 'bout',3 , 'run',2
    %WriteBoutMovies(expt{x}, catInfo(x), Tscan{x}, periBout{x}, movParam, loco{x});
end

%% Generate movies of fibers during registered loco bouts
movParam.Tperi = [2,2]; %[15,10]; % time before/after bout to show
movParam.Toffset = 0;
movParam.zProj = [];
movParam.binT = 1;
movParam.displayPct = [4,99.9];
movParam.sbx = [];
movParam.fmtSpec = '%2.1f';
movParam.boutType = 'loco';
movParam.aviRate = 10; % frames per second
movParam.level = 'cat';
for x = xPresent 
    movParam.dir = strcat('D:\MATLAB\LevyLab\Figures\Fibers\', expt{x}.name, '\' );  mkdir(movParam.dir);
    movParam.edges = segParams{x}.edges;
    if expt{x}.Nplane > 1
        movParam.sourceSbx = 'sbx_interp';  
    else
        movParam.sourceSbx = 'sbx_affine'; %'sbxz'; %
    end
    movParam.scalebar = MakeScaleBar( round(expt{x}.umPerPixel*[50,50]), {[0,expt{x}.Ncol-segParams{x}.edges(1)-segParams{x}.edges(2)]+0.5, [0,expt{x}.Nrow-segParams{x}.edges(3)-segParams{x}.edges(4)]+0.5},...
        [0.1,0.95], [0,0], 'label',false, 'color','w', 'show',false );
    movParam.regType = sprintf('aff_'); %'affine_axons'; 
    %WriteBoutMovies(expt{x}, catInfo(x), Tscan{x}, loco{x}, periBout{x}, movParam, ROI{x});  % 
    
    for f = 8 %1%:Nfiber(x)
        if expt{x}.Nplane > 1
            movParam.zProj = unique(round(fiber{x}(f).cent(:,3)))'; %1; % %segParams{x}.zProj; %segPlanes;
        else
            movParam.zProj = 1;
        end
        movParam.regType = sprintf('aff_fiber%02.0f_', f); %'affine_axons'; 
        WriteBoutMovies(expt{x}, catInfo(x), Tscan{x}, loco{x}, periBout{x}, movParam, ROI{x}(fiber{x}(f).ROI), 'run',1 ); % , 'overwrite',true (1)  , ROI{x}   periBout{x}(1:expt{x}.csd-1)  , 'run',1:expt{x}.csd-1
    end
end



%% Generate plane movies
for x = xPresent
    tifDir = strcat( catInfo(x).dir, 'Ztifs\' ); %'D:\2photon\DL68\170712\Ztifs\';
    mkdir(tifDir);
    zWrite = [5,10,15,20,25,29]; %1:2:13; %[5,10,15,20,25,29]; % [16];
    %affPath = strcat(catInfo(x).dir,catInfo(x).exptName,'.sbx_affine'); % 'D:\2photon\DL68\170712\DL68_170712.sbx_affine';
    for z = zWrite
        %WriteSbxPlaneTif(catInfo(x).path, catInfo(x), z, 'dir',tifDir, 'name',catInfo(x).exptName, 'type','zinterp', 'edge',[0,0,0,0], 'verbose',true, 'bin',1, 'overwrite',true ); % , 'scale',2
        WriteSbxPlaneTif(expt{x}.sbx, catInfo(x), z, 'dir',tifDir, 'name',expt{x}.name, 'type','interp', 'edge',segParams{x}.edges, 'scale',1, 'verbose',true, 'chan',1, 'overwrite',false ); % [0,0,0,0]
        %WriteSbxPlaneTif(strcat(catInfo(x).dir,catInfo(x).exptName,'.sbx_affine'), catInfo(x), z, 'dir',tifDir, 'name',catInfo(x).exptName, 'type','aff', 'edge',[0,0,0,0], 'verbose',true, 'scale',4, 'overwrite',true ); % , 'scale',2

    end
end

%tifDir = 'D:\2photon\DL118\180628_DL118\180628_DL118_run3\Ztifs\';
%WriteSbxPlaneTif(runInfo{x}(expt{x}.csd).path, runInfo{x}(expt{x}.csd), z, 'dir',tifDir, 'name',sbxInfo.exptName, 'type','raw', 'edge',[0,0,0,0], 'scale',2, 'verbose',true, 'chan',1, 'overwrite',true );

%%
figure;
for x = intersect(xPresent, xCSD)
    plot( Tscan{x}{expt{x}.csd}-csdBout{x}(expt{x}.csd).Tstart(1), 100*fluor{x}(expt{x}.csd).dFF.vol ); hold on;
end
xlim([-Inf, 600]); ylim([-Inf,Inf]);
xlabel('Peri-CSD Time (s)'); ylabel('dF/F (whole volume, %)');
set(gca,'Xtick',[-60:60:600], 'TickDir','out', 'FontSize',24)