clear; clc; close all;
dataDir = 'D:\2photon\';  %'C:\2photon';
dataSet = 'Afferents'; % 'Pollen'; %'Astrocyte'; %    'NGC'; % 'Neutrophil'; % 
% Parse data table
dataTablePath = 'R:\Levy Lab\2photon\ImagingDatasets.xlsx'; %'D:\MATLAB\ImagingDatasets.xlsx'; % 'D:\MATLAB\NGCdata.xlsx'; 
dataTable = readcell(dataTablePath, 'sheet',dataSet);  % 'NGC', ''
colNames = dataTable(1,:); dataTable(1,:) = [];
dataCol = struct('mouse',find(contains(colNames, 'Mouse')), 'date',find(contains(colNames, 'Date')), 'FOV',find(contains(colNames, 'FOV')), ...
    'volume',find(contains(colNames, 'Volume')), 'run',find(contains(colNames, 'Runs')), 'csd',find(contains(colNames, 'CSD')), 'ref',find(contains(colNames, 'Ref')), 'done',find(contains(colNames, 'Done')));
Nexpt = size(dataTable, 1);
dataTable(:,dataCol.date) = cellfun(@num2str, dataTable(:,dataCol.date), 'UniformOutput',false);

% Initialize variables
expt = cell(1,Nexpt); runInfo = cell(1,Nexpt); catInfo = cell(1,Nexpt); Tscan = cell(1,Nexpt); loco = cell(1,Nexpt);  fluor = cell(1,Nexpt); 
deform = cell(1,Nexpt); scaleCat = cell(1,Nexpt); ROIdeform = cell(1,Nexpt);
segParams = cell(1,Nexpt); ROI = cell(1,Nexpt); Nauto = nan(1,Nexpt); Nmanual = nan(1,Nexpt); axon = cell(1,Nexpt); stillCosSim = cell(1,Nexpt); stillCorr = cell(1,Nexpt); 
fluorEvent = cell(1,Nexpt); NfluorEvent = cell(1,Nexpt); fluorEventRate = cell(1,Nexpt); eventRaster = cell(1,Nexpt);
periBout = cell(1,Nexpt); periStat = cell(1,Nexpt); csdBout = cell(1,Nexpt); csdStat = cell(1,Nexpt); %periParam = cell(1,Nexpt); 
stillEpoch = cell(1,Nexpt); stillSumm = cell(1,Nexpt);
defVars = {'transAP', 'transML', 'transMag', 'scaleAP', 'scaleML', 'scaleMag', 'stretchAP', 'stretchML', 'shearAP', 'shearML', 'shearMag', 'shiftZ', 'DshiftZ'}; %
NdefVars = numel( defVars ); %, 'dShiftZ'
allVars = [{'fluor'}, defVars, {'velocity'}]; NallVars = numel(allVars);
%deformLim = struct('trans',[-Inf, Inf], 'scale',[-Inf, Inf], 'shear',[-Inf, Inf], 'shift',[-Inf, Inf]); %struct('trans',[-30, 30], 'scale',[0.97, 1.03], 'shear',[-0.02, 0.02], 'shift',[-3, 3]); %   
deformLim = struct('trans',[-20, 20], 'scale',[0.95, 1.05], 'stretch',100*[-1,1], 'shear',[-0.03, 0.03], 'shift',[-3.5, 3.5]);
viewLims = struct('trans',[-Inf,Inf], 'scale',[-Inf,Inf], 'stretch',[-Inf,Inf], 'shear',[-Inf,Inf], 'shift',[-3, 3], 'velocity',[-3, 15], 'fluor',[-1, 3]); 
GetRunNumber = @(x)(str2double(x(strfind(x, 'run')+3:end)));
% Various useful subsets of the data to choose from
xDone = find([dataTable{:,dataCol.done}] > 0);
x2D = [10, 12, 15, 21, 22, 23, 27, 28, 33, 34, 35, 36, 39, 43]; %intersect(xDone, find( [dataTable{:,dataCol.volume}] == 0 ));
x3D = [1, 3, 4, 5, 6, 8, 9, 11, 13, 14, 30, 31, 37, 38, 40, 41]; %intersect(xDone, find([dataTable{:,dataCol.volume}] == 1 ));
xCSD = [5, 8, 9, 11, 13, 14, 28, 30, 31, 36, 37, 41]; %intersect(find(cellfun(@mean, cellfun(@ismissing, dataTable(:,dataCol.csd), 'uniformoutput',false) ) ~= 1), xDone )'
xSham = [];
x3Dcsd = [5, 9, 13, 30, 37, 41]; % intersect( x3D, xCSD ); only include experiments where CSD wave was well-registered
x2Dcsd = intersect( x2D, xCSD );  
%xBase = setdiff(xDone, xCSD);
xPNS = 45:54;
xFiber = [1, 3, 5, 9, 10, 12, 14, 27, 28, 30, 31, 33, 35, 36, 37, 38, 39, 41, 43]; %37; % 11,13
xSlope = [5,9,30,31,37,41,28];
xSuper = [37, 41, 31];
xNGC = [62:65];
xCtrl = [66, 87, 91, 97, 98, 99]; % 80,   89,
% Choose which subset to  use
xPresent = xCSD; %x3D; % % x2D; 30; % %xCtrl; %x3Dcsd; %;5; %xCtrl; %x3Dcsd; %[x2D, xCtrl]; % xCtrl; %97; % x3D; %xCtrl; % 
Npresent = numel(xPresent);
overwriteROI = false;
for x = xPresent
    % Parse data table
    [expt{x}, runInfo{x}] = ParseDataTable(dataTable, x, dataCol, dataDir );
    expt{x}.Nroi = NaN; expt{x}.Naxon = NaN;
    if isnan(expt{x}.csd), expt{x}.preRuns = expt{x}.runs; else, expt{x}.preRuns = 1:expt{x}.csd-1; end
    % Get basic run-level metadata and fluorescence data
    for r = flip(expt{x}.runs)
        loco{x}(r) = GetLocoData( runInfo{x}(r), 'show',false ); % DuraDataPath(expt{x}.mouse, expt{x}.date, r, 'quad'), 
        fluor{x}(r) = GetFluor3D(runInfo{x}(r), 'overwrite',false);  % expt{x}, r
    end
    [Tscan{x}, runInfo{x}] = GetTime(runInfo{x}); 

    % Fill in experiment-level metadata
    expt{x}.Nrow = runInfo{x}(1).sz(1); expt{x}.Ncol = runInfo{x}(1).sz(2); 
    if all([runInfo{x}.otlevels]==runInfo{x}(1).otlevels) && all([runInfo{x}.nchan]==runInfo{x}(1).nchan)
        expt{x}.Nplane = runInfo{x}(1).otlevels;
        expt{x}.Nchan = runInfo{x}(1).nchan;
    else
        error('All runs must have the same number of planes and channels')
    end   
    expt{x}.Nscan = floor([runInfo{x}.nframes]/expt{x}.Nplane); expt{x}.totScan = sum(expt{x}.Nscan); expt{x}.totFrame = sum([runInfo{x}.totFrame]); 
    expt{x}.scanLims = [0, cumsum(expt{x}.Nscan)];
    expt{x}.frameRate = runInfo{x}(1).framerate; 
    expt{x}.scanRate = runInfo{x}(1).framerate/expt{x}.Nplane; 
    expt{x}.sbx.cat = strcat(expt{x}.dir, expt{x}.name, '.sbxcat '); %'.sbx_interp'
    catInfo{x} = ConcatenateRunInfo(expt{x}, runInfo{x}, 'suffix','sbxcat'); % Get concatenated metadata
    try
        expt{x}.zoom = str2double(catInfo{x}.config.magnification_list(catInfo{x}.config.magnification,:));
    catch
        fprintf('\nx = %i: %s - magnification is not clear, assumed to be 2.4', x, expt{x}.name)
        expt{x}.zoom = 2.4; 
    end
    expt{x}.umPerPixel = (1/0.53)/expt{x}.zoom; 
    expt{x}.refRun = floor(median(expt{x}.runs));
    expt{x}.sbx.cat = strcat(expt{x}.dir, expt{x}.name, '.sbxcat');
    expt{x}.sbx.reg = strcat(expt{x}.dir, expt{x}.name, '.sbxreg ');
    % {
    % Get deformation data
    [~, deform{x}, regParams, badInd] = GetDeformCat3D( catInfo{x}, deformLim, 'show',false, 'overwrite',false, 'window',find(Tscan{x}{1}<=32,1,'last') );  %  deformCat, affParams
    
    % Get locomotion state
    try
        loco{x} = GetLocoState(expt{x}, loco{x}, 'dir',strcat(dataDir, expt{x}.mouse,'\'), 'name',expt{x}.mouse, 'var','velocity', 'show',false); %
    catch
        fprintf('\n%s: GetLocoState failed!\n', expt{x}.name)
    end
    
    % Segment movies
    %{
    zProj = 1;    
    %badMat = false( expt{x}.Nscan, expt{x}.Nplane ); badMat(badInd) = true;  % size(deformCat{x}.scaleAP)
    censScans = []; % 1:30; %  find(any(badMat(:,zProj), 2))';
    segEdges = [80,80,40,40];
    SegmentCat3D(catInfo{x}, 'chan','green', 'zProj',zProj, 'censScans',censScans, 'chunkSize',30, 'overwrite',false, 'minFoot',50, 'edges',segEdges); % 'edges', segParams{x}.edges,  affParams(x).edges , 'censScans',segParams{x}.censScans segParams{x}.zProj
    %}
    segParams{x} = GetSegParams( catInfo{x} );

    % Load or generate and save mean and max projections
    maxProjPath = [expt{x}.dir, expt{x}.name, '_maxProj.tif']; %[expt{x}.dir, expt{x}.name, '_affineProj.tif'];
    meanProjPath = [expt{x}.dir, expt{x}.name, '_meanProj.tif'];
    if ~exist(maxProjPath, 'file')
        [~,zprojPath] = FileFinder(expt{x}.dir, 'contains','zproj', 'type','tif'); zprojPath = zprojPath{1};
        %zprojPath = sprintf('%s%s_zproj.tif',expt{x}.dir, expt{x}.name ); %double(pipe.io.read_tiff(zprojPath));
        fprintf('\n%s does not exist! Loading %s', maxProjPath, zprojPath); % expt{x}.sbx
        zprojMovie = imread_big(zprojPath); %readSBX(expt{x}.sbx, catInfo{x}, 100, expt{x}.scanLims(expt{x}.preRuns(end)+1)-100 );
        expt{x}.maxProj = max(zprojMovie, [] , 3);
        saveastiff(expt{x}.maxProj, maxProjPath);
        expt{x}.meanProj = mean(zprojMovie, 3);
        saveastiff(uint16(expt{x}.meanProj), meanProjPath);
        clearvars zprojMovie;
    else
        expt{x}.maxProj = loadtiff( maxProjPath );
        expt{x}.meanProj = loadtiff( meanProjPath );
    end

    % Load ROI
    try
        ROI{x} = MakeROI3D(expt{x}, 'overwrite',overwriteROI, 'corrPrct',90, 'minFoot',50);  % , 'corrPrct',90, [ROI{x}, preROI{x}] overwriteROI
        ROI{x} = rmfield(ROI{x}, 'correlation'); % correlation field is memory-intensive and not useful
        Nauto(x) = numel(ROI{x});
    catch
        Nauto(x) = 0;
    end
    %{
    try
        [ROI{x}, ~ ] = LoadManualROI(expt{x}, [3,3,1], [8,8,1; 21,21,2], segParams{x}.zProj ); %^LoadManualROI( expt{x}, 3 );
        Nmanual(x) = numel(ROI{x});
    catch
        Nmanual(x) = 0;
    end
    expt{x}.Nroi = numel(ROI{x});
    %}
   
    % Get fluor signals from ROIs
    if expt{x}.Nroi > 0
        % Visualize ROIs in FOV
        %[expt{x}.roiProj, ~, expt{x}.roiLabel] = VisualizeSegmentation(expt{x}, ROI{x}, 'overwrite',overwriteROI);   %  
        %ROI{x} = WriteROIproj(expt{x}, catInfo{x}, ROI{x}, 'edges',segParams{x}.edges, 'overwrite',overwriteROI); % ROI = , 'rSet',32  overwriteROI
        
        % Get ROI-specific fluor and deformation signals
        try
            fluor{x} = LoadROIfluor(expt{x}, fluor{x}, loco{x});
        catch
            fluor{x} = GetROIfluor(expt{x}, catInfo{x}, ROI{x}, fluor{x}, deform{x}, loco{x}, 'window',find(Tscan{x}{1}<=32,1,'last'), 'lp',0, 'deconvolve',false, 'overwrite',overwriteROI); %  
        end
        
        %if expt{x}.Nplane > 1,  ROIdeform{x} = GetROIdeform(ROI{x}, deform{x}, defVars);  end
        for run = flip(1:expt{x}.Nruns)
            [fluorEvent{x}{run}, NfluorEvent{x}(run,:), fluorEventRate{x}(run,:), eventRaster{x}{run}] = ...
                DetectEvents(Tscan{x}{run}, fluor{x}(run).z.ROI, fluor{x}(run).dFF.ROI, 'minMag',0.05, 'minDur',1, 'minZ',1, 'show',false); % expt{x}, 
            %[~, ~, loco{x}(run)] = PeriLoco3D(expt{x}, Tscan{x}{run}, loco{x}(run), deform{x}(run), fluor{x}(run).F.ROI, defVars, 'base',20, 'run',2, 'iso',[0,0]); % 'Nshuff',1000 periParam(x) periBout{x}(run)
        end
        
        if expt{x}.csd > 0
            [csdBout{x}(expt{x}.csd), csdStat{x}, csdParam(x)] = PeriCSD(expt{x}, Tscan{x}{expt{x}.csd}, loco{x}(expt{x}.csd), deform{x}(expt{x}.csd), fluor{x}(expt{x}.csd), defVars, 'show',false, 'overwrite',false);
        end
        %[stillEpoch{x}, stillSumm{x}] = GetStillEpochs(expt{x}, Tscan{x}, loco{x}, fluor{x}, deform{x}, csdBout{x}); % find long periods of stillness between locomotive bouts

        % Merge ROIs into putative axons and get axonal signals
        %[axon{x}, expt{x}, stillCosSim{x}, stillCorr{x}] = MergeROI3D(expt{x}, Tscan{x}, loco{x}, ROI{x}, fluor{x}, SAFE{x}, 'method','cluster', 'mergeThresh',1.5, 'show',false);  
        %fluor{x} = GetAxonFluor(expt{x}, catInfo{x}, axon{x}, fluor{x}, deform{x}, 'window',find(Tscan{x}{1}<=32,1,'last'), 'overwrite',false); 
        %fluor{x} = GetAxonFluor(expt{x}, catInfo{x}, fiber{x}, fluor{x}, deform{x}, 'window',find(Tscan{x}{1}<=32,1,'last'), 'overwrite',true); 
        %ROI{x} = WriteROIproj(expt{x}, catInfo{x}, ROI{x}, axon{x}, 'edges',segParam(x).edges, 'overwrite',false);
        %WriteROIproj(expt{x}, catInfo{x}, ROI{x}, fiber{x}, 'edges',segParams{x}.edges, 'overwrite',false);
    else
        fprintf('\n%s: No ROI found!', expt{x}.name)
    end
    %}
end
% May need to run BoutResponse, SpontVsLocoEvent, FiberAnalyis, CSDspread

%%
for x = xPresent
    % Merge ROIs into putative axons and get axonal signals
    [axon{x}, expt{x}, stillCosSim{x}, stillCorr{x}] = MergeROI3D(expt{x}, Tscan{x}, loco{x}, ROI{x}, fluor{x}, SAFE{x}, 'method','cluster', 'mergeThresh',2, 'show',true);    
    %fluor{x} = GetAxonFluor(expt{x}, catInfo{x}, axon{x}, fluor{x}, deform{x}, 'window',find(Tscan{x}{1}<=32,1,'last'), 'overwrite',false); 
    %fluor{x} = GetAxonFluor(expt{x}, catInfo{x}, fiber{x}, fluor{x}, deform{x}, 'window',find(Tscan{x}{1}<=32,1,'last'), 'overwrite',true); 
    %ROI{x} = WriteROIproj(expt{x}, catInfo{x}, ROI{x}, axon{x}, 'edges',segParam(x).edges, 'overwrite',false);
    %WriteROIproj(expt{x}, catInfo{x}, ROI{x}, fiber{x}, 'edges',segParams{x}.edges, 'overwrite',false);

end


%% check SNR
%{
SNR = cell(1,Nexpt);
for x = xPresent
    fluorCat = [fluor{x}(1:2).dFF];
    fluorCat = vertcat(fluorCat.ROI);
    for roi = flip(1:expt{x}.Nroi) %run = 1:expt{x}.Nruns
        dFFtemp = fluorCat(:,roi)';
        dFFtemp(isnan(dFFtemp)) = [];
        SNR{x}(roi) = pipe.proc.dffsnr(dFFtemp);
    end
    [~,SNRsort{x}] = sort(SNR{x}, 'descend');
    
    %plot( fluorCat(:,SNRsort{x}(end)) )
    
end


for x = xPresent
    plot( SNR{x}(boutResponse(x).neut), boutResponse(x).meanEffect(boutResponse(x).neut), 'k.' ); hold on;
    plot( SNR{x}(boutResponse(x).exc), boutResponse(x).meanEffect(boutResponse(x).exc), 'b.' ); 
    plot( SNR{x}(boutResponse(x).inh), boutResponse(x).meanEffect(boutResponse(x).inh), 'r.' ); hold off;
    xlabel('Signal-noise ratio'); ylabel('Mean Bout Effect (dF/Fstill)');
    title( sprintf('x = %i', x) );
    pause;
end
%}

%% Are measures of deformation comparable between 2D and 3D datasets?
sp(2) = subplot(2,1,2); sp(1) = subplot(2,1,1);
scaleLims = [0,5];
for x = xPresent
    scaleCat = deform{x}(1:2);
    scaleCat = vertcat(scaleCat.scaleMag);
    if expt{x}.Nplane > 1
        %scaleCat = max(scaleCat(:, segParams{x}.zProj), [], 2, 'omitnan');
        scaleCat = mean(scaleCat(:, segParams{x}.zProj), 2, 'omitnan'); %
        %scaleCat = scaleCat(:, segParams{x}.zProj);
        subplot(sp(2));
    else
        subplot(sp(1));
    end
    [fScale, xScale] = ecdf( scaleCat(:) );   
    plot(xScale, fScale); hold on;   
end
linkaxes(sp,'x');
xlim([0,5]);
subplot(sp(1)); title('2D'); line(scaleLims, 0.5*[1,1], 'color','k', 'linestyle','--');
subplot(sp(2)); title('3D'); line(scaleLims, 0.5*[1,1], 'color','k', 'linestyle','--');

%% Make a mean z-projection for each frame from concatenated data, pre-affine registration
catProjPath = sprintf('%s%s_cat_Zproj.tif', expt{x}.dir, expt{x}.name);
if ~exist(catProjPath,'file')
    zProj = 12:14;
    out = WriteSbxZproj(catInfo{x}.path, catInfo{x}, catProjPath, 1, expt{x}.totScan, 1, zProj);
end

