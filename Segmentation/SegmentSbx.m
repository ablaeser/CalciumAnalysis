function segParams = SegmentSbx(sbxInfo, segParams, varargin) % [zproj_interp, mask, trace_raw, trace_hp, segmentation_data, sbxPath] =
IP = inputParser;
addRequired( IP, 'sbxInfo', @isstruct )
addRequired( IP, 'segParams', @isstruct )
addParameter( IP, 'overwrite', false, @islogical )
parse( IP, sbxInfo, segParams, varargin{:} );
overwrite = IP.Results.overwrite;

%SETUP SEGMENTATION PARAMETERS
% Determine file paths
if sbxInfo.Nplane == 1
    sbxSuff = 'sbxreg'; % 'sbx_affine'; %sprintf('%s%s.sbx_interp', sbxInfo.dir, sbxInfo.exptName );
else
    sbxSuff = 'sbx_interp';
end
Nz = numel(segParams.z);
segParams.path.params = sprintf('%s%s_seg_params%s.mat', sbxInfo.dir, sbxInfo.exptName, segParams.name );
segParams.path.sbx = sprintf('%s%s.%s', sbxInfo.dir, sbxInfo.exptName, sbxSuff);
%segParams.path.proj = sprintf('%s%s_segProj.tif', sbxInfo.dir, sbxInfo.exptName);
for z = 1:Nz
    segParams.path.ica_raw{z} = sprintf('%s%s_2D_ica%s_z%i.mat', sbxInfo.dir, sbxInfo.exptName, segParams.name, z); %strcat(fdir,filesep,fname,'_2D_ica.mat');
    segParams.path.ica_filt{z} = sprintf('%s%s_2D_ica_filt%s_z%i.mat', sbxInfo.dir, sbxInfo.exptName, segParams.name, z);
    segParams.path.prelim{z} = sprintf('%s%s_seg_data_prelim%s_z%i.mat', sbxInfo.dir, sbxInfo.exptName, segParams.name, z);
    segParams.path.final{z} = sprintf('%s%s_seg_data_final%s_z%i.mat', sbxInfo.dir, sbxInfo.exptName, segParams.name, z);
end
% Limits of the data to be segmented
if isempty(segParams.z), segParams.z = 1:sbxInfo.Nplane; end
if isempty(segParams.edges)
    warning('Edges are undefined, defaulting to [80,80,40,40]');
    segParams.edges = [80,80,40,40];
end
hp_filt = @(x)(x - movmedian(x, segParams.hp_cutoff)); % high-pass filtering

% censor scans if needed
if isempty(segParams.seg_scan), segParams.seg_scan = 1:sbxInfo.totScan; end
if ~isempty(segParams.cens_scans)
    [~,cens_scan_ind] = intersect(segParams.seg_scan, segParams.cens_scans);
    segParams.seg_scan(cens_scan_ind) = [];
end
if sbxInfo.Nplane == 1, segParams.corr_thresh_pct = 0; end
%segParams.xyproj_width = 100; %width of stretched XY projection (pixels)
%segParams.min_vol = 100;

% make sbx_interp, if it's necessary and doesn't exist
if sbxInfo.Nplane > 1 && (~exist(segParams.path.sbx, 'file')  || overwrite )
    InterpTemporalCat( sbxInfo, 'edges',segParams.edges, 'pmt',-1, 'chunkSize',segParams.chunk_size );
end
% Generate/load projections to be used for segmentation
[seg_proj,~,seg_proj_path] = WriteSbxZproj(segParams.path.sbx, sbxInfo, 'z',segParams.z, 'chan','green', 'monochrome',true, 'name',[sbxInfo.exptName,'_segProj'], 'overwrite',overwrite);
segParams.path.proj = seg_proj_path;

% Save parameters used for segmentation
if ~exist(segParams.path.params,'file') || overwrite
    fprintf('\nSaving %s', segParams.path.params);
    save(segParams.path.params, 'segParams');
else
    fprintf('\nLoading %s', segParams.path.params);
    load(segParams.path.params, 'segParams');
end
% Run PCA/ICA on each projection
seg_data = cell(1,Nz);
for z = 1:Nz
    fprintf('\nz = %i of %i: Segmenting planes %i-%i\n', z, Nz, segParams.z{z}(1), segParams.z{z}(end))
    seg_proj{z} = seg_proj{z}(:,:,segParams.seg_scan); % restrict segmentation to seg_scans

    if exist(segParams.path.ica_filt{z}, 'file')
        fprintf('\nLoading %s... ', segParams.path.ica_filt{z});
        load( segParams.path.ica_filt{z} ); toc % load segmentation_data variable
    else
        % Run PCA/ICA on projection to generate initial 2D masks
        if ~exist(segParams.path.ica_raw{z}, 'file') || overwrite
            fprintf('\nPCA/ICA on projection');
            fix(clock)
            mask = PCAICA_2D_extract(seg_proj{z}, 'axons',true, 'hp_cutoff',segParams.hp_cutoff, 'edges',segParams.edges, 'save',segParams.path.ica_raw{z});   % trace_raw, trace_hp
            fix(clock)
        else
            fprintf('\nLoading %s... ', segParams.path.ica_raw{z});
            load(segParams.path.ica_raw{z});
        end

        % Filter and break up masks
        ROI = FilterMasks(mask, 'minFoot',segParams.min_foot, 'edge',segParams.edges, 'break',false, 'show',false); % false
        ROIlabel = zeros( size(mask{1}));
        for m = 1:numel(ROI)
            ROIlabel( ROI(m).mask ) = m;
        end
        imshow( label2rgb(ROIlabel) );
        mask = {ROI.mask}; Nmask = numel(mask);

        % Extract mask traces from full data and high-pass filter them
        trace_raw = cell(1,Nmask);
        fprintf('\nExtracting traces');
        w = waitbar(0,'Extracting PCA mask traces...'); % w =
        tic
        for t = 1:size(seg_proj{z},3)
            tempVol = seg_proj{z}(:,:,t);
            for m = 1:Nmask
                trace_raw{m}(1,t) = mean( tempVol(mask{m}(:)) );
            end
            waitbar(t/size(seg_proj{z},3));
        end
        delete(w);
        trace_hp = cellfun(hp_filt, trace_raw, 'UniformOutput',false);
        fprintf('\nSaving extracted traces to %s... ', segParams.path.ica_filt{z});
        save(segParams.path.ica_filt{z}, 'ROI', 'mask','Nmask', 'trace_raw', 'trace_hp', '-nocompression'); % , '-append'
    end

    % Preliminary 3D segmentation
    if ~exist(segParams.path.prelim{z},'file') || overwrite
        fprintf('\nExtracting 3D masks (n = %i masks)', numel(mask));
        seg_data{z} = MakeMask3D(segParams.path.sbx, sbxInfo, mask, trace_hp, trace_raw, segParams, segParams.z{z}, 'save',segParams.path.prelim{z}); %
    elseif ~exist(segParams.path.final{z}, 'file') || overwrite
        fprintf('\nLoading %s... ', segParams.path.prelim{z});
        load( segParams.path.prelim{z}, 'segmentation_data' ); toc % load segmentation_data variable
        seg_data{z} = segmentation_data;
    end
    toc
end
% Curate
for z = 1:Nz
    SegmentationROI3D( segParams.path.final{z}, seg_data{z});  %
end


%{
% Run PCA/ICA on 2D projection
tic
if ~exist(segParams.path.ica_filt, 'file') || overwrite
    % Generate a projection, if one doesn't exist
    if sbxInfo.Nplane > 1
        if exist(segParams.path.proj, 'file')
            zproj_interp = loadtiff(segParams.path.proj);
        else
            fprintf('\nProjecting planes %i - %i', segParams.z(1), segParams.z(end) );
            zproj_interp = WriteSbxZproj(segParams.path.sbx, sbxInfo, 'z',segParams.z, 'chan','green', 'monochrome',false); % segParams.seg_scan(1) numel(segParams.seg_scan) , 'firstScan',1 , 'type','zproj'
            zproj_interp = zproj_interp(:,:,segParams.seg_scan); % retain only good scans
            fprintf('\nWriting %s', segParams.path.proj); tic
            WriteTiff(uint16(zproj_interp), segParams.path.proj); toc
        end
    else
        if sbxInfo.totScan < 20000
            binT = 1;
        else
            binT = 10;
        end
        % NOTE - NEED TO DEAL WITH SUPPRESSED SCANS WHEN BINNING
        zproj_interp = WriteSbxPlaneTif(segParams.path.sbx, sbxInfo, 1, 'chan','green', 'binT',binT, 'firstScan',segParams.seg_scan(1), 'Nscan',numel(segParams.seg_scan), 'type','zproj', 'overwrite',true ); % , 'name',sbxInfo.exptName
    end

    % Run PCA/ICA on projection to generate initial 2D masks
    if ~exist(segParams.path.ica_raw, 'file') || overwrite
        fprintf('\nPCA/ICA on projection');
        fix(clock)
        mask = PCAICA_2D_extract(zproj_interp, 'axons',true, 'hp_cutoff',segParams.hp_cutoff, 'edges',segParams.edges, 'save',segParams.path.ica_raw);   % trace_raw, trace_hp
        fix(clock)
    else
        fprintf('\nLoading %s... ', segParams.path.ica_raw);
        load(segParams.path.ica_raw);
    end

    % Filter and break up masks
    ROI = FilterMasks(mask, 'minFoot',segParams.min_foot, 'edge',segParams.edges, 'break',false, 'show',false); % false
    ROIlabel = zeros( size(mask{1}));
    for m = 1:numel(ROI)
        ROIlabel( ROI(m).mask ) = m;
    end
    imshow( label2rgb(ROIlabel) );
    mask = {ROI.mask}; Nmask = numel(mask);

    % Extract mask traces from full data and high-pass filter them
    trace_raw = cell(1,Nmask);
    hp_filt = @(x)(x - movmedian(x, segParams.hp_cutoff));

    fprintf('\nExtracting traces');
    w = waitbar(0,'Extracting PCA mask traces...'); % w =
    tic
    for t = 1:size(zproj_interp,3)
        tempVol = zproj_interp(:,:,t);
        for m = 1:Nmask
            trace_raw{m}(1,t) = mean( tempVol(mask{m}(:)) );
        end
        waitbar(t/size(zproj_interp,3));
    end
    delete(w);
    %end
    trace_hp = cellfun(hp_filt, trace_raw, 'UniformOutput',false);
    fprintf('\nSaving extracted traces to %s... ', segParams.path.ica_filt);
    save(segParams.path.ica_filt, 'ROI', 'mask','Nmask', 'trace_raw', 'trace_hp', '-nocompression'); % , '-append'
else %if ~exist(segParams.path.prelim,'file') || overwrite
    fprintf('\nLoading %s... ', segParams.path.ica_filt);
    load(segParams.path.ica_filt);
    %FilterMasks(mask, 'minFoot',segParams.min_foot, 'edge',segParams.edges, 'break',false, 'show',true); % false
end
toc

tic
if ~exist(segParams.path.prelim,'file') || overwrite
    % Preliminary 3D segmentation
    fprintf('\nExtracting 3D masks (n = %i masks)', numel(mask));
    segmentation_data = Extract_3D_Mask(segParams.path.sbx, sbxInfo, mask, trace_hp, trace_raw, segParams.seg_scan, segParams.cens_scans, 'save',segParams.path.prelim, 'z',segParams.z, ...
        'blur',segParams.blur, 'prctile_cutoff',segParams.corr_thresh_pct, 'chunkSize',segParams.chunk_size); %
    SegmentationROI3D( segParams.path.final, segmentation_data);  %
elseif ~exist(segParams.path.final, 'file') || overwrite
    % Curate
    fprintf('\nLoading %s... ', segParams.path.prelim);
    load( segParams.path.prelim ); toc % load segmentation_data variable
    SegmentationROI3D( segParams.path.final, segmentation_data);  %
end
toc
%}
end