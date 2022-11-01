function segmentation_data = MakeMask3D(sbxPath, sbxInfo, mask, trace_hp, trace_raw, segParams, zProj, varargin)
IP = inputParser;
addRequired( IP, 'sbxPath', @ischar ) % path to the .sbx+ file to be projected
addRequired( IP, 'sbxInfo', @isstruct ) % @isstruct
addRequired( IP, 'mask', @iscell )
addRequired( IP, 'trace_hp', @iscell )
addRequired( IP, 'trace_raw', @iscell )
addRequired( IP, 'segParams', @isstruct )
addRequired( IP, 'z', @isnumeric )
addParameter( IP, 'save', '', @ischar ) 
addParameter( IP, 'overwrite', false, @islogical ) 
parse( IP, sbxPath, sbxInfo, mask, trace_hp, trace_raw, segParams, zProj, varargin{:} ); 
savePath = IP.Results.save;
tempPath = strcat(fileparts(sbxPath), '_corrTemp.mat');
filter_trace = @(x)(double(x - movmedian(x, segParams.hp_cutoff, 2)));
% extract each pixel/voxel's (HP-filtered) trace, and calculate correlation with the mask's (HP-filtered) mean trace
Nmask = numel(mask);
mask_Npix = cellfun(@sum, mask, repmat({[1,2]},1,Nmask)); 
c_vol = repmat( {zeros(sbxInfo.sz(1), sbxInfo.sz(2))}, Nmask, sbxInfo.Nplane ); % Nx
fix(clock)
tic
if sbxInfo.Nplane == 1
    if segParams.corr_thresh_pct > 0
        [chunkLims, Nchunk] = MakeChunkLims(projScans(1), projScans(end), projScans(end), 'size',segParams.chunkSize );
        trace_px_hp_chunk = cell(Nchunk, Nmask); % trace_px_chunk = cell(Nchunk,Nmask);
        w = waitbar(0,'Getting pixel traces...  ');
        % Extract pixel traces in chunks
        for c = 1:Nchunk %par
            fprintf('\nc = %i', c);
            chunkScans = chunkLims(c,1):chunkLims(c,2);
            NchunkScans = numel(chunkScans);
            [~,chunkCensScans,~] = intersect(chunkScans, censScans);
            NchunkCensScans = numel(chunkCensScans);
            plane_data = pipe.io.sbxRead(sbxPath, chunkScans(1), NchunkScans, 1); % pipe.io.sbxRead(path_interp,1,Nt,1,z);
            if NchunkCensScans > 0 %~isempty(censScans)
                fprintf('\nCensoring scans %i - %i', chunkScans(chunkCensScans(1)), chunkScans(chunkCensScans(end)) );
                plane_data(:,:,chunkCensScans) = []; % remove censored scans
                NchunkScans = NchunkScans - NchunkCensScans;
            end
            plane_data = imgaussfilt(plane_data, segParams.blur); %blur the loaded plane
            for m = 1:Nmask % par
                trace_px_chunk = plane_data(repmat(mask{m},[1,1,NchunkScans])); %mask the plane with the 2D mask {c,m}
                trace_px_chunk = reshape(trace_px_chunk,[],NchunkScans); % and reshape to Npx-by-t {c,m}
                trace_px_hp_chunk{c,m} = trace_px_chunk - movmedian(trace_px_chunk, segParams.hp_cutoff, 2); %high-pass filter the pixel traces
            end
            toc
            clearvars plane;
            waitbar(c/Nchunk,w)
        end
        delete(w);
        fprintf('\nSaving %s...', tempPath);
        save(tempPath, 'trace_px_hp_chunk', '-v7.3');  toc%  'c','trace_px_hp_chunk',
        %for each pixel in the plane, calculate the Pearson's correlation coefficient with the high-pass zproj trace
        fprintf('\nCalculating correlations...  ');
        for m = 1:Nmask %
            trace_px_hp = cat(2, trace_px_hp_chunk{:,m});
            Ctemp = zeros(1,sum(mask{m}(:)));
            for j = 1:size(trace_px_hp,1)
                Ctemp(j) = corr(double(trace_px_hp(j,:))', trace_hp{m}');
            end
            % reshape correlation into an image
            c_im = zeros(size(mask{m}));
            c_im(mask{m}) = Ctemp;
            c_vol{m,1} = c_im;
        end
        toc
    else
        fprintf('\nBypassing correlation extraction (corr threshold set to 0)');
        for m = 1:Nmask %
            c_im = zeros(size(mask{m}));
            c_im(mask{m}) = 1;
            c_vol{m} = c_im;
        end
    end
else
    fprintf('\nCalculating voxel-trace correlations')
    trace_length = numel(trace_hp{1});
    %mask_pix = cellfun(@repmat, mask, repmat({[1,1,trace_length]}, 1, Nmask), 'UniformOutput',false ); % too memory intensive
    tic
    for z = zProj %   par is faster, might be too memory intense
        if isempty(segParams.cens_scans) %#ok<PFBNS> 
            plane_data = readSBX(sbxPath, sbxInfo, segParams.seg_scan(1), trace_length, 1, z); 
        else
            plane_data = readSBX(sbxPath, sbxInfo, 1, sbxInfo.totScan, 1, z);  
            plane_data = plane_data(:,:,segParams.seg_scan);   % remove censored scans
        end
        plane_data = imgaussfilt(plane_data, segParams.blur);  %blur the loaded plane
        % Extract pixel-level traces (tried using logical indexing approach, but it was slower)
        mask_pix_trace = cell(1,Nmask);
        for scan = 1:trace_length
            tempScan = plane_data(:,:,scan);
            for m = 1:Nmask
                mask_pix_trace{m}(:,scan) = tempScan(mask{m});
            end
        end      
        % Calculate each pixel's correlation with the mask's trace
        %figure; subplot(1,2,1); imagesc(mask_pix_trace{1} ); subplot(1,2,2); imagesc(mask_pix_trace_hp{1} )
        mask_pix_trace_hp = cellfun(filter_trace, mask_pix_trace, 'UniformOutput',false);
        for m = 1:Nmask % par
            %for each pixel in the plane, calculate the Pearson's correlation coefficient with the high-pass zproj trace
            Ctemp = zeros(1,mask_Npix(m));
            for j = 1:mask_Npix(m)
                Ctemp(j) = corr(mask_pix_trace_hp{m}(j,:)', trace_hp{m}');
            end
            % reshape correlation into an image
            c_im = zeros(size(mask{m}));
            c_im(mask{m}) = Ctemp;
            c_vol{m,z} = c_im;
        end
    end
    toc
end

% Threshold the correlation volumes and make the 3D masks
fprintf('\nConstructing segmentation_data...');
segmentation_data = repmat( struct('mask_2d',[], 'mask_3d',[], 'mean_trace_raw',[], 'mean_trace_hp',[], 'correlation',[], 'xproj',[], 'yproj',[], 'zproj',[], 'keep',false), 1, Nmask );
parfor m = 1:Nmask %
    %reshape plane-by-plane correlations into volume
    cc = [c_vol{m,:}];
    cc = reshape(cc, sbxInfo.sz(1), sbxInfo.sz(2), []);
    %replace 0s with NaNs for projection
    c_nan = cc;
    c_nan(cc == 0) = NaN;
    %do X Y and Z projections
    c_nan_xproj = squeeze(nanmax(c_nan,[],1));
    c_nan_xproj = imresize(c_nan_xproj,[size(c_nan_xproj,1),segParams.xyproj_width]);
    c_nan_yproj = squeeze(nanmax(c_nan,[],2));
    c_nan_yproj = imresize(c_nan_yproj,[size(c_nan_yproj,1),segParams.xyproj_width]);
    c_nan_zproj = nanmax(c_nan,[],3);
    %threshold at segParams.corr_thresh_pct (typically 90th prctile)
    threshold = prctile(c_nan, segParams.corr_thresh_pct, 'all');
    c_thresh = cc > threshold;
    %label objects in binarized image
    c_L = bwlabeln(c_thresh);
    %calculate volume of each object
    rp = regionprops3(c_L,'Volume');
    %remove objects below volume cutoff
    F = find(rp.Volume > segParams.min_vol);
    c_filt = zeros(size(c_thresh)); %zeros(sbxInfo.sz(1),sbxInfo.sz(2),sbxInfo.Nplane); %
    for j = 1:numel(F)
        c_filt = c_filt | (c_L == F(j));
    end
    %figure, isosurface(c_filt);
    
    % save results to seg_data struct
    segmentation_data(m).mask_2d = mask{m};
    segmentation_data(m).mask_3d = c_filt;
    segmentation_data(m).mean_trace_raw = trace_raw{m};
    segmentation_data(m).mean_trace_hp = trace_hp{m};
    segmentation_data(m).correlation = c_nan;
    segmentation_data(m).xproj = c_nan_xproj;
    segmentation_data(m).yproj = c_nan_yproj;
    segmentation_data(m).zproj = c_nan_zproj;
    segmentation_data(m).keep = false;
end
fix(clock)

if ~isempty(savePath)
    fprintf('\nSaving %s...   ', savePath);
    save(savePath, 'segmentation_data', 'segParams', 'zProj', '-v7.3');
end
end