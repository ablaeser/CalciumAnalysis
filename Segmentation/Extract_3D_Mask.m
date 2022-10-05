function segmentation_data = Extract_3D_Mask(sbxPath, sbxInfo, mask, trace_hp, trace_raw, projScans, censScans, varargin)
p = inputParser;
addOptional(p,'prctile_cutoff',90); %correlation coefficient cutoff prctile
addOptional(p,'volume_cutoff',100); %size filter for objects in 3d mask (voxels)
addOptional(p,'xyproj_width',100); %width of stretched XY projection (pixels)
addOptional(p,'blur',1); %gaussian blur width for trace extraction
addOptional(p,'hp_cutoff',20); %high pass filter cutoff for traces
addOptional(p,'chunkSize',1000)
addParameter(p, 'z', 1:sbxInfo.Nplane )
addParameter(p,'save','',@ischar)
parse(p,varargin{:});
p = p.Results;
tempPath = strcat(fileparts(sbxPath), '_corrTemp.mat');
Nmask = numel(mask);

% extract each pixel/voxel's (HP-filtered) trace, and calculate correlation with the mask's (HP-filtered) mean trace
c_vol = repmat( {zeros(sbxInfo.sz(1), sbxInfo.sz(2))}, Nmask, sbxInfo.Nplane ); % Nx
fix(clock)
tic
if sbxInfo.Nplane == 1
    if p.prctile_cutoff > 0
        [chunkLims, Nchunk] = MakeChunkLims(projScans(1), projScans(end), projScans(end), 'size',p.chunkSize );
        trace_px_hp_chunk = cell(Nchunk, Nmask); % trace_px_chunk = cell(Nchunk,Nmask);
        w = waitbar(0,'Getting pixel traces...  ');
        % Extract pixel traces in chunks
        for c = 1:Nchunk %par
            fprintf('\nc = %i', c);
            chunkScans = chunkLims(c,1):chunkLims(c,2);
            NchunkScans = numel(chunkScans);
            [~,chunkCensScans,~] = intersect(chunkScans, censScans);
            NchunkCensScans = numel(chunkCensScans);
            plane_chunk = pipe.io.sbxRead(sbxPath, chunkScans(1), NchunkScans, 1); % pipe.io.sbxRead(path_interp,1,Nt,1,z);
            if NchunkCensScans > 0 %~isempty(censScans)
                fprintf('\nCensoring scans %i - %i', chunkScans(chunkCensScans(1)), chunkScans(chunkCensScans(end)) );
                plane_chunk(:,:,chunkCensScans) = []; % remove censored scans
                NchunkScans = NchunkScans - NchunkCensScans;
            end
            plane_chunk = imgaussfilt(plane_chunk, p.blur); %blur the loaded plane
            for m = 1:Nmask % par
                trace_px_chunk = plane_chunk(repmat(mask{m},[1,1,NchunkScans])); %mask the plane with the 2D mask {c,m}
                trace_px_chunk = reshape(trace_px_chunk,[],NchunkScans); % and reshape to Npx-by-t {c,m}
                trace_px_hp_chunk{c,m} = trace_px_chunk - movmedian(trace_px_chunk, p.hp_cutoff, 2); %high-pass filter the pixel traces
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
    [~,projCensScans] = intersect(projScans, censScans );
    trace_length = numel(trace_hp{1});
    parfor z = p.z %  
        plane_chunk = readSBX(sbxPath, sbxInfo, projScans(1), numel(projScans), 1, z); %pipe.io.sbxRead(sbxPath,projScans(1),numel(projScans)-1,1,z); % pipe.io.sbxRead(path_interp,1,Nt,1,z);
        plane_chunk(:,:,projCensScans) = []; % remove censored scans
        % make sure the length of plane_chunk matches trace_hp
        if size(plane_chunk,3) > trace_length
            warning('plane_chunk is too long- trimming to match trace_hp')
            plane_chunk = plane_chunk(:,:,1:trace_length);
        end
        plane_chunk = imgaussfilt(plane_chunk, p.blur); %blur the loaded plane
        for m = 1:Nmask % par
            trace_px = plane_chunk(repmat(mask{m},[1,1,size(plane_chunk,3)])); %mask the plane with the 2D mask
            trace_px = reshape(trace_px,[],size(plane_chunk,3)); % and rehsape to Npx-by-t
            trace_px_hp = trace_px - movmedian(trace_px, p.hp_cutoff, 2); %high-pass filter the pixel traces
            %for each pixel in the plane, calculate the Pearson's correlation coefficient with the high-pass zproj trace
            Ctemp = zeros(1,sum(mask{m}(:)));
            for j = 1:size(trace_px,1)
                Ctemp(j) = corr(double(trace_px_hp(j,:))',trace_hp{m}');
            end
            % reshape correlation into an image
            c_im = zeros(size(mask{m}));
            c_im(mask{m}) = Ctemp;
            c_vol{m,z} = c_im;
        end
    end
end

fprintf('\nConstructing segmentation_data...');
segmentation_data = repmat( struct('mask_2d',[], 'mask_3d',[], 'mean_trace_raw',[], 'mean_trace_hp',[], 'correlation',[], 'xproj',[], 'yproj',[], 'zproj',[], 'keep',false), 1, Nmask );
% 
parfor m = 1:Nmask
    %reshape plane-by-plane correlations into volume
    cc = [c_vol{m,:}];
    cc = reshape(cc, sbxInfo.sz(1), sbxInfo.sz(2), []);
    %replace 0s with NaNs for projection
    c_nan = cc;
    c_nan(cc == 0) = NaN;
    %do X Y and Z projections
    c_nan_xproj = squeeze(nanmax(c_nan,[],1));
    c_nan_xproj = imresize(c_nan_xproj,[size(c_nan_xproj,1),p.xyproj_width]);
    c_nan_yproj = squeeze(nanmax(c_nan,[],2));
    c_nan_yproj = imresize(c_nan_yproj,[size(c_nan_yproj,1),p.xyproj_width]);
    c_nan_zproj = nanmax(c_nan,[],3);
    %threshold at p.prctile_cutoff (typically 90th prctile)
    threshold = prctile(c_nan,p.prctile_cutoff,'all');
    c_thresh = cc > threshold;
    %label objects in binarized image
    c_L = bwlabeln(c_thresh);
    %calculate volume of each object
    rp = regionprops3(c_L,'Volume');
    %remove objects below volume cutoff
    F = find(rp.Volume > p.volume_cutoff);
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
if ~isempty(p.save)
    fprintf('\nSaving %s...   ', p.save);
    save(p.save, 'segmentation_data', 'p', 'projScans', 'censScans', '-v7.3');
end
end