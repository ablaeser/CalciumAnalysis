function ExtractTrace3D(path,varargin)
    p = inputParser;
    addOptional(p,'hp_cutoff',51); %rolling hp filter cutoff
    addOptional(p,'prctile_cutoff',10); %rolling hp filter prctile cutoff
    addOptional(p,'t_max',20); %# timeframes to make mask-images
    addOptional(p,'z_max',5); %# zplanes to make mask-images
    addOptional(p,'summary_plots',true);
    addOptional(p,'ROI_plots',true); %decide whether or not to save ROI images
    
    parse(p,varargin{:});
    p = p.Results;
    
    [fdir, fname, ~] = fileparts(path);
    [~, Nx, Ny, Nz, Nt] = GetDimensions(path,[],[]);

    seg_data_path_final = strcat(fdir,filesep,fname,'_seg_data_final.mat');
    load(seg_data_path_final,'-mat'); % load segmentation_data struct
    keep_idx = [seg_data(:).keep];
    seg_data_keep = seg_data(keep_idx);
    clear seg_data;
    
    Nroi = numel(seg_data_keep);
    if Nroi == 0, error('No filters!'); end

    icaguidata = load(strcat(fdir,filesep,fname,'.icadata'),'-mat');
    icaguidata = icaguidata.icaguidata;


    %% shape all masks into Nroi x (Nx*Ny*Nz) array for easy matrix multiplication
    mask_3d_all = [seg_data_keep(:).mask_3d];
    mask_3d_all = reshape(mask_3d_all,Nx,Ny,[],Nz);
    mask_3d_all = permute(mask_3d_all,[3,1,2,4]);
    zExclude = find(squeeze(~any(mask_3d_all, [1,2,3])))'; % determine which planes were excluded from analysis
    mask_3d_all = reshape(mask_3d_all,Nroi,[]);
    %divide by the total number of voxels (mask_weight, i.e. size) so that matrix multiplying
    %by mask will give you an average, not a sum.
    mask_weight = sum(mask_3d_all,2);
    mask_3d_all = mask_3d_all ./ repmat(mask_weight,[1,Nx*Ny*Nz]);

    %% make 3d neuropils
    %{
    allMask = false(Nx,Ny,Nz);
    for r = 1:Nroi, allMask = allMask | seg_data_keep(r).mask_3d; end 
    %imshow( max(allMask,[],3) );
    allMaskInd = find(allMask); 
    %figure; 
    for i = 1:Nroi
        mask_3d = seg_data_keep(i).mask_3d;
        mask_3d_np = imdilate(mask_3d,strel('cuboid',[12,12,2])) & ~imdilate(mask_3d,strel('cuboid',[8,8,1]));
        mask_3d_np(:,:,zExclude) = false; % suppress planes that were excluded from segmentation 
        % suppress NP pixels that are also in an ROI itself
        tempNPind = find(mask_3d_np);
        [~,badNPind,~] = intersect( tempNPind, allMaskInd );
        if ~isempty(badNPind) 
            fprintf('\nROI %i: Suppressing %i neuropil voxels, out of %i (%2.1f pct), that overlap with an ROI', ...
                i, numel(badNPind), numel(tempNPind), 100*numel(badNPind)/numel(tempNPind) );
            mask_3d_np(badNPind) = false;
        end
        % Show the process
        %{
        sp(1) = subplot(1,3,1); imshow( max( mask_3d, [], 3), [] ); title(sprintf('ROI %i', i) );
        sp(2) = subplot(1,3,2); imshow( max(mask_3d_np, [], 3), [] ); title('Neuropil');
        labelMat = uint16( max( mask_3d + 2*mask_3d_np, [], 3) );
        sp(3) = subplot(1,3,3); imshow( label2rgb(labelMat)); title('ROI | Neuropil'); % imshow( max( mask_3d | mask_3d_np, [], 3), [] );
        impixelinfo; linkaxes(sp,'xy'); pause;
        %}
        
        % Save result
        seg_data_keep(i).mask_3d_np = mask_3d_np;
    end
    
    %repeat reshaping for for neuropil
    mask_3d_np_all = [seg_data_keep(:).mask_3d_np];
    mask_3d_np_all = reshape(mask_3d_np_all,Nx,Ny,[],Nz);
    mask_3d_np_all = permute(mask_3d_np_all,[3,1,2,4]);
    mask_3d_np_all = reshape(mask_3d_np_all,Nroi,[]);
    np_weight = sum(mask_3d_np_all,2);
    mask_3d_np_all = mask_3d_np_all ./ repmat(np_weight,[1,Nx*Ny*Nz]);
    %% extract 3d traces from 3d masks calculated in InterpolationChunked.m
    trace = cell(1,Nt); trace_np = cell(1,Nt);
    %don't let the parpool get too big, since it's a memory intensive operation
    delete(gcp);
    parpool(6);
    w = parfor_progressbar(Nt,'extracting 3d traces');
    parfor t = 1:Nt
        vol = double(pipe.io.sbxRead(path,Nz*(t-1)+1,Nz,1,[]));
        trace{t} = mask_3d_all * vol(:);
        trace_np{t} = mask_3d_np_all * vol(:);
        w.iterate(1);
    end
    close(w);
    trace = [trace{:}];
    trace_np = [trace_np{:}];
    %% make lowpass filter using rolling low-percentile 
    trace_sub = trace - trace_np;
    trace_sub_lp = movprctile(trace_sub', p.prctile_cutoff, p.hp_cutoff, 1)';
    trace_sub_hp = trace_sub - trace_sub_lp;
    %% add traces to final dataset
    for i = 1:Nroi
        seg_data_keep(i).trace_raw_3d = trace(i,:);
        seg_data_keep(i).trace_np_3d = trace_np(i,:);
        seg_data_keep(i).trace_sub_hp = trace_sub_hp(i,:);
    end
    %}
    seg_data_path_final = strcat(fdir,filesep,fname,'_seg_traces.mat');
    save(seg_data_path_final,'seg_data_keep', 'p','-v7.3');
    %% show all "good" masks in 2D
    if p.summary_plots
        all_masks_2d = zeros(Nx,Ny);
        for i = 1:Nroi
            mask_2d = seg_data_keep(i).mask_2d;
            all_masks_2d(mask_2d) = i;
        end
        all_masks_rgb = label2rgb(all_masks_2d,'prism','k');
        %
        trace_norm = rescale(trace_sub_hp);

        figure('Position',[700,400,500,800]);
        subplot(2,1,1);
        imshow(all_masks_rgb);
        title('ROI masks');
        subplot(2,1,2);
        c = prism(Nroi);
        hold on
        for i = 1:Nroi
            plot(trace_norm(i,:)+i,'LineWidth',1);

        end
        hold off
        xlim([1,Nt]);
        ylim([0.5,Nroi+1]);
        title('all ROI traces');
        xlabel('Volume #');
        ylabel('Axon #');
        saveas(gcf,strcat(fdir,filesep,fname,'_summary.pdf'));
        close;
    end
    %% plot each ROI and trace
    if p.ROI_plots
        ROIdir = strcat(fdir,filesep,'seg_data');
        try
            rmdir(ROIdir, 's');
        end
        mkdir(ROIdir);

        for i = 1:Nroi
            temptrace = trace_sub_hp(i,:);
            [~,t_idx] = maxk(temptrace,p.t_max);

            mask_3d = seg_data_keep(i).mask_3d;
            mask_3d_z = max(mask_3d,[],3);

            mask_size = reshape(mask_3d,[],Nz);
            mask_size = sum(mask_size,1);

            [~,z_idx] = maxk(mask_size,p.z_max);


            maxk_im = zeros(Nx,Ny,Nz,p.t_max);

            for t = 1:p.t_max
                maxk_im(:,:,:,t) = pipe.io.sbxRead(path,Nz*(t_idx(t)-1)+1,Nz,1,[]);
            end
            maxk_im = max(maxk_im,[],4);
            %
            %do zproj of top 5 most populous planes
            im_zproj = maxk_im(:,:,z_idx);
            im_zproj = max(im_zproj,[],3);

            %do x and y proj
            mask_max = max(mask_3d,[],3);
            mask_max_hull = bwconvhull(mask_max);
            mask_proj = repmat(mask_max_hull,[1,1,Nz]);
            maxk_im_proj = maxk_im .* mask_proj;

            im_xproj = squeeze(max(maxk_im_proj,[],1));
            im_xproj_stretch = imresize(im_xproj',[100,Ny]);
            mask_3d_x = squeeze(max(mask_3d,[],1));
            mask_3d_x_stretch = imresize(mask_3d_x',[100,Ny]);

            im_yproj = squeeze(max(maxk_im_proj,[],2));
            im_yproj_stretch = imresize(im_yproj,[Nx,100]);
            mask_3d_y = squeeze(max(mask_3d,[],2));
            mask_3d_y_stretch = imresize(mask_3d_y,[Nx,100]);
            %
            %make montage of x- y- and z-projected volumes
            montage = zeros(Nx+100,Ny+100);
            montage(1:Nx,1:Ny) = im_zproj;
            montage(Nx+1:end,1:Ny) = im_xproj_stretch;
            montage(1:Nx,Ny+1:end) = im_yproj_stretch;

            montage = rescale(montage);
            montage = imadjust(montage,[0,0.5]);
            montage = im2uint16(montage);

            %make montage of the corresponding 3d masks
            mask_montage = zeros(Nx+100,Ny+100);
            mask_montage(1:Nx,1:Ny) = mask_3d_z;
            mask_montage(Nx+1:end,1:Ny) = mask_3d_x_stretch;
            mask_montage(1:Nx,Ny+1:end) = mask_3d_y_stretch;

            % combine and display the montages
            montage_rgb = repmat(uint16(montage),[1,1,3]);
            mask_outline= zeros(size(montage_rgb));
            mask_outline = imdilate(mask_montage,strel('disk',7)) & ~imdilate(mask_montage,strel('disk',6));
            mask_outline = logical(mask_outline);
            montage_rgb(mask_outline) = intmax('uint16');

            % make plot and save it
            F = figure('Position',[500,400,800,800]);
            subplot(5,1,1:4);
            imshow(montage_rgb,[]);
            title(strcat("ROI number ",num2str(i)));
            subplot(5,1,5);
            plot(temptrace);
            xlabel('volume #');
            title('high-pass fluorescence');

            saveas(F,strcat(ROIdir,filesep,fname,'_ROI_',sprintf('%03d.pdf',i)));
            close(F);
        end
    end
end