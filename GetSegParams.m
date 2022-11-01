function segParams = GetSegParams(sbxInfo)

paramsPath = sprintf('%s%s_seg_params.mat', sbxInfo.dir, sbxInfo.exptName ); 
if exist(paramsPath, 'file')
    fprintf('\nLoading %s', paramsPath);
    load(paramsPath); % , 'IP'
    if ~exist('segParams','var') && exist('IP','var')
        segParams = IP.Results;
    end
else
    error('%s does not exist', paramsPath);
end


end

