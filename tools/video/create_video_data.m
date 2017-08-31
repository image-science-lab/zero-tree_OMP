function Y = create_video_data(video_mat_cell, h, w, t, npatches)
    % Function to create video patches from video matrices.
    %
    % Input:
    %   video_mat_cell: Cell of video matrices.
    %   h: Height of patch
    %   w: Width of patch
    %   t: Number of frames per patch
    %   npatches: Total number of patches.
    %
    % Output:
    %   Y: h*w*t X npatches patches matrix.
    
    nvideos = numel(video_mat_cell);
    Ycell = cell(nvideos);
    npatches_per_vid = ceil(npatches/nvideos);
    
    parfor idx=1:nvideos
        Ycell{idx} = get_video_patches(video_mat_cell{idx}, h, w, t, npatches_per_vid);
    end
    Y = cell2mat(Ycell);
end