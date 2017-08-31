function video = grid_video(video_cell, nrows, ncols, rgap, cgap)
    % Function to create a stack of video patches over a 2D grid.
    % 
    % Inputs:
    %   video_cell: Cell of individual videos.
    %   nrows: Number of rows of videos.
    %   ncols: Number of columns of videos.
    %   rgap: Gap between rows, as a fraction of number of rows in video.
    %   cgap: Gap between columns, as a fraction of number of columns in
    %   video.
    %
    % Output:
    %   im: Grid of videos.
    
    nimg = numel(video_cell);
    
    vid0 = video_cell{1};
    if ndims(vid0) == 3
        [H, W, T] = size(vid0);
        nchan = 1;
    else
        [H, W, T, nchan] = size(vid0);
    end
    
    video = ones(round((1+rgap)*H)*nrows + round(rgap*H),...
                 round((1+cgap)*W)*ncols + round(cgap*W), nchan, T); 
    
    for idx = 1:nimg
        [x, y] = ind2sub([nrows, ncols], idx);
        row_idx = (x-1)*round((1+rgap)*H) + 1 + round(rgap*H);
        col_idx = (y-1)*round((1+cgap)*W) + 1 + round(rgap*W);
        try
            video(row_idx:row_idx+H-1, col_idx:col_idx+W-1, :, :) = ...
                reshape(video_cell{idx}, H, W, nchan, T);
        catch error
            keyboard
        end
    end
    
    % Just to be sure, squeeze it
    video = squeeze(video);
end