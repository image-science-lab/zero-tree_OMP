function im = grid_images(imcell, nrows, ncols, rgap, cgap)
    % Function to plot a stack of images over a 2D grid.
    % 
    % Inputs:
    %   imcell: Cell of individual images.
    %   nrows: Number of rows of images.
    %   ncols: Number of columns of images.
    %   rgap: Gap between rows, as a fraction of number of rows in image.
    %   cgap: Gap between columns, as a fraction of number of columns in
    %   image.
    %
    % Output:
    %   im: Grid of images.
    
    nimg = numel(imcell);
    
    im0 = imcell{1};
    if ndims(im0) == 2
        [H, W] = size(im0);
        nchan = 1;
    else
        [H, W, nchan] = size(im0);
    end
    
    im = ones(round((1+rgap)*nrows*H), round((1+cgap)*ncols*W), nchan); 
    
    for idx = 1:nimg
        [y, x] = ind2sub([ncols, nrows], idx);
        row_idx = (x-1)*(1+rgap)*H + 1;
        col_idx = (y-1)*(1+cgap)*W + 1;
        im(row_idx:row_idx+H-1, col_idx:col_idx+W-1, :) = ...
            reshape(imcell{idx}, H, W, nchan);
    end
    
    % Just to be sure, squeeze it
    im = squeeze(im);
end