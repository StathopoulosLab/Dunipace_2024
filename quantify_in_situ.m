function [im, data] = quantify_in_situ(varargin)
%QUANTIFY_IN_SITU Compute mean intensity
% 
%   Inputs
%       varargin:
%           1): data structure returned from this function
% 
%   Outputs
%       data: a structure containing the following fields
%           folder: the folder path containing the opened file
%           name: part of the file name
%           condition: name of experimental condition to group
%               similarly treated experiments
%           max_projection: maximum intensity projection of a z-stack
%           sum_projection: sum of intensities projection of a z-stack
%           pix_len: physical length of a pixel in microns
%           mask_max_sig: mask of the signal from the max projection
%           mask_max_bg: mask of the background from the max projection
%           mask_sum_sig: mask of the signal from the sum projection
%           mask_sum_bg: mask of the background from the sum projection
%           sum_I: total sum of all intensities in signal
%           norm_sum_I: normalized total sum of all intensities, calculated
%               by subtracting the average background from each pixel,
%               using sum projection
%           norm_I: normalized average intensity from max projection,
%               calculated by subtracting the average background intensity
%           norm_A: normalized area of the signal from max projection,
%               calculated by dividing by the area of the embryo
%           norm_Maj: normalized major axis of the signal from max
%               projection, calculated by dividing by the major axis of the
%               embryo
%           norm_Min: normalized minor axis of the signal from max
%               projection, calculated by dividing by the major axis of the
%               embryo
%           Maj_im_max: major axis of the signal from max projection, found
%               by fitting an ellipse to the segmented signal
%           Min_im_max: minor axis of the signal from max projection, found
%               by fitting an ellipse to the segmented signal
%           A_im_max: area of the signal from max projection, in square
%               microns
%           I_im_max: mean intensity of the signal from max projection
%           A_bg_max: area of the background signal from max projection, in
%               square microns
%           I_bg_max: mean intensity of the background signal from max
%               projection
%           Maj_em_max: major axis of the whole embryo from max projection,
%               found by fitting an ellipse to the whole embryo
%           Min_em_max: minor axis of the whole embryo from max projection,
%               found by fitting an ellipse to the whole embryo
%           A_em_max: area of the whole embryo from max projection, in
%               square microns
%           I_em_max: mean intensity of the whole embryo from max
%               projection
%           Maj_im_sum: major axis of the signal from sum projection, found
%               by fitting an ellipse to the segmented signal
%           Min_im_sum: minor axis of the signal from sum projection, found
%               by fitting an ellipse to the segmented signal
%           A_im_sum: area of the signal from sum projection, in square
%               microns
%           I_im_sum: mean intensity of the signal from sum projection
%           A_bg_sum: area of the background signal from sum projection, in
%               square microns
%           I_bg_sum: mean intensity of the background signal from sum
%               projection
%           Maj_em_sum: major axis of the whole embryo from sum projection,
%               found by fitting an ellipse to the whole embryo
%           Min_em_sum: minor axis of the whole embryo from sum projection,
%               found by fitting an ellipse to the whole embryo
%           A_em_sum: area of the whole embryo from sum projection, in
%               square microns
%           I_em_sum: mean intensity of the whole embryo from sum
%               projection
%           n_sig_obj_max: number of unconnected objects detected in the
%               signal of max projection
%           n_bg_obj_max: number of unconnected objects detected in the
%               background signal of max projection
%           n_sig_obj_sum: number of unconnected objects detected in the
%               signal of sum projection
%           n_bg_obj_sum: number of unconnected objects detected in the
%               background signal of sum projection
%           T_bg_max: threshold for segmenting the background for the max
%               projection
%           T_sig_max: threshold for segmenting the signal for the max
%               projection
%           sd_max: standard deviation for performing a gaussian blur for
%               the max projection
%           r_bg_max: radius of a disk structuring element to
%               morphologically close the image after thresholding the
%               background for the max projection
%           r_sig_max: radius of a disk structuring element to
%               morphologically close the image after thresholding the
%               signal for the max projection
%           min_size_obj_max: size, in pixels, to keep in image (removes
%               objects less than this from mask) for the max projection
%           T_bg_sum: threshold for segmenting the background for the sum
%               projection
%           T_sig_sum: threshold for segmenting the signal for the sum
%               projection
%           sd_sum: standard deviation for performing a gaussian blur for
%               the sum projection
%           r_bg_sum: radius of a disk structuring element to
%               morphologically close the image after thresholding the
%               background for the sum projection
%           r_sig_sum: radius of a disk structuring element to
%               morphologically close the image after thresholding the
%               signal for the sum projection
%           min_size_obj_sum: size, in pixels, to keep in image (removes
%               objects less than this from mask) for the sum projection
%           use_bg_chnl_max: enter channel number if background
%               segmentation should use another channel for the max
%               projection
%           use_bg_chnl_sum: enter channel number if background
%               segmentation should use another channel for the sum
%               projection
%           
% 
%   Overview
%       This function makes max and sum projections; segments the signal,
%       whole embryo, and background; and reports the mean intensity or sum
%       intensity, area, major axis length, and minor axis length as well
%       as the normalized values. All lengths are in microns and areas are
%       in square microns. The whole embryo can be segmented using either
%       the signal channel or another channel (specified using
%       use_bg_chnl_max or use_bg_chnl_sum). If using the signal channel, a
%       lower threshold (T_bg_max or T_bg_sum) should be used. The signal
%       is segmented using T_sig_max and T_sig_sum. The background signal
%       is determined by removing the signal from the whole embryo.
%       Normalized mean intensity (norm_I) is calculated by subtracting the
%       mean intensity of the background from the mean intensity of the
%       signal. Normalized sum of intensities (norm_sum_I) is calculated by
%       dividing the sum of intensities by the average background intensity
%       mutiplied by the number of pixels in the sum. Normalized area
%       (norm_A) is calulated by dividing the area of the signal by the
%       area of the embryo. Normalized length (norm_Maj) is calulated by
%       dividing the length of the signal by the length of the embryo.
%       Normalized width (norm_Min) is calculated by dividing the width of
%       the signal by the width of the embryo.
    
    % Determine if first input is a structure and calculate n
    if ~isempty(varargin) && isstruct(varargin{1})
        n = size(varargin{1}, 2);
        im = varargin{1};
        data = varargin{2};
    else
        % Use menu to select files
        [name, folder] = uigetfile({'*.czi', 'CZI files (*.czi)'},...
                'Select the microscope images', 'Multiselect', 'on');

%         % Use menu to select files
%         [name, folder] = uigetfile({'*.tif', 'TIF files (*.tif)'},...
%                 'Select the microscope images', 'Multiselect', 'on');
        
        % Calculate n
        n = size(name, 2);
        
        % If parameters were inputted use those
        if ~isempty(varargin)
            params = varargin{1};
        % Else use the defaults
        else
            params = [0.04, 0.050, 3, 10, 5, 1000, 3, 1;...
                      0.02, 0.005, 3, 10, 5, 1000, 3, 1];
%             params = [0.04, 0.050, 3, 10, 5, 1000, 4, 1;...
%                       0.02, 0.005, 3, 10, 5, 1000, 4, 1];
        end
        
        % Make empty structure
        im = struct('folder', cell(1,n),...
                    'name', [],...
                    'condition', 'empty',...
                    'max_projection', [],...
                    'avg_projection', [],...
                    'n_z', [],...
                    'pix_len', [],...
                    'mask_max_sig', [],...
                    'mask_max_bg', [],...
                    'mask_avg_sig', [],...
                    'mask_avg_bg', []);

        data = struct('folder', cell(1,n),...
                      'name', [],...
                      'condition', 'empty',...
                      'n_z', [],...
                      'pix_len', [],...
                      'sum_I_avg', [],...
                      'norm_sum_I_avg', [],...
                      'norm_I_avg', [],...
                      'norm_A_avg', [],...
                      'norm_Maj_avg', [],...
                      'norm_Min_avg', [],...
                      'norm_MaxFeret_avg', [],...
                      'norm_MinFeret_avg', [],...
                      'Maj_im_max', [],...
                      'Min_im_max', [],...
                      'MaxFeret_im_max', [],...
                      'MinFeret_im_max', [],...
                      'A_im_max', [],...
                      'I_im_max', [],...
                      'A_bg_max', [],...
                      'I_bg_max', [],...
                      'Maj_em_max', [],...
                      'Min_em_max', [],...
                      'MaxFeret_em_max', [],...
                      'MinFeret_em_max', [],...
                      'A_em_max', [],...
                      'I_em_max', [],...
                      'Maj_im_avg', [],...
                      'Min_im_avg', [],...
                      'MaxFeret_im_avg', [],...
                      'MinFeret_im_avg', [],...
                      'A_im_avg', [],...
                      'I_im_avg', [],...
                      'A_bg_avg', [],...
                      'I_bg_avg', [],...
                      'Maj_em_avg', [],...
                      'Min_em_avg', [],...
                      'MaxFeret_em_avg', [],...
                      'MinFeret_em_avg', [],...
                      'A_em_avg', [],...
                      'I_em_avg', [],...
                      'n_sig_obj_max', [],...
                      'n_bg_obj_max', [],...
                      'n_sig_obj_avg', [],...
                      'n_bg_obj_avg', [],...
                      'T_bg_max', params(1,1),...
                      'T_sig_max', params(1,2),...
                      'sd_max', params(1,3),...
                      'r_bg_max', params(1,4),...
                      'r_sig_max', params(1,5),...
                      'min_size_obj_max', params(1,6),...
                      'use_bg_chnl_max', params(1,7),...
                      'sig_chnl_max', params(1,8),...
                      'T_bg_avg', params(2,1),...
                      'T_sig_avg', params(2,2),...
                      'sd_avg', params(2,3),...
                      'r_bg_avg', params(2,4),...
                      'r_sig_avg', params(2,5),...
                      'min_size_obj_avg', params(2,6),...
                      'use_bg_chnl_avg', params(2,7),...
                      'sig_chnl_avg', params(2,8),...
                      'C_em_max', [],...
                      'C_em_avg', [],...
                      'C_im_max', [],...
                      'C_im_avg', [],...
                      'angle_em_max', [],...
                      'angle_em_avg', [],...
                      'angle_im_max', [],...
                      'angle_im_avg', [],...
                      'varI_im_max', [],...
                      'varI_im_avg', [],...
                      'disp_ind_im_max', [],...
                      'disp_ind_im_avg', [],...
                      'stdI_im_max', [],...
                      'stdI_im_avg', [],...
                      'coeff_var_im_max', [],...
                      'coeff_var_im_avg', []);
    end
    
    % For each file
    for i = 1:n
        % If inputs weren't structures
        if isempty(varargin) || ~isstruct(varargin{1})
            % Open image
            [im(i).folder, im(i).name, im(i).condition,...
                im(i).max_projection, im(i).avg_projection,...
                im(i).pix_len, im(i).n_z] = open_im(name{i}, folder);

            data(i).folder = im(i).folder;
            data(i).name = im(i).name;
            data(i).condition = im(i).condition;
            data(i).pix_len = im(i).pix_len;
            data(i).n_z = im(i).n_z;
        end
        
        % Segment the maximum and sum projections
        [im(i), data(i)] = segment_projection(im(i), data(i), 'max_projection', 1);
        [im(i), data(i)] = segment_projection(im(i), data(i), 'avg_projection', 2);
        
        % Calculate normalized values for max and sum projections
        data(i) = normalize_mean_I(im(i), data(i));
    end
end

function [path, embryo_number, condition, im_max, im_avg,...
    pix_len, Z] = open_im(name, folder)
%OPEN_IMG Open a czi with a z-stack, a time series, and channels
% 
%   Inputs
%       name: file name of the image
%       folder: folder with all the image files
% 
%   Outputs
%       path: the folder path containing the opened file
%       embryo_number: part of the file name before the first space
%       condition: condition
%       img_max: raw max z-projection of images
%       img_sum: raw sum z-projection of images
%       pix_len: length of a pixel in microns
% 
%   Overview
%       A z-projection is made and the red color channel is selected.
%       Images can be any number of channels
%       If the wrong channel is selected adjust the 1 (third index to img)
%       in img = squeeze(max(squeeze(img(:,:,1,:,:)), [], 3)); to choose
%       the correct channel.
            
    % Construct full path
    path = fullfile(folder, name);

    % Split and save part of file name before first space as unique
    % identifier
    file_ext = strsplit(name, '.');

    if strcmp(file_ext{end}, 'czi')
        file_name_parts = strsplit(name, {' ', '.czi'});
    elseif strcmp(file_ext{end}, 'tif')
        file_name_parts = strsplit(name, {' ', '.tif'});
    end

    embryo_number = file_ext{1};
    condition = file_name_parts{1};
    
    % Opens images using Bio-Formats for MATLAB
    % https://docs.openmicroscopy.org/bio-formats/6.1.0/users/matlab/index.html
    im = bfopen(path);
    
    % Save sizes of images in all dimensions, including time and color
    % channels
    X = im{1,4}.getPixelsSizeX(0).getValue();
    Y = im{1,4}.getPixelsSizeY(0).getValue();
    Z = im{1,4}.getPixelsSizeZ(0).getValue();
    T = im{1,4}.getPixelsSizeT(0).getValue();
    C = im{1,4}.getPixelsSizeC(0).getValue();
    
    % The physical length of a pixel
    pix_len = im{1,4}.getPixelsPhysicalSizeX(0).value(...
                      ome.units.UNITS.MICROMETER); % in µm
    pix_len = pix_len.doubleValue();
    
    if strcmp(file_ext{end}, 'czi')
        % Reshape image data to match dimensions, X, Y, channels, z, time
        im = permute(reshape(cat(3, im{1,1}{:,1}), Y, X, C, Z, T),[1,2,4,5,3]);
    elseif strcmp(file_ext{end}, 'tif')
        im = reshape(cat(3, im{1,1}{:,1}), Y, X, C, Z, T);
        temp_im = zeros(Y,X,1,Z,T);
        for j =1:T
            for i = 1:Z
                temp_im(:,:,1) = rgb2gray(im(:,:,:,i,j));
            end
        end

        im = uint8(permute(temp_im, [1,2,4,5,3]));
    end

    % Select the channel (usually red channel), remove the channel
    % dimension and make a z-projection
    im_max = max(im, [], 3);
    im_avg = mean(im, 3);
end

function [im_in, data] = segment_projection(im_in, data, projection, i)
%SEGMENT_MAX_PROJECTION Segment image using max projection.
%   
%   Input
%       data: data from parent function
%   
%   Output
%       data: data updated with image segmentation data
%
%   Overview
%       Segments an image of a drosophila embryo, using thresholding
    
    % Field names for storing max or sum projection data
    field = {'n_bg_obj_max', 'n_bg_obj_avg';...
             'T_bg_max', 'T_bg_avg';...
             'T_sig_max', 'T_sig_avg';...
             'sd_max', 'sd_avg';...
             'r_bg_max', 'r_bg_avg';...
             'r_sig_max', 'r_sig_avg';...
             'min_size_obj_max', 'min_size_obj_avg';...
             'use_bg_chnl_max', 'use_bg_chnl_avg';...
             'sig_chnl_max', 'sig_chnl_avg';...
             'n_sig_obj_max', 'n_sig_obj_avg'};
    
    % Save projection image from channel specified
    im = im_in.(projection)(:,:,:,:, data.(field{9,i}));
    
    % If nonzero, use channel specified for the background
    if data.(field{8,i})
        bg_im = im_in.(projection)(:,:,:,:, data.(field{8,i}));
    % Else use the same projection with the signal
    else
        bg_im = im;
    end
    
    % If sum projection
    if i >= 2
        % Scale by maximum for 16 bit
        im_sc = im./(2.^16 - 1);
        bg_im_sc = bg_im./(2.^16 - 1);
    % Else don't scale
    else
        im_sc = im;
        bg_im_sc = bg_im;
    end
    
    % Gaussian blur with standard deviation specified
    im_blur_bg = imgaussfilt(bg_im_sc, data.(field{4,i}));
    
    % Segment entire embryo from blurred image using thresholding
    mask_bg = imbinarize(im_blur_bg, data.(field{2,i}));
    
    % Morophlogically close image to try to get one solid shape, remove any
    % small objects that are not connected to the entire embryo
    se = strel('disk', data.(field{5,i}));
    mask_bg = imopen(mask_bg, se);
    mask_bg = imfill(mask_bg, "holes");
    mask_bg = bwareaopen(mask_bg, data.(field{7,i}));
    
    % Gaussian blur with standard deviation specified
    im_blur = imgaussfilt(im_sc, data.(field{4,i}));
    
    % Segment signal from blurred image using thresholding
    mask_sig = imbinarize(im_blur, data.(field{3,i}));
    
    % Morophlogically close image to try to get solid shapes, remove any
    % small objects that are not connected to the signal
    se = strel('disk', data.(field{6,i}));
    mask_sig = imclose(mask_sig,se);
    mask_sig = imfill(mask_sig, "holes");
    mask_sig = bwareaopen(mask_sig, data.(field{7,i}));
    
    % Get properties of the segmented object. Use a logical mask to get
    % properties for individual (not connected) objects, or use
    % double(mask) to get properties where all segmented regions are
    % treated as one object
    stats_embryo = regionprops(mask_bg, im, 'Area',...
                               'MajorAxisLength', 'MinorAxisLength',...
                               'MeanIntensity', 'Centroid','Orientation',...
                               'MaxFeretProperties', 'MinFeretProperties');
    
    % Save the number of objects segmented from background, indicative of
    % additional embryos being in the field of view
    data.(field{1,i}) = size(stats_embryo, 1);
    
    % If more than one object is detected in the background, pick the
    % embryo closest to the center
    if size(stats_embryo, 1) > 1
        % Find the closest object/embryo to the center of the image
        [~, ind_bg] = min(calc_dist(cat(1,stats_embryo.Centroid),...
                                   (size(im) + 1)./2));
        
        % Keep only the data for the object/embryo closest to the center of
        % the image
        stats_embryo = stats_embryo(ind_bg);
        
        % Remove other objects/embryos from the mask
        L_bg = labelmatrix(bwconncomp(mask_bg));
        mask_bg(~(L_bg == ind_bg)) = false;
    end

    % Remove any segmented signal from additional objects/embryos
    mask_sig(~mask_bg) = false;
    mask_sig = bwareaopen(mask_sig, data.(field{7,i}));

    stats_sig = regionprops(double(mask_sig), im, 'Area',...
                        'MajorAxisLength', 'MinorAxisLength',...
                        'MeanIntensity', 'Centroid', 'Orientation',...
                        'MaxFeretProperties', 'MinFeretProperties');

    data.(field{10,i}) = size(stats_sig, 1);
    
    % Get image properties
    [im_in, data] = get_im_props(im_in, data, im, stats_embryo,...
        mask_bg, stats_sig, mask_sig, i);
end

function d = calc_dist(x, y)
%CALC_DIST Calculates the distance between two points in n-dimensions
%   
%   Input
%       x: corrdinates for an array of points (can be mxn in size)
%       y: corrdinates for a point (should be mxn in size)
%   
%   Output
%       d: array of distances bewteen the points in x and the point in y
%
%   Overview
%       Calculates the distance in n-dimensial space by subtracting x and y
%       by applying element-wise operation to the two arrays with implicit
%       expansion enabled. This will subtract y from each row of x if x is
%       mxn and y is 1xn, where n is the number of dimensions. These valued
%       are then squared and summed upon the second dimension of the array.
%       Finally the square root is taken. This gives the distance formula,
%       d = sqrt((x1-x2)^2+(y1-y2)^2) but for n-dimensions and for
%       mutiple points in x from one point y.

    d = sqrt(sum(bsxfun(@minus, x, y).^2,2));
end

function [im_in, data] = get_im_props(im_in, data, im, stats_embryo,...
    mask_bg, stats_sig, mask_sig, i)
%GET_IMAGE_PROPERTIES Get image properties
%   
%   Input
%       data: the data structure from the main function
%       im: a z-projection (max or sum)
%       stats_embryo: a structure with image properties for the segmented
%           embryo
%       mask_bg: a logical mask for the segmented background signal
%       stats_sig: a structure with image properties for the segmented
%           signal
%       mask_sig: a lgocal mask for the segmented signal
%   
%   Output
%       data: the data structure from the main function, updated with the
%           image properties
%
%   Overview
%       This function takes the structure returned from region props, and
%       saves the values of those properties into the data structure from
%       the main function
    
    % Names of fields for saving image properties
    field = {'mask_max_sig', 'mask_avg_sig';...
             'mask_max_bg', 'mask_avg_bg';...
             'A_bg_max', 'A_bg_avg';...
             'I_bg_max', 'I_bg_avg';...
             'A_em_max', 'A_em_avg';...
             'Maj_em_max', 'Maj_em_avg';...
             'Min_em_max', 'Min_em_avg';...
             'MaxFeret_em_max', 'MaxFeret_em_avg';...
             'MinFeret_em_max', 'MinFeret_em_avg';...
             'I_em_max', 'I_em_avg';...
             'A_im_max', 'A_im_avg';...
             'Maj_im_max', 'Maj_im_avg';...
             'Min_im_max', 'Min_im_avg';...
             'MaxFeret_im_max', 'MaxFeret_im_avg';...
             'MinFeret_im_max', 'MinFeret_im_avg';...
             'I_im_max', 'I_im_avg';...
             'C_em_max', 'C_em_avg';...
             'C_im_max', 'C_im_avg';...
             'angle_em_max', 'angle_em_avg';...
             'angle_im_max', 'angle_im_avg';...
             'varI_im_max', 'varI_im_avg';...
             'disp_ind_im_max', 'disp_ind_im_avg';...
             'stdI_im_max', 'stdI_im_avg';...
             'coeff_var_im_max', 'coeff_var_im_avg'};

    % If signal and background were detected
    if (size(stats_embryo, 1) ~= 0) && (size(stats_sig, 1) ~= 0)
        % Remove signal from segmented embryo to get background
        mask_bg(mask_sig) = false;

        % Save the masks
        im_in.(field{1,i}) = mask_sig;
        im_in.(field{2,i}) = mask_bg;
        
        % Calculate area, convert to microns squared, and calculate mean
        % intensity for background
        data.(field{3,i}) = sum(mask_bg(:)) .* data.pix_len.^2;
        data.(field{4,i}) = mean(im(mask_bg));
        
        % Convert area and lengths/widths to microns squared or microns,
        % and save the values, including mean intensity, for the embryo
        data.(field{5,i}) = stats_embryo.Area .* data.pix_len.^2;
        data.(field{6,i}) = stats_embryo.MajorAxisLength .* data.pix_len;
        data.(field{7,i}) = stats_embryo.MinorAxisLength .* data.pix_len;
        data.(field{8,i}) = stats_embryo.MaxFeretDiameter .* data.pix_len;
        data.(field{9,i}) = stats_embryo.MinFeretDiameter .* data.pix_len;

        data.(field{10,i}) = stats_embryo.MeanIntensity;
        
        % Convert area and lengths/widths to microns squared or microns,
        % and save the values including mean intensity for the signal
        data.(field{11,i}) = stats_sig.Area .* data.pix_len.^2;
        data.(field{12,i}) = stats_sig.MajorAxisLength .* data.pix_len;
        data.(field{13,i}) = stats_sig.MinorAxisLength .* data.pix_len;
        data.(field{14,i}) = stats_sig.MaxFeretDiameter .* data.pix_len;
        data.(field{15,i}) = stats_sig.MinFeretDiameter .* data.pix_len;

        data.(field{16,i}) = stats_sig.MeanIntensity;

        % Save the centers and orientation
        data.(field{17,i}) = stats_embryo.Centroid;
        data.(field{18,i}) = stats_sig.Centroid;
        data.(field{19,i}) = stats_embryo.Orientation;
        data.(field{20,i}) = stats_sig.Orientation;
        data.(field{21,i}) = var(double(im(mask_sig)));
        data.(field{22,i}) = data.(field{17,i}) / data.(field{10,i});
        data.(field{23,i}) = std(double(im(mask_sig)));
        data.(field{24,i}) = data.(field{19,i}) / data.(field{10,i});
    else
        im_in.(field{1,i}) = [];
        im_in.(field{2,i}) = [];
        for j = 3:size(field,1)
            data.(field{j,i}) = [];
        end
    end
end

function [data] = normalize_mean_I(im_in, data)
%NORMALIZE_MEAN Subtract background levels from signal levels.
% 
%   Input
%       data: the data structure in the main function
% 
%   Output
%       data: the data structure from the main function that has been
%             updated with normalized data
% 
%   Overview
%       This function takes the data structure from the main function and
%       normalizes the data. Specifically, this function subtracts out the
%       background signal for maximum and sum projections. In addition, the
%       length, width, and area of the signal is normalized by the length,
%       width, or area of the entire embryo.

            
    % Calculate normalized intensity by subtracting and/or dividing by
    % background
    data.norm_I_avg = (data.I_im_avg - data.I_bg_avg);
    
    im = im_in.avg_projection(:,:, data.('sig_chnl_avg')) * data.n_z;
    % Calculate sum of intensities of all pixels
    data.sum_I_avg = sum(im(im_in.mask_avg_sig));
    
    % Calculate number of pixels for normalizing sum
    n = sum(im_in.mask_avg_sig(:));
    
    % Subtract the average background intensity for each pixel
    data.norm_sum_I_avg = data.sum_I_avg - (n .* data.I_bg_avg * data.n_z);
    
    % Normalize area of signal by area of the embryo
    data.norm_A_avg = data.A_im_avg ./ data.A_em_avg;
    
    % Normalize length of signal by length of the embryo
    data.norm_Maj_avg = data.Maj_im_avg ./ data.Maj_em_avg;
    
    % Normalize width of signal by width of the embryo
    data.norm_Min_avg = data.Min_im_avg ./ data.Min_em_avg;

    % Normalize length of signal by length of the embryo
    data.norm_MaxFeret_avg = data.MaxFeret_im_avg ./ data.MaxFeret_em_avg;
    
    % Normalize width of signal by width of the embryo
    data.norm_MinFeret_avg = data.MinFeret_im_avg ./ data.MinFeret_em_avg;
end