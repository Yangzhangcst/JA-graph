function class_number = AP_segmentation(img, KDthreshold, n, m, bins, smoothing)
% AP_segmentation takes four inputs for segmenting PET images based on a
% similarity function that calculates the affinity (similarity) between the
% data points on the histogram. 
%
% PARAMETERS     Value
%   img          A 2D or 3D image. The histogram will be estimated from this image
%
%   bins         The number of bins for the histogram estimation,
%                Default = 256
%
%   KDthreshold  The threshold for removing noise based on the kernel
%                desity estimation via diffusion. Using the maximum density
%                of the data points on the image, any points that are less
%                than the threshold times this value is most likely an
%                outlier and is removed for more robust histogram estimation.
%
%   n,m          The parameters for calculating the affinity function of the . 
%                data points n,m must be greater than or equal to zero (and 
%                at least one must be nonzero). 
%
%                Default values are n = 3, m = 1.
%   smoothing    The window size of the exponential smoothing. Default = 20
%
%
%   Outputs:
%    class_number The number of the clustering.
%
    if nargin < 5 bins = 256; smoothing = 20; end

    %%Estimate the histogram of the 2D or 3D image
    H_total = zeros(1,bins)';
    bands = size(img, 3);
    two_D_histogram = zeros(bins,1,bands);
    
    for i = 1:bands        
        temp_img = img(:,:,i);
        two_D_histogram(:,:,i) = hist(temp_img(:), bins)';
        %%%% Use 3D Histogram if more than one slice is given%%%%%
        H_total = H_total + two_D_histogram(:,:,i);          
    end
    %%%%% Histogram Smoothing Step %%%%%
    [H] = histogram_smoothing(H_total, KDthreshold, smoothing);
    
    %%%%% Affinity Calculation Step %%%%%
    binary = affinity_calculation(H, n, m);
    class_number = length(binary(2:end));
end



