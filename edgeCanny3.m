function [eout, dx, dy, dz, dxthin, dythin, dzthin] = edgeCanny3(image_scan, ratios, thresh, sigma)

if nargin < 4
    sigma = sqrt(2); 
if nargin < 2
    % magic numbers 
    pNotEdges = .7;   % Used for selecting thresholds
    threshRatio = .4; % Low thresh is this fraction of the high.
else
    pNotEdges = ratios(1); 
    threshRatio = ratios(2);
end
end

% Compute smoothed numerical gradient of image I along y (vertical)
% direction. GY corresponds to dG/dy, where G is the Gaussian Smoothed
% version of image I.

[smoothing_kernel, derivative_kernel] = getGaussKernel(sigma);

dx = imfilter(image_scan, permute(smoothing_kernel, [2,1,3]), 'conv', 'replicate');
dx = imfilter(dx, derivative_kernel, 'conv', 'replicate');
dy = imfilter(image_scan, smoothing_kernel, 'conv', 'replicate'); 
dy = imfilter(dy, permute(derivative_kernel, [2,1,3]), 'conv', 'replicate');
dz = imfilter(image_scan, permute(smoothing_kernel, [3,2,1]), 'conv', 'replicate');
dz = imfilter(dz, permute(derivative_kernel, [3,1,2]), 'conv', 'replicate');

% magnitude of gradient 
magGrad = sqrt( dx.^2 + dy.^2 + dz.^2 );
% Normalize for threshold selection
magmax = max(magGrad(:));
if magmax > 0
    magGrad = magGrad / magmax;
end

% determine Hysteresis thresholds 
if nargin < 3
    [~, highThresh, lowThresh] = histThreshold(magGrad, pNotEdges, threshRatio);
else
    lowThresh = thresh(1); highThresh = thresh(2);
end

% Perform Non-Maximum Suppression Thining and Hysteresis Thresholding of Edge
% Strength
eout = thinAndThreshold(dx,dy,dz, magGrad, lowThresh, highThresh);
% thresh = [lowThresh highThresh]; % output thresholds if requested 
dxthin = dx(eout); dythin = dy(eout); dzthin = dz(eout);

% ----------- local functions ------------------------------------------

    function H = thinAndThreshold(dx,dy,dz, magGrad, lowThresh, highThresh)
        % Perform Non-Maximum Suppression Thining and Hysteresis Thresholding of
        % Edge Strength
        
        % We will accrue indices which specify ON pixels in strong edgemap
        % The array e will become the weak edge map.
        
        E = false(size(magGrad)); % or false(size(image_scan)) ?
        idxStrong = [];
        for dir = 1:24
            idxLocalMax = cannyFindLocalMaxima(dir, dx,dy,dz, magGrad);
            idxWeak = idxLocalMax(magGrad(idxLocalMax) > lowThresh);
            E(idxWeak) = 1;
            idxStrong = [idxStrong; idxWeak(magGrad(idxWeak) > highThresh)];
        end
        
        if ~isempty(E)
            marker = find(magGrad>highThresh & E);
            
            if ~isempty(rstrong) % result is all zeros if idxStrong is empty
                H = imreconstruct(marker, E, 26);
            else
                H = false(size(E));
            end
        else
            H = false(size(E));
        end
    end

    function [gaussKernel3, derivGaussKernel, gaussKernel] = getGaussKernel(sigma)
% Determine filter length
filterExtent = ceil(4*sigma);
x = -filterExtent:filterExtent;

% Create 1-D Gaussian Kernel
c = 1/(sqrt(2*pi)*sigma);
gaussKernel = c * exp(-(x.^2)/(2*sigma^2));

% Normalize to ensure kernel sums to one
gaussKernel = gaussKernel/sum(gaussKernel);

% Create 1-D Derivative of Gaussian Kernel
derivGaussKernel = gradient(gaussKernel);

% Normalize to ensure kernel sums to zero
negVals = derivGaussKernel < 0;
posVals = derivGaussKernel > 0;
derivGaussKernel(posVals) = derivGaussKernel(posVals)/sum(derivGaussKernel(posVals));
derivGaussKernel(negVals) = derivGaussKernel(negVals)/abs(sum(derivGaussKernel(negVals)));  

% make the gauss kernel suitable for 3D
[X,Y] = meshgrid(gaussKernel); 
Z = X.*Y; 
gaussKernel3 = permute(Z, [3,1,2]);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Local Function : cannyFindLocalMaxima
%
function idxLocalMax = cannyFindLocalMaxima(direction,ix,iy,iz,mag)
%
% This sub-function helps with the non-maximum suppression in the Canny
% edge detector.  The input parameters are:
%
%   direction - the index of which direction the gradient is pointing,
%               read from the diagram below. direction is 1, 2, 3, or 4.
%   ix        - input image filtered by derivative of gaussian along x
%   iy        - input image filtered by derivative of gaussian along y
%   mag       - the gradient magnitude image
%
%    there are 4 cases: [24 cases in 3D]
%
%                         The X marks the pixel in question, and each
%         3     2         of the quadrants for the gradient vector
%       O----0----0       fall into two cases, divided by the 45
%     4 |         | 1     degree line.  In one case the gradient
%       |         |       vector is more horizontal, and in the other
%       O    X    O       it is more vertical.  There are eight
%       |         |       divisions, but for the non-maximum suppression
%    (1)|         |(4)    we are only worried about 4 of them since we
%       O----O----O       use symmetric points about the center pixel.
%        (2)   (3)


[m,n,l] = size(mag);

% Find the indices of all points whose gradient (specified by the
% vector (ix,iy)) is going in the direction we're looking at.

I = {-iz, -iy, -ix, 0, ix, iy, iz};

dirIdx1 = [1,2,3; 1,3,2; 2,1,3; 2,3,1; 3,1,2; 3,2,1]; 
dirIdx2 = repmat([-1,1,1], [6,1]) .* dirIdx1; 
dirIdx3 = repmat([1,-1,1], [6,1]) .* dirIdx1;
dirIdx4 = repmat([1,1,-1], [6,1]) .* dirIdx1;
dirIdx = [dirIdx1; dirIdx2; dirIdx3; dirIdx4]; 
dirIdx = [zeros(24,1), dirIdx]; dirIdx = dirIdx + 4;

di = dirIdx(direction, :); 

idx = find( ((I{di(1)} <= I{di(2)}) & (I{di(2)} <= I{di(3)})) | ...
    ((I{di(1)} >= I{di(2)}) & (I{di(2)} >= I{di(3)})) );

%{
switch direction
    case 1
        idx = find((iy<=0 & ix>-iy) | (iy>=0 & ix<-iy));
    case 2
        idx = find((ix>0 & -iy>=ix) | (ix<0 & -iy<=ix));
    case 3
        idx = find((ix<=0 & ix>iy)  | (ix>=0 & ix<iy));
    case 4
        idx = find((iy<0 & ix<=iy)  | (iy>0 & ix>=iy));
end
    %}

% Exclude the exterior pixels
if ~isempty(idx)
    v = mod(idx, m); u = mod(idx, m*n);
    extIdx = (v==1 | v==0 | u<m | (u>(n-1)*m) | idx<=m*n | idx>m*n*(l-1));
    idx(extIdx) = [];
end

ixv = ix(idx);
iyv = iy(idx);
izv = iz(idx);
gradmag = mag(idx);

% Do the linear interpolations for the interior pixels

%{
switch direction
    case 1
        d = abs(iyv./ixv);
        gradmag1 = mag(idx+m).*(1-d) + mag(idx+m-1).*d;
        gradmag2 = mag(idx-m).*(1-d) + mag(idx-m+1).*d;
    case 2
        d = abs(ixv./iyv);
        gradmag1 = mag(idx-1).*(1-d) + mag(idx+m-1).*d;
        gradmag2 = mag(idx+1).*(1-d) + mag(idx-m+1).*d;
    case 3
        d = abs(ixv./iyv);
        gradmag1 = mag(idx-1).*(1-d) + mag(idx-m-1).*d;
        gradmag2 = mag(idx+1).*(1-d) + mag(idx+m+1).*d;
    case 4
        d = abs(iyv./ixv);
        gradmag1 = mag(idx-m).*(1-d) + mag(idx-m-1).*d;
        gradmag2 = mag(idx+m).*(1-d) + mag(idx+m+1).*d;
end
%}
idxLocalMax = idx(gradmag>=gradmag1 & gradmag>=gradmag2);
end

end