function [polyFilt, polyFiltVal, polyFS, polyFSVal, polyVal, X] = ...
    piecepolyfilter(wcutoff, xBreaks, yCoefs, xRng, xN)
% Applies rectangular filters to piecewise polynomial e.g. spline. 
% Input: 
%   wcutoff: frequency cutoff, in units of 1/[x units] 
%   xBreaks: bounds of x intervals as (L+1)-by-1 vector 
%   yCoefs: order-m polynomial coefficients corresponding to each pair of 
%           bounds as rows of an L-by-m matrix, going from zeroth-order
%           coefficient (first column) to mth-order coefficient (last
%           column). 
%   xRng: range of X values to evaluate 
%   xN: # of points per interval. If unspecified, default = 5
% Output: 
%   result of rectangular filter using convolution in [y units]:
%     polyFilt: L-by-1 cell array of functions corresponding to each pair
%               of bounds 
%     polyFiltVal: vector evaluated as function of X 
%   result of curtailed Fourier Series in [y units]
%     polyFS: L-by-1 cell array of functions corresponding to each pair
%             of bounds
%     polyFSVal: vector evaluated as function of X 
%   polyVal: unfiltered piecewise polynomial vector evaluated as function 
%            of X in [y units]
%   X: vector independent variable of evaluated functions in [x units]

if nargin < 5
    xN = 5; % default to 5 x points per interval
end

% Select Range 
inRng = (xBreaks > xRng(1))&(xBreaks < xRng(2));
xBreaks = [xRng(1), xBreaks(inRng), xRng(2)];
inRng = (inRng(1:(end-1)))&(inRng(2:end));
yCoefs = yCoefs(inRng, :);

ply_helper = @(p, x) p(1) + p(2)*x + p(3)*x.^2 + p(4)*x.^3;

ply = @(p, x, b) ply_helper(p, x) .* (x>=b(1) & x<=b(2));

filtply_helper = @(w, p, x, t) (1/(pi*w^3))*(...
    w^3 * sinint(w*(t-x)) .* (ply_helper(p, x)) -...
    cos(w*(t-x)) .* (w^2 * (p(2) + p(3)*(t+x)) + p(4)*(t*t*w*w + t*w*w*x + w*w.*x.*x - 2)) +...
    w * (p(3) + p(4)*(2*t + x)) .* sin(w*(t-x)));

filtply = @(w, p, x, b) filtply_helper(w, p, x, b(2)) - filtply_helper(w, p, x, b(1));

FScoef_helper = @(p, x, n, T) ((exp(-2*j*pi*n*x/T))/(8 * pi^4 * n^4)) .* (...
    4*j*p(1) * pi^3 * n^3 + 2*p(2) * pi^2 * n^2 *(T + 2*j*pi*n*x) + ...
    4*j*p(3) * pi^3 * n^3 * x.^2 + 4*p(3)*T*x * pi^2 * n^2 -...
    2*j*pi*p(3)*n * T^2 + 4*j*p(4) * pi^3 * n^3 * x.^3 +...
    6*p(4)*T * pi^2 * n^2 * x.^2 - 6*j*pi*p(4)*n*x * T^2 - 3*p(4) * T^3);

FScoef0_helper = @(p, x, T) (p(1)*x + (p(2)*x.^2)/2 + (p(3)*x.^3)/3 + (p(4)*x.^4)/4)/T;

FScoef = @(p, n, a, b) FScoef_helper(p, b, n, b-a) - FScoef_helper(p, a, n, b-a);
FScoef0 = @(p, a, b) FScoef0_helper(p, b, b-a) - FScoef0_helper(p, a, b-a);

FSn = @(p, x, n, b) exp(j*2*pi*n*x/(b(2)-b(1))) * FScoef(p, n, b(1), b(2));
FS0 = @(p, x, b) FScoef0(p, b(1), b(2));

FSfilt = @(w, p, x, b) sum(cell2mat(arrayfun(@(n) FSn(p,x,n,b)', ...
    [-ceil(w*(b(2)-b(1))/(2*pi)):-1, 1:ceil(w*(b(2)-b(1))/(2*pi))], ...
    'UniformOutput', false))') + FS0(p, x, b);

L = size(yCoefs,1);
polyFilt = cell(L, 1);
polyFS = cell(L, 1);
X = zeros(xN*L, 1);
polyFiltVal = zeros(size(X));
polyFSVal = zeros(size(X));
polyVal = zeros(size(X));

for i = 1:L

    P = yCoefs(i,:);
    B = [xBreaks(i), xBreaks(i+1)];
    
    x = linspace(B(1), B(2), xN);
    PLY = @(x) ply(P, x, B);
    PLYfilt = @(x) filtply(wcutoff, P, x, B);
    PLYfiltFS = @(x) real(FSfilt(wcutoff, P, x, B));
    
    idx = (i-[1,0])*xN + [1,0];
    X(idx(1):idx(2)) = x;
    polyVal(idx(1):idx(2)) = PLY(x); 
    polyFiltVal(idx(1):idx(2)) = PLYfilt(x);
    polyFSVal(idx(1):idx(2)) = PLYfiltFS(x);
    
    polyFilt{i} = PLYfilt; 
    polyFS{i} = PLYfiltFS;

end

end