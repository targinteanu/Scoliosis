fn = 'C:\Users\Toren\Desktop\scoliosis\EOS patient 36\selected DICOM\sag';
imsag = dicomread(fn);
imsag = double(imsag); imsag = imsag - min(imsag(:)); imsag = imsag/max(imsag(:));
ifo = dicominfo(fn);
xscl = ifo.PixelSpacing(1); zscl = ifo.PixelSpacing(2);

figure; imshow(imsag); hold on;

%[outln, vertx, verty] = roipoly; 
load('EOStestROI.mat');
[y1,idx1, y2,idx2] = minmax2(@(v) max(v), verty);
x1 = vertx(idx1); x2 = vertx(idx2);

% img is a double img [0, 1]; outln is a BW outline

% filtering 
im = imreconstruct(double(~outln), imsag); 
imgFiltered = imsag - im; 
imgFiltered = imgFiltered .* outln;
imgFiltered = imgFiltered/max(imgFiltered(:));

% spline fitting 
[r, c] = find(imgFiltered);
[splineObj, gof] = fit(r*zscl, c*xscl, 'smoothingspline', 'Weights', imgFiltered(find(imgFiltered(:))));
z = min(r):max(r); 
z = z*zscl;
splineSample = [z/zscl; ppval(z, splineObj.p)/xscl]';
% display results on current axes 
% splineSample is in pixels, but splineObj is in mm

plot([x1,x2],[y1,y2],'o', 'LineWidth',2);
x1=x1*xscl; x2=x2*xscl; y1=y1*zscl; y2=y2*zscl;
m = (y2-y1)/(x2-x1); k = y1-m*x1;

imgFiltered = imgFiltered.*(imgFiltered > .2);
[r, c] = find(imgFiltered);
N = 10;
coeff1 = WLS(c*xscl, r*zscl, imgFiltered(find(imgFiltered(:))), N);
coeff2 = WLStangent(c*xscl, r*zscl, imgFiltered(find(imgFiltered(:))), N, -1/m, k, 7.3672*xscl);
coeff2fun = @(coeff, x) arrayfun(@(i) sum(coeff.*(x(i).^(0:(length(coeff)-1)))), 1:length(x));

visboundaries(outln); 
plot(splineSample(:,2), splineSample(:,1), 'b');
plot(coeff2fun(coeff1',z)/xscl, z/zscl, 'r');
plot(coeff2fun(coeff2',z)/xscl, z/zscl, 'm');

% helper functions --------------------------------------------------

function [v1,idx1, v2,idx2] = minmax2(fun, vals)
    vals = unique(vals);
    [v1, idx1] = fun(vals);
    vals2 = vals([1:(idx1-1), (idx1+1):end]);
    v2 = fun(vals2);
    idx2 = find(vals == v2);
end

function beta = WLS(x,y,w,N)
    X = cell2mat(arrayfun(@(xi) xi.^(0:N), x, 'UniformOutput', false));
    beta = ((X'*diag(w)*X)^-1)*X'*diag(w)*y;
end

function beta = WLStangent(x,y,w,N, m0,c0,xstart)
    X = cell2mat(arrayfun(@(xi) xi.^([0,2:N]), x, 'UniformOutput', false));
    b = ((X'*diag(w)*X)^-1)*X'*diag(w)*y;
    
    i = 2:N;
    eqn = @(x0) b(1) - c0 + x0.*(m0 + 1/m0) + sum((i-1).*b(2:end).*(x0.^i));
    %options = optimoptions('fsolve','Display','none','PlotFcn',@optimplotfirstorderopt);
    solx = fsolve(eqn, xstart);
    b1 = m0 - sum(i.*b(2:end)'.*(solx.^(i-1)));
    beta = [b(1); b1; b(2:end)];
end