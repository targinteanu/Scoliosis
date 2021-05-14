%% load patients 
%%{
clear;
[fn, fp] = uigetfile('*.mat', 'Select Scan Directory File'); 
load([fp, '\', fn]); 

patients_avail = [];
for p = patient_list 
    fnp = [base_fp, num2str(p), img_fp];
    if exist([fnp, 'patient',num2str(p),' EOSoutline data.mat'], 'file')
        patients_avail = [patients_avail, p];
    end
end
%}

%% collect metrics for all patients 
%%{
q = 10;
ntot = length(patients_avail);
varnames = {'Cobb', 'Wr', '|Wr|', 'K', '|K|', 'T', '|T|', 'SVA', 'CVA', '3VA', 'n', 'apex', 'neut', 'Lcurv', 'Ltors'};
vars = zeros(ntot, length(varnames));
% vars key: 
%   1 - max coronal cobb angle 
%   2 - writhe
%   3 - |Writhe|
%   4 - max curvature 
%   5 - |max curvature|
%   6 - max torsion 
%   7 - |max torsion|
%   8 - SVA (Sag Vertical Alignment)
%   9 - CVA (Cor Vertical Alignment)
%   10 - 3D Vertical Alignment Norm
%   11 - # of coronal curves 
%   12 - largest apex loc
%   13 - neutral loc
%   14 - max curvature loc
%   15 - max torsion loc
for i = 1:ntot
    p = patients_avail(i);
    load([base_fp, num2str(p), img_fp, 'patient',num2str(p),' EOSoutline data.mat']);
    load([base_fp, num2str(p), img_fp, 'patient',num2str(p),' filtered data.mat']);
    
    [SVA, CVA, VA3D] = SCVA(splfilt);
    vars(i,8) = SVA; vars(i,9) = CVA; vars(i,10) = VA3D;
    
    XYZH = PlumblineDistance(splfilt, 2);
    [idxMin, idxMax] = localMinMax(XYZH);
    thetas = cobbAngleMinMax(XYZH, idxMin, idxMax);
    thetasCor = thetas(:,3);
    [vars(i,1), idxMaxTheta] = max(thetasCor);
    vars(i,11) = length(thetasCor);
    
    idxMaxApex = idxMax(idxMaxTheta);
    [~,idxClosestNeutral] = min(abs(idxMin - idxMaxApex));
    idxClosestNeutral = idxMin(idxClosestNeutral);
    vars(i,12) = idxMaxApex/size(XYZH,1);
    if isempty(idxMin)
        vars(i,13) = 0;
    else
        vars(i,13) = idxClosestNeutral/size(XYZH,1);
    end
    
    vars(i,2) = (getWrithe(splfilt));
    vars(i,3) = abs(vars(i,2));
    
    [d1, d2, d3, tau, kappa] = LewinerQuantity(splfilt, q);
    [vars(i,4), iK] = max(kappa); vars(i,5) = abs(vars(i,4));
    [vars(i,6), iT] = max(tau); vars(i,7) = abs(vars(i,6));
    
    vars(i,14) = (iK+q)/(2*q + length(kappa));
    vars(i,15) = (iT+q)/(2*q + length(tau));

    i/ntot
end

[R, P] = corr(vars);
[r,c] = find(P <= .05);
%figure('Position', [20 100 1450 600]); 
%subplot(121); heatmap(varnames, varnames, R); title('Correlation'); 
%subplot(122); heatmap(varnames, varnames, P); title('p-value');
figure; heatmap(varnames, varnames, R); title('Correlation');
figure; heatmap(varnames, varnames, P); title('p-value');

figure; 
plot(vars(:,1), abs(vars(:,2)), '.'); hold on; grid on;
text(vars(:,1), abs(vars(:,2)), ...
    arrayfun(@(n) num2str(n), vars(:,5), 'UniformOutput', false));
xlabel('Max Coronal Cobb Angle'); ylabel('|Writhe|');

%}

%% more details for each curve 
%{
ncurves = unique(vars(:,5));
for ncurve = ncurves'
    pp = patients_avail(vars(:,5)'==ncurve);
    cobbs = zeros(length(pp), ncurve);
    subWr = zeros(size(cobbs));
    wr = zeros(length(pp),1);
    for i = 1:length(pp)
        p = pp(i);
        load([base_fp, num2str(p), img_fp, 'patient',num2str(p),' EOSoutline data.mat']);
        load([base_fp, num2str(p), img_fp, 'patient',num2str(p),' filtered data.mat']);
        
        XYZH = PlumblineDistance(splfilt, 2);
        %%{
        [idxMin, idxMax] = localMinMax(XYZH);
        thetas = cobbAngleMinMax(XYZH, idxMin, idxMax);
        [thetas, thetaIdx] = sort(thetas(:,3)');
        cobbs(i,:) = thetas;
        
        boundIdx = getBoundIdx(idxMin, idxMax, XYZH(:,3));
        sw = zeros(size(boundIdx,1),1);
        for j = 1:size(boundIdx,1)
            Rng = (max(boundIdx(j,1),2)):(min(boundIdx(j,2),size(XYZH,1)));
            sw(j) = getWrithe(XYZH(:,1:3), Rng);
        end
        subWr(i,:) = sw(thetaIdx);
        wr(i) = getWrithe(XYZH(:,1:3));
        %%}
        
        i/length(pp)
    end
    figure; 
    subplot(211); plot(cobbs); title(num2str(ncurve)); grid on;
    subplot(212); plot(subWr); title(num2str(ncurve)); grid on;
    
    [R, P] = corr([abs(wr), cobbs, abs(subWr)]);
    figure('Position', [20 100 1450 600]); 
    subplot(121); heatmap(R); title('Correlation'); 
    subplot(122); heatmap(P); title('p-value');
    
end        

%}

%% functions 

function [theta, n1, n2] = numericCobbAngle(R, t1, t2)
dR = diff(R); 
Ri = .5 * ( R(2:end,:) + R(1:(end-1),:) ); dRi = diff(Ri);
ds = sum(dR.^2, 2).^.5;
dsi = sum(dRi.^2, 2).^.5;
dR = dR./ds;
ddR = diff(dR);
ddR = ddR./dsi;
xy1 = R(t1,:); xy2 = R(t2,:);
t1 = t1-1; t2 = t2-1;
t1 = min(t1, size(ddR,1)); t1 = max(t1,1);
t2 = min(t2, size(ddR,1)); t2 = max(t2,1);
n1 = ddR(t1,:); n1 = n1/norm(n1);
n2 = ddR(t2,:); n2 = n2/norm(n2);
if length(xy1) < 3
    [~,sgn] = intersect2d(xy1', n1', xy2', n2'); sgn=sign(sgn); 
    n1=sgn(1)*n1; n2=sgn(2)*n2;
end
theta = acos(n1 * n2') * 180/pi;
%theta = min(theta, 180-theta); % come up with a better way of picking the right angle 
end

function [xy3,c] = intersect2d(xy1, a, xy2, b)
c = ([a, b])^-1 * (xy1 - xy2); c(1) = -c(1);
xy3 = c(2)*b + xy2;
end

function [thta3, thtaSag, thtaCor] = getCobbAngle(R, t1, t2)
% x - sag
xy = R(:,[1,3]); 
[thtaSag, a, b] = numericCobbAngle(xy, t1, t2);
% y - cor
xy = R(:,[2,3]); 
[thtaCor, a, b] = numericCobbAngle(xy, t1, t2);
% 3d 
[thta3, n1, n2] = numericCobbAngle(R, t1, t2);
end

function idx = minSgn(vals, sgn)
origIdx = 1:length(vals);
vals = vals*sgn; 
sel = vals > 0;
origIdx = origIdx(sel); vals = vals(sel);
[~,i] = min(vals); idx = origIdx(i);
end

function boundIdx = getBoundIdx(lc_min, lc_max, z)
neutIdx = [1; lc_min; length(z)]; apIdx = lc_max;
boundIdx = zeros(length(apIdx),2);
for i = 1:length(apIdx)
    d = z(neutIdx)-z(apIdx(i));
    boundIdx(i,1) = neutIdx(minSgn(d, -1));
    boundIdx(i,2) = neutIdx(minSgn(d, 1));
end
end

function thetas = cobbAngleMinMax(SplSclHt, lc_min, lc_max)
x = SplSclHt(:,1); y = SplSclHt(:,2); z = SplSclHt(:,3); R = [x,y,z];
boundIdx = getBoundIdx(lc_min, lc_max, z);
thetas = getCobbAngles(R, boundIdx);
end

function thetas = getCobbAngles(R, boundIdx)
thetas = zeros(size(boundIdx,1), 3);
for i = 1:size(boundIdx,1)
    t1 = boundIdx(i,1); t2 = boundIdx(i,2);
    [thetas(i,1), thetas(i,2), thetas(i,3)] = getCobbAngle(R, t1, t2);
end
end



function [lc_min, lc_max] = localMinMax(SplSclHt)
x = SplSclHt(:,1); y = SplSclHt(:,2); z = SplSclHt(:,3); h = SplSclHt(:,4);
[pk_min,lc_min,pw_min,pp_min] = findpeaks(-h); 
[pk_max,lc_max,pw_max,pp_max] = findpeaks(h);
end



function SplSclHt = PlumblineDistance(R, varToShow)
projuv = @(u,v) ((u*v')/(v*v'))*v;

if (nargin < 2)|(varToShow==3)
    R0 = R;
else
    R0 = R(:,[varToShow,3]);
end

plumb = R0([1,end],:); 
vPlumb = diff(plumb); 
R0 = R0 - plumb(1,:);

dist = arrayfun(@(i) norm( R0(i,:) - projuv(R0(i,:), vPlumb) ), 1:size(R0,1))';
SplSclHt = [R,dist];
end



function [d1, d2, d3, tau, kappa] = LewinerQuantity(p, q)
    vertebrae = (1+q):(size(p,1)-q); 
    tau = zeros(size(vertebrae)); % local torsion at each point on the spine 
    kappa = zeros(size(vertebrae)); % local curvature at each point on the spine 
    d2 = zeros(size(vertebrae)); d1 = d2; d3 = d2;
    for vertebra = vertebrae
        [tau(vertebra-q), d,dd,ddd, kappa(vertebra-q)] = ...
            lewinerTorsion(p, vertebra, q);
        d2(vertebra-q) = norm(dd); d1(vertebra-q) = norm(d); d3(vertebra-q) = norm(ddd);
    end
end



function [SVA, CVA, VA3D] = SCVA(XYZ)

p1 = XYZ(1,:); p2 = XYZ(end,:); 
p3 = [p1(1), p1(2), p2(3)];

% x - sag
SVA = -(p1(1)-p2(1));

% y - cor
CVA = (p1(2)-p2(2));

% 3D
VA3D = norm(p3-p2);

end