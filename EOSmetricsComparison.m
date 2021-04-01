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

ntot = length(patients_avail);
vars = zeros(ntot, 5);
% vars key: 
%   1 - max coronal cobb angle 
%   2 - writhe 
%   3 - max curvature 
%   4 - max torsion 
%   5 - # of coronal curves 
for i = 1:ntot   
    p = patients_avail(i);
    load([base_fp, num2str(p), img_fp, 'patient',num2str(p),' EOSoutline data.mat']);
    load([base_fp, num2str(p), img_fp, 'patient',num2str(p),' filtered data.mat']);
    
    XYZH = PlumblineDistance(splfilt, 2);
    [idxMin, idxMax] = localMinMax(XYZH);
    thetas = cobbAngleMinMax(XYZH, idxMin, idxMax);
    thetasCor = thetas(:,3);
    vars(i,1) = max(thetasCor);
    vars(i,5) = length(thetasCor);
    
    vars(i,2) = getWrithe(splfilt);
    
    [d1, d2, d3, tau, kappa] = LewinerQuantity(splfilt, 10);
    vars(i,3) = max(kappa); 
    vars(i,4) = max(tau);

    i/ntot
end

[R, P] = corr(vars);
figure('Position', [100 100 1100 400]); 
subplot(121); heatmap(R); title('Correlation'); 
subplot(122); heatmap(P); title('p-value');

figure; 
plot(vars(:,1), abs(vars(:,2)), '.'); hold on; grid on;
text(vars(:,1), abs(vars(:,2)), ...
    arrayfun(@(n) num2str(n), vars(:,5), 'UniformOutput', false));
xlabel('Max Coronal Cobb Angle'); ylabel('|Writhe|');

figure; 
plot(vars(:,3), (vars(:,2)), '.'); hold on; grid on;
text(vars(:,3), (vars(:,2)), ...
    arrayfun(@(n) num2str(n), vars(:,5), 'UniformOutput', false));
xlabel('Max Curvature'); ylabel('Writhe');

figure; 
plot(vars(:,4), (vars(:,2)), '.'); hold on; grid on;
text(vars(:,4), (vars(:,2)), ...
    arrayfun(@(n) num2str(n), vars(:,5), 'UniformOutput', false));
xlabel('Max Torsion'); ylabel('Writhe');

%%

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

function thetas = cobbAngleMinMax(SplSclHt, lc_min, lc_max)
x = SplSclHt(:,1); y = SplSclHt(:,2); z = SplSclHt(:,3); R = [x,y,z];
neutIdx = [1; lc_min; length(z)]; apIdx = lc_max;
boundIdx = zeros(length(apIdx),2);
for i = 1:length(apIdx)
    d = z(neutIdx)-z(apIdx(i));
    boundIdx(i,1) = neutIdx(minSgn(d, -1));
    boundIdx(i,2) = neutIdx(minSgn(d, 1));
end
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

if (nargin < 4)|(varToShow==3)
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