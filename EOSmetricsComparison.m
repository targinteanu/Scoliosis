%% load patients 
clear;

patients_loaded = [];

qst = questdlg('Add a new directory?', 'Directory Selection', 'Yes', 'No', 'No');

while strcmp(qst, 'Yes')

    [fn, fp] = uigetfile('*.mat', 'Select Scan Directory File'); 
    load([fp, '\', fn]); 

    patients_avail = [];
    for p = patient_list 
        fnp = [base_fp, num2str(p), img_fp];
        if exist([fnp, 'patient',num2str(p),' EOSoutline data.mat'], 'file')
            patients_avail = [patients_avail, p];
        end
    end
    
    blankPatient.splfilt = []; blankPatient.femheadsScl = [];
    PL = repmat(blankPatient, size(patients_avail));
    for i = 1:length(patients_avail)
        p = patients_avail(i);
        load([base_fp, num2str(p), img_fp, 'patient',num2str(p),' filtered data.mat']);
        load([base_fp, num2str(p), img_fp, 'patient',num2str(p),' EOSoutline data.mat']);
        PL(i).splfilt = splfilt;
        PL(i).femheadsScl = femheadsScl;
    end
    
    patients_loaded = [patients_loaded, PL];
    
qst = questdlg('Add a new directory?', 'Directory Selection', 'Yes', 'No', 'No');
end


%% collect metrics for all patients 
q = 10;
ntot = length(patients_loaded);

    z = zeros(ntot, 1); % empty 
    CobbCor = z; % Max Coronal Cobb Angle 
    apex = z; % apex location of largest Coronal curve (by Cobb angle)
    neut = z; % neutral location closest to 'apex'
    LL = z; % Lumbar Lordosis Cobb Angle 
    SS = z; % Sacral Slope 
    PT = z; % Pelvic Tilt
    Wri = z; % Writhe (numeric integration) 
    K = z; % Max Curvature (Lewiner) 
    KL = z; % location of 'K'
    T = z; % Max Torsion (Lewiner) 
    TL = z; % location of 'T'
    SVA = z; % Sagittal Vertical Alignment 
    CVA = z; % Coronal Vertical Alignment
    VA3D = z; % 3D Vertical Alignment norm
    nCor = z; % number of Coronal curves 
    nSag = z; % number of Sagittal curves 
    
for i = 1:ntot
    splfilt = patients_loaded(i).splfilt;
    femheadsScl = patients_loaded(i).femheadsScl;
    
    [SVA(i), CVA(i), VA3D(i)] = SCVA(splfilt); 
    [SS(i), PT(i)] = getPelvicParams(splfilt, femheadsScl);
    
    XYZH = PlumblineDistance(splfilt, 1);
    [idxMin, idxMax] = localMinMax(XYZH);
    thetas = cobbAngleMinMax(XYZH, idxMin, idxMax);
    thetasSag = thetas(:,2);
    LL(i) = thetasSag(end); % mostly correct BUT some patients show only one sag angle ????
    nSag(i) = length(thetasSag);
    
    XYZH = PlumblineDistance(splfilt, 2);
    [idxMin, idxMax] = localMinMax(XYZH);
    thetas = cobbAngleMinMax(XYZH, idxMin, idxMax);
    thetasCor = thetas(:,3);
    [CobbCor(i), idxMaxTheta] = max(thetasCor);
    nCor(i) = length(thetasCor);
    
    idxMaxApex = idxMax(idxMaxTheta);
    [~,idxClosestNeutral] = min(abs(idxMin - idxMaxApex));
    idxClosestNeutral = idxMin(idxClosestNeutral);
    apex(i) = idxMaxApex/size(XYZH,1);
    if isempty(idxMin)
        neut(i) = 0;
    else
        neut(i) = idxClosestNeutral/size(XYZH,1);
    end
    
    % check fix sign convention ?????
    
    Wri(i) = getWrithe(splfilt); 
    
    [d1, d2, d3, tau, kappa] = LewinerQuantity(splfilt, q);
    [K(i), iK] = max(kappa); iK = iK+q; % only positive values possible 
    [~, iT] = max(abs(tau)); T(i) = tau(iT); iT = iT+q;
    
    KL(i) = (iK)/(2*q + length(kappa));
    TL(i) = (iT)/(2*q + length(tau));

    i/ntot
end
%% construct table and view correlations
VarTable = table(CobbCor, LL, SS, PT, SVA, CVA, VA3D, K, T, Wri, apex, neut, KL, TL, nCor, nSag);
VarTable.Properties.DimensionNames = {'Patient', 'Variables'};
VarTable.PI = VarTable.SS + VarTable.PT; % Pelvic Incidence 
VarTable.absWri = abs(VarTable.Wri); 
VarTable.absT = abs(VarTable.T);
VarTable.absSVA = abs(VarTable.SVA); 
VarTable.absCVA = abs(VarTable.CVA);
VarTable.ODI = 0.089 * VarTable.SVA + 0.253 * (VarTable.PI - VarTable.LL) + 25.332; % ODI predictor (Schwab et al)
varnames = VarTable.Properties.VariableNames;

[R, blankPatient] = corr(VarTable.Variables);
[r,c] = find(blankPatient <= .05);
figure; heatmap(varnames, varnames, R); title('Correlation');
figure; heatmap(varnames, varnames, blankPatient); title('p-value');

%% split into surgical rates groups 
SVAcutoff = 47; % mm
PTcutoff = 22; % deg
PILLcutoff = 11; % deg

gSVA = VarTable.SVA < SVAcutoff; 
gPT = VarTable.PT < PTcutoff; 
gPILL = VarTable.PI - VarTable.LL < PILLcutoff; 
gSVA_PT = gSVA&gPT; gSVA_PILL = gSVA&gPILL; gPT_PILL = gPT&gPILL; gSVA_PT_PILL = gSVA&gPT&gPILL;

groups = {gSVA, gPT, gPILL, gSVA_PT, gSVA_PILL, gPT_PILL, gSVA_PT_PILL}; 
groupnames = {'SVA', 'PT', 'PI-LL', 'SVA & PT', 'SVA & PI-LL', 'PT & PI-LL', 'all'};
pVarDiff = zeros(length(groups), length(varnames));

for gi = 1:length(groups)
    g = groups{gi};
    vars1 = VarTable{g,:}; vars2 = VarTable{~g,:};
    for vi = 1:length(varnames)
        [~,pVarDiff(gi,vi)] = ttest2(vars1(:,vi), vars2(:,vi), 'Vartype', 'unequal');
    end
end

showComparison = sum((pVarDiff < .1)) > 0;
pVarSel = pVarDiff(:, showComparison); varnamesSel = varnames(showComparison);
figure; heatmap(varnamesSel, groupnames, pVarSel);
%figure; heatmap(varnames, groupnames, pVarDiff);

%% More details for each coronal curve 
ncurves = unique(VarTable.n);
for ncurve = ncurves'
    pp = patients_loaded(VarTable.n'==ncurve);
    cobbs = zeros(length(pp), ncurve);
    subWr = zeros(size(cobbs));
    wr = zeros(length(pp),1);
    for i = 1:length(pp)
        splfilt = pp(i).splfilt;
        
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
    
    subVarTable = table(cobbs, wr, subWr); 
    subVarTable.absWr = abs(subVarTable.wr);
    subVarTable.absSubWr = abs(subVarTable.subWr);
    vars = subVarTable.Variables; 
    varnames = subVarTable.Properties.VariableNames;
    varnames2 = cell(1, size(vars,2));
    v2 = 1;
    for v = 1:length(varnames)
        vn = size(subVarTable(:,v).Variables,2);
        for v3 = 1:vn
            varnames2{v2} = [varnames{v},' ',num2str(v3)];
            v2 = v2+1;
        end
    end
    
    [R, blankPatient] = corr(vars);
    figure; heatmap(varnames2, varnames2, R); title([num2str(ncurve),'Curves Correlation']); 
    figure; heatmap(varnames2, varnames2, blankPatient); title([num2str(ncurve),'Curves p-value']);
    
end        

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



function [SS, PT] = getPelvicParams(splfilt, femheadsScl)
projuv = @(u,v) ((u*v')/(v*v'))*v;
fhXYZ = femheadsScl;
femaxis = fhXYZ(:,2) - fhXYZ(:,1);
R = splfilt; R = R(:,[1,3]);

% Sacral Slope
dR = diff(R); 
Ri = .5 * ( R(2:end,:) + R(1:(end-1),:) ); dRi = diff(Ri);
ds = sum(dR.^2, 2).^.5;
dsi = sum(dRi.^2, 2).^.5;
dR = dR./ds;
ddR = diff(dR);
ddR = ddR./dsi;
n_sacralPlate = ddR(end,:); n_sacralPlate = -n_sacralPlate/norm(n_sacralPlate);
n2 = [-1,0];
SS = acos(n_sacralPlate * n2') * 180/pi;

% Pelvic Tilt
a1 = fhXYZ(:,1)'; a12 = femaxis'; b = splfilt(end,:) - a1; 
b_parl = projuv(b, a12); b_perp = b - b_parl; b_perp = -[b_perp(1), b_perp(3)];
n3 = [0,1]; 
PT = acos(b_perp/norm(b_perp) * n3') * 180/pi;
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