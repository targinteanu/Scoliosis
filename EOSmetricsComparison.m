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
    
    blankPatient.splfilt = []; blankPatient.femheadsScl = []; blankPatient.splSclBnd = [];
    blankPatient.num = 0;
    PL = repmat(blankPatient, size(patients_avail));
    for i = 1:length(patients_avail)
        p = patients_avail(i);
        load([base_fp, num2str(p), img_fp, 'patient',num2str(p),' filtered data.mat']);
        load([base_fp, num2str(p), img_fp, 'patient',num2str(p),' EOSoutline data.mat']);
        PL(i).splfilt = splfilt;
        PL(i).femheadsScl = femheadsScl;
        PL(i).splSclBnd = splSclBnd;
        PL(i).num = p;
        i/length(patients_avail)
    end
    
    patients_loaded = [patients_loaded, PL];
    
qst = questdlg('Add a new directory?', 'Directory Selection', 'Yes', 'No', 'No');
end


%% collect metrics for all patients 
q = 10;
ntot = length(patients_loaded);

    z = zeros(ntot, 1); % empty 
    CobbCor = z; % Coronal Cobb Angle with max magnitude
    apex = z; % apex location of largest Coronal curve (by Cobb angle)
    neut = z; % neutral location closest to 'apex'
    LL = z; % Lumbar Lordosis Cobb Angle 
    TK = z; % Thoracic Kyphosis Cobb Angle
    SS = z; % Sacral Slope 
    PT = z; % Pelvic Tilt
    Wri = z; % Writhe (numeric integration) 
    K = z; % Max Curvature (Lewiner) 
    KL = z; % location of 'K'
    Ka = z; % Curvature (Lewiner) at 'apex'
    T = z; % Max Torsion (Lewiner) 
    TL = z; % location of 'T'
    Tn = z; % Torsion (Lewiner) at 'neut'
    SVA = z; % Sagittal Vertical Alignment 
    CVA = z; % Coronal Vertical Alignment
    VA3D = z; % 3D Vertical Alignment norm
    nCor = z; % number of Coronal curves 
    nSag = z; % number of Sagittal curves 
    
for i = 1:ntot
    splfilt = patients_loaded(i).splfilt;
    femheadsScl = patients_loaded(i).femheadsScl;
    splSclBnd = patients_loaded(i).splSclBnd;
    
    [SVA(i), CVA(i), VA3D(i)] = SCVA(splfilt); 
    [SS(i), PT(i)] = getPelvicParams(splfilt, femheadsScl);
    
%    XYZH = PlumblineDistance(splfilt, 1);
%    [idxMin, idxMax] = localMinMax(XYZH);
%    thetas = cobbAngleMinMax(XYZH, idxMin, idxMax);
    thetas = getLL(splfilt, splSclBnd);
    thetasSag = thetas(:,2);
    LL(i) = thetasSag(end); 
    TK(i) = thetasSag(1);
    nSag(i) = length(thetasSag);
    
    XYZH = PlumblineDistance(splfilt, 2);
    [idxMin, idxMax] = localMinMax(XYZH);
    thetas = cobbAngleMinMax(XYZH, idxMin, idxMax);
    thetasCor = thetas(:,3);
    [~, idxMaxTheta] = max(abs(thetasCor));
    nCor(i) = length(thetasCor);
    CobbCor(i) = thetasCor(idxMaxTheta);
    
    idxMaxApex = idxMax(idxMaxTheta);
    [~,idxClosestNeutral] = min(abs(idxMin - idxMaxApex));
    idxClosestNeutral = idxMin(idxClosestNeutral);
    apex(i) = dist3frac(splfilt, idxMaxApex);
    if isempty(idxMin)
        neut(i) = 0;
    else
        neut(i) = dist3frac(splfilt, idxClosestNeutral);
    end
        
    Wri(i) = getWrithe(splfilt); 
    
    [d1, d2, d3, tau, kappa] = LewinerQuantity(splfilt, q);
    [K(i), iK] = max(kappa); iK = iK+q; % only positive values possible 
    [~, iT] = max(abs(tau)); T(i) = tau(iT); iT = iT+q;
    
    Ka(i) = kappa(idxMaxApex - q);
    if isempty(idxMin)
        Tn(i) = tau(1);
    else
        Tn(i) = tau(idxClosestNeutral - q); 
    end
    
%    KL(i) = (iK)/(2*q + length(kappa));
%    TL(i) = (iT)/(2*q + length(tau));
    KL(i) = dist3frac(splfilt, iK);
    TL(i) = dist3frac(splfilt, iT);

    i/ntot
end

% construct table 
VarTable = table(CobbCor, LL, TK, SS, PT, SVA, CVA, K, T, Wri, apex, neut, KL, TL, Ka, Tn);
VarTable.Properties.DimensionNames = {'Patient', 'Variables'};
VarTable.PI = VarTable.SS + VarTable.PT; % Pelvic Incidence 
VarTable.absWri = abs(VarTable.Wri); 
VarTable.absT = abs(VarTable.T);
VarTable.absTn = abs(VarTable.Tn);
VarTable.absCobbCor = abs(VarTable.CobbCor);
VarTable.absLL = abs(VarTable.LL);
VarTable.absTK = abs(VarTable.TK);
VarTable.absSVA = abs(VarTable.SVA); 
VarTable.absCVA = abs(VarTable.CVA);
VarTable.ODI = 0.089 * VarTable.SVA + 0.253 * (VarTable.PI - VarTable.LL) + 25.332; % ODI predictor (Schwab et al)
varnames = VarTable.Properties.VariableNames;

%% view table variable ranges 
% histograms 
mtot = size(VarTable,2);
ncol = 5; nrow = ceil((mtot+2)/ncol);
figure;
for i = 1:mtot
    subplot(nrow, ncol, i)
    dta = VarTable.(i);
    histogram(dta,10);
    ttl = VarTable(:,i).Properties.VariableNames{1};
    title([ttl,': ',...
        num2str(mean(dta),3),'\pm',num2str(std(dta),3),...
        ' [',num2str(min(dta),3),',',num2str(max(dta),3),']']);
end
i = i+1;
subplot(nrow, ncol, i)
histogram(nCor); title('nCor');
i = i+1;
subplot(nrow, ncol, i)
histogram(nSag); title('nSag');

% box plots 
figure('Color', 'white'); 
boxplot([-VarTable.CobbCor, VarTable.LL, -VarTable.TK, VarTable.SS, VarTable.PT, VarTable.PI],...
    'Labels', {'(i)Cobb', '(ii)L.L.', '(iii)T.K.', '(iv)S.S.', '(v)P.T.', '(vi)P.I.'},...
    'PlotStyle','compact');
ylabel('Angle (degrees)'); title('A'); grid on;

figure('Color', 'white');
boxplot([VarTable.SVA, -VarTable.CVA], 'Labels', {'(i)S.V.A.', '(ii)C.V.A.'},...
    'PlotStyle','compact');
ylabel('Alignment (millimeters)'); title('B'); grid on;

figure('Color', 'white');
subplot(141);
boxplot(VarTable.Wri, 'Labels', {'Writhe'},'PlotStyle','compact');
ylabel('Writhe'); grid on; title('C');
subplot(142);
boxplot([VarTable.K, VarTable.Ka], 'Labels', ...
    {'(i)Cur.Max.', '(ii)Cur.Ap.'},'PlotStyle','compact');
ylabel('[millimeters]^{-1}'); grid on;
subplot(143);
boxplot(VarTable.T, 'Labels', {'(iii)Tor.Max.'},'PlotStyle','compact');
ylabel('[millimeters]^{-1}'); grid on; title('D');
subplot(144);
boxplot(VarTable.Tn, 'Labels', {'(iv)Tor.Nt.'},'PlotStyle','compact');
ylabel('[millimeters]^{-1}'); grid on;

%% view correlations
% matrix 
[R, blankPatient] = corr(VarTable.Variables);
[r,c] = find(blankPatient <= .05);
figure; heatmap(varnames, varnames, R); title('Correlation');
figure; heatmap(varnames, varnames, blankPatient); title('p-value');

% featured linear relationships 
lw = 1;
figure('Color', 'white', 'Position', [10, 100, 1450, 500]);
subplot(131);
plotwithfit(-VarTable.CobbCor, VarTable.Wri, lw); grid on;
title('A: Cobb Angle and Writhe'); xlabel('Cobb Angle (degrees)'); ylabel('Writhe');
subplot(132);
plotwithfit(-VarTable.CobbCor, -VarTable.CVA, lw); grid on;
title('B: Cobb Angle and C.V.A.'); xlabel('Cobb Angle (degrees)'); ylabel('C.V.A. (millimeters)');
subplot(133); 
plotwithfit(VarTable.absCVA, VarTable.absWri, lw); grid on;
title('C: C.V.A. and Writhe'); xlabel('C.V.A. Magnitude (mm)'); ylabel('Writhe Magnitude');

figure('Color', 'white', 'Position', [10, 100, 1450, 500]);
subplot(131);
plotwithfit(VarTable.Ka, VarTable.absCobbCor, lw); grid on; 
title('A: Apical Curvature and Cobb Angle'); 
xlabel('Curvature (mm^{-1})'); ylabel('Cobb Angle Magnitude (deg.)');
subplot(132);
plotwithfit(VarTable.Ka, VarTable.absTn, lw); grid on; 
title('B: Apical Curvature and Neutral Torsion'); 
xlabel('Curvature (mm^{-1})'); ylabel('Torsion Magnitude (mm^{-1})');
subplot(133);
plotwithfit(VarTable.Ka, VarTable.absT, lw); grid on;
title('C: Apical Curvature and Max Torsion');
xlabel('Curvature (mm^{-1})'); ylabel('Torsion Magnitude (mm^{-1})');

figure('Color', 'white', 'Position', [10, 100, 1450, 500]);
subplot(121);
[~,gof1] = plotwithfit(VarTable.K, VarTable.LL, lw, 'sb', '--b', false); grid on; hold on;
[~,gof2] = plotwithfit(VarTable.K, VarTable.SS, lw, 'or', '-.r', false); 
R2LL = gof1.rsquare; R2SS = gof2.rsquare;
title('A: Max Curvature, Lumbar Lordosis, and Sacral Slope'); 
xlabel('Curvature (mm^{-1})'); ylabel('Angle (deg.)');
legend('Lumbar Lordosis', ['L.L. fit (R^2=',num2str(R2LL),')'], ...
       'Sacral Slope',    ['S.S. fit (R^2=',num2str(R2SS),')'], ...
       'Location', 'southeast');
subplot(122); 
plotwithfit(VarTable.K, VarTable.SVA, lw); grid on;
title('B: Max Curvature and S.V.A.'); 
xlabel('Curvature (mm^{-1})'); ylabel('S.V.A. (mm)');

figure('Color', 'white', 'Position', [10, 100, 1450, 500]);
subplot(121);
[~,gof1] = plotwithfit(VarTable.absTn, VarTable.LL, lw, 'sb', '--b', false); grid on; hold on;
[~,gof2] = plotwithfit(VarTable.absTn, -VarTable.TK, lw, 'or', '-.r', false); 
R2LL = gof1.rsquare; R2TK = gof2.rsquare;
title('A: Neutral Torsion, Lumbar Lordosis, and Thoracic Kyphosis'); 
xlabel('Torsion Magnitude (mm^{-1})'); ylabel('Angle (deg.)');
legend('Lumbar Lordosis',   ['L.L. fit (R^2=',num2str(R2LL),')'], ...
       'Thoracic Kyphosis', ['T.K. fit (R^2=',num2str(R2TK),')'], ...
       'Location', 'northeast');
subplot(122); 
plotwithfit(VarTable.absTn, VarTable.absSVA, lw); grid on;
title('B: Neutral Torsion and S.V.A.'); 
xlabel('Torsion Magnitude (mm^{-1})'); ylabel('S.V.A. Magnitude (mm)');

figure('Color', 'white', 'Position', [10, 100, 1450, 500]);
subplot(131);
plotwithfit(VarTable.ODI, VarTable.K, lw); grid on; 
title('A: ODI model and Max Curvature'); 
xlabel('ODI prediction'); ylabel('Curvature (mm^{-1})');
subplot(132);
plotwithfit(VarTable.ODI, VarTable.absTn, lw); grid on; 
title('B: ODI model and Neutral Torsion'); 
xlabel('ODI prediction'); ylabel('Torsion Magnitude (mm^{-1})');
subplot(133);
plotwithfit(VarTable.ODI, VarTable.absT, lw); grid on;
title('C: ODI model and C.V.A.');
xlabel('ODI prediction'); ylabel('C.V.A. Magnitude (mm)');

%% split into surgical rates groups 
SVAcutoff = 47; % mm
PTcutoff = 22; % deg
PILLcutoff = 11; % deg

gSVA = VarTable.SVA < SVAcutoff; 
gPT = VarTable.PT < PTcutoff; 
gPILL = VarTable.PI - VarTable.LL < PILLcutoff; 
gSVA_PT = gSVA|gPT; gSVA_PILL = gSVA|gPILL; gPT_PILL = gPT|gPILL; gSVA_PT_PILL = gSVA|gPT|gPILL;

groups = {gSVA, gPT, gPILL, gSVA_PT, gSVA_PILL, gPT_PILL, gSVA_PT_PILL}; 
groupnames = {'SVA', 'PT', 'PI-LL', 'SVA & PT', 'SVA & PI-LL', 'PT & PI-LL', 'all'};
pVarDiff = zeros(length(groups), length(varnames));
nGroups = zeros(length(groups),1);

for gi = 1:length(groups)
    g = groups{gi};
    nGroups(gi) = sum(g);
    vars1 = VarTable{g,:}; vars2 = VarTable{~g,:};
    for vi = 1:length(varnames)
        [~,pVarDiff(gi,vi)] = ttest2(vars1(:,vi), vars2(:,vi), 'Vartype', 'unequal');
    end
end

showComparison = sum((pVarDiff < .1)) > 0;
pVarSel = pVarDiff(:, showComparison); varnamesSel = varnames(showComparison);
figure; heatmap(varnamesSel, groupnames, pVarSel);
nGroups
%figure; heatmap(varnames, groupnames, pVarDiff);

%% box plots
figure('Color', 'white'); 
subplot(221);
boxplot(VarTable.absWri, gSVA,'PlotStyle','compact','Orientation','horizontal',...
    'Labels', {['SVA',hex2dec('2265'),num2str(SVAcutoff),'mm'], ...
    ['SVA<',num2str(SVAcutoff),'mm']}); 
grid on;
ylim([min(VarTable.absWri)-.25*std(VarTable.absWri), max(VarTable.absWri)+.25*std(VarTable.absWri)]); 
ylabel('Writhe'); title('A'); 
subplot(222); 
boxplot(VarTable.absWri, gPILL,'PlotStyle','compact','Orientation','horizontal',...
    'Labels', {['    PI-LL',hex2dec('2265'),num2str(PILLcutoff),hex2dec('00B0')], ...
    ['    PI-LL<',num2str(PILLcutoff),hex2dec('00B0')]}); 
grid on;
ylim([min(VarTable.absWri)-.25*std(VarTable.absWri), max(VarTable.absWri)+.25*std(VarTable.absWri)]); 
ylabel('Writhe'); title('B'); 

subplot(223);
boxplot(VarTable.K, gSVA,'PlotStyle','compact','Orientation','horizontal',...
    'Labels', {['SVA',hex2dec('2265'),num2str(SVAcutoff),'mm'], ...
    ['SVA<',num2str(SVAcutoff),'mm']}); 
grid on;
ylim([min(VarTable.K)-.25*std(VarTable.K), max(VarTable.K)+.25*std(VarTable.K)]); 
ylabel('Curvature'); title('C'); 
subplot(224); 
boxplot(VarTable.K, gPILL,'PlotStyle','compact','Orientation','horizontal',...
    'Labels', {['    PI-LL',hex2dec('2265'),num2str(PILLcutoff),hex2dec('00B0')], ...
    ['    PI-LL<',num2str(PILLcutoff),hex2dec('00B0')]}); 
grid on;
ylim([min(VarTable.K)-.25*std(VarTable.K), max(VarTable.K)+.25*std(VarTable.K)]); 
ylabel('Curvature'); title('D'); 

%% More details for each coronal curve 
ncurves = unique(nCor);
for ncurve = ncurves'
    pp = patients_loaded(nCor'==ncurve);
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

function [fo, gof] = plotwithfit(x, y, lineW, lineSpc, fo_lineSpc, showR2, fitTyp)
if nargin < 7
    fitTyp = 'poly1';
    if nargin < 6
        showR2 = true;
        if nargin < 4
            lineSpc = 'xk';
            fo_lineSpc = '--k';
            if nargin < 3
                lineW = 1;
            end
        else
            fo_lineSpc = ['--', lineSpc(end)];
        end
    end
end

plot(x, y, lineSpc, 'LineWidth', lineW);
[fo, gof] = fit(x, y, fitTyp);
hold on; plot(fo, fo_lineSpc);
legend('off');

if showR2
    txtDst = .1;
    txtX = min(x) + txtDst*(max(x) - min(x));
    txtY = max(y) - txtDst*(max(y) - min(y));
    text(txtX, txtY, ['R^2 = ',num2str(gof.rsquare)], ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end
end


function s = dist3(R, t1, t2)
Ri = .5 * ( R(2:end,:) + R(1:(end-1),:) ); dRi = diff(Ri);
dsi = sum(dRi.^2, 2).^.5;
t1 = max(1, t1-1); t2 = max(1, t2-1);
t1 = min(size(dsi,1), t1); t2 = min(size(dsi,1), t2);
s = sum(dsi(t1:t2));
end

function s = dist3frac(R, t)
s1 = dist3(R, 1, t); s2 = dist3(R, t, size(R,1));
s = s1/(s1+s2);
end


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
sgn = 1;
if length(xy1) < 3
    [xy3,sgn] = intersect2d(xy1', n1', xy2', n2'); sgn=sign(sgn); 
    n1=sgn(1)*n1; n2=sgn(2)*n2;
    if xy3(1) > xy1(1)
        sgn = 1;
    else
        sgn = -1;
    end
end
theta = acos(n1 * n2') * 180/pi; theta = sgn*theta; 
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


function thetas = getLL(splfilt, splSclBnd)
R = splfilt; 
t0 = 1;
t1 = splSclBnd(3) - splSclBnd(2); 
t2 = size(R,1);
boundIdx = [t0 t1; t1 t2];
thetas = getCobbAngles(R, boundIdx); 
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
n_sacralPlate = ddR(end,:); n_sacralPlate = n_sacralPlate/norm(n_sacralPlate);
n2 = [-1,0];
SS = 180; 
while SS > 90
    n_sacralPlate = -n_sacralPlate; 
    SS = acos(n_sacralPlate * n2') * 180/pi;
end

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