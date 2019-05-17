N = 33;
Wr = zeros(1, N); WrA = zeros(1, N); Tor = zeros(1, N); 

noisesz = [1,5,10]; J = length(noisesz);
noiseWr = zeros(J, N); noiseWrA = zeros(J, N); noiseTor = zeros(J, N); 
noiseWrStd = zeros(J, N); noiseWrAStd = zeros(J, N); noiseTorStd = zeros(J, N); 
rng('shuffle');

%%
for idx = 1:N

%%
load('spines_XYZ.mat');

    xyz = spinesXYZ{idx, 2:end};
    x = xyz(1:3:end); 
    y = xyz(2:3:end); 
    z = xyz(3:3:end); 
    cm = [x;y;z]';

%%
writhe_noiseless = levittWrithe(cm); 
writheabs_noiseless = levittWritheAbs(cm); 
    % get torsion 
    q = 2; % 2 vertebrae above, 2 vertebrae below, current vertebra -> 5 points to fit each cubic
    vertebrae = (1+q):(length(x)-q);
    tau = arrayfun(@(v) lewinerTorsion(cm, v, q), vertebrae);
    [torsion_noiseless, tauIdx] = max(abs(tau));

%rng('shuffle');
numtrials = 1000; 
for j = 1:J
writhe = zeros(1, numtrials); 
writheabs = zeros(1, numtrials);
torsion = zeros(1, numtrials);

for trial = 1:numtrials
    noise = noisesz(j)*randn(size(cm)); 
    writhe(trial) = levittWrithe(cm + noise); 
    writheabs(trial) = levittWritheAbs(cm + noise);
    
    torsion(trial) = abs(lewinerTorsion(cm + noise, tauIdx+q, q));
end

writhe = abs(writhe); writhe_noiseless = abs(writhe_noiseless);

Wr(j,idx) = writhe_noiseless; WrA(j,idx) = writheabs_noiseless; Tor(j,idx) = torsion_noiseless;
noiseWr(j,idx) = mean(writhe); noiseWrStd(j,idx) = std(writhe); 
noiseWrA(j,idx) = mean(writheabs); noiseWrAStd(j,idx) = std(writheabs); 
noiseTor(j,idx) = mean(torsion); noiseTorStd(j,idx) = std(torsion); 

%{
figure; 
bar([writhe_noiseless, mean(writhe), writheabs_noiseless, mean(writheabs)]); 
hold on;
errorbar(...
    [writhe_noiseless, mean(writhe), writheabs_noiseless, mean(writheabs)], ...
    [0, std(writhe), 0, std(writheabs)], 'k', 'LineStyle', 'none');
grid on; title(['added noise variance = ' num2str(noisesz)]); 
ylabel('writhe'); 
xticklabels({'writhe', 'noisy writhe', 'writhe abs', 'noisy writhe abs'});
%}

end
idx/N
end
%%

ptstyle = {'o', '^', 's'};
lnstyle = {'--', '-.', ':'};
lncolor = {'b', 'r', 'k'};
lw = [1.5 1 .5];

figure; 

subplot(1,3,1); 
lgd = cell(1, 2*J);
for j = fliplr(1:J)
    errorbar(Wr, noiseWr(j,:), noiseWrStd(j,:), [ptstyle{j} lncolor{j}], 'LineWidth', lw(j)); 
    hold on;
    [fo, gof] = fit(Wr', noiseWr(j,:)', 'poly1'); 
    plot([min(Wr) max(Wr)], fo([min(Wr) max(Wr)]), [lnstyle{j} lncolor{j}]); 
    lgd{2*j} = ['noise var. = ' num2str(noisesz(j)^2)];
    lgd{2*j-1} = [num2str(fo.p1) '*x + ' num2str(fo.p2) ': R^2=' num2str(gof.rsquare)];
end
lgd = fliplr(lgd);
xlabel('Writhe Magnitude'); ylabel('Writhe Magnitude, added noise'); 
grid on; title('A) Writhe'); legend(lgd, 'Location', 'southoutside');

subplot(1,3,2); 
for j = fliplr(1:J)
    errorbar(WrA, noiseWrA(j,:), noiseWrAStd(j,:), [ptstyle{j} lncolor{j}], 'LineWidth', lw(j)); 
    hold on;
    [fo, gof] = fit(WrA', noiseWrA(j,:)', 'poly1'); 
    plot([min(WrA) max(WrA)], fo([min(WrA) max(WrA)]), [lnstyle{j} lncolor{j}]); 
    lgd{2*j} = ['noise var. = ' num2str(noisesz(j)^2)];
    lgd{2*j-1} = [num2str(fo.p1) '*x + ' num2str(fo.p2) ': R^2=' num2str(gof.rsquare)];
end
lgd = fliplr(lgd);
xlabel('Absolute Writhe'); ylabel('Absolute Writhe, added noise'); 
grid on; title('B) Absolute Writhe'); legend(lgd, 'Location', 'southoutside'); 

subplot(1,3,3); 
for j = fliplr(1:J)
    errorbar(Tor, noiseTor(j,:), noiseTorStd(j,:), [ptstyle{j} lncolor{j}], 'LineWidth', lw(j)); 
    hold on;
    [fo, gof] = fit(Tor', noiseTor(j,:)', 'poly1'); 
    plot([min(Tor) max(Tor)], fo([min(Tor) max(Tor)]), [lnstyle{j} lncolor{j}]); 
    lgd{2*j} = ['noise var. = ' num2str(noisesz(j)^2)];
    lgd{2*j-1} = [num2str(fo.p1) '*x + ' num2str(fo.p2) ': R^2=' num2str(gof.rsquare)];
end
lgd = fliplr(lgd);
xlabel('Max Torsion Magnitude'); ylabel('Max Torsion Magnitude, added noise'); 
grid on; title('C) Torsion'); legend(lgd, 'Location', 'southoutside'); 