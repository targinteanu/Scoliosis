N = 33;
Wr = zeros(1, N); WrA = zeros(1, N); Tor = zeros(1, N); 
noiseWr = zeros(1, N); noiseWrA = zeros(1, N); noiseTor = zeros(1, N); 
noiseWrStd = zeros(1, N); noiseWrAStd = zeros(1, N); noiseTorStd = zeros(1, N); 
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
noisesz = 5;
writhe = zeros(1, numtrials); 
writheabs = zeros(1, numtrials);
torsion = zeros(1, numtrials);

for trial = 1:numtrials
    noise = noisesz*randn(size(cm)); 
    writhe(trial) = levittWrithe(cm + noise); 
    writheabs(trial) = levittWritheAbs(cm + noise);
    
    torsion(trial) = abs(lewinerTorsion(cm + noise, tauIdx+q, q));
end

writhe = abs(writhe); writhe_noiseless = abs(writhe_noiseless);

Wr(idx) = writhe_noiseless; WrA(idx) = writheabs_noiseless; Tor(idx) = torsion_noiseless;
noiseWr(idx) = mean(writhe); noiseWrStd(idx) = std(writhe); 
noiseWrA(idx) = mean(writheabs); noiseWrAStd(idx) = std(writheabs); 
noiseTor(idx) = mean(torsion); noiseTorStd(idx) = std(torsion); 

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

idx/N
end
%%

figure; 

subplot(1,3,1); errorbar(Wr, noiseWr, noiseWrStd, 'o'); 
xlabel('Writhe'); ylabel('Writhe, added noise'); 
grid on; title(['added noise variance = ' num2str(noisesz)]); 

subplot(1,3,2); errorbar(WrA, noiseWrA, noiseWrAStd, 'o'); 
xlabel('Abs Writhe'); ylabel('Abs Writhe, added noise'); 
grid on; title(['added noise variance = ' num2str(noisesz)]); 

subplot(1,3,3); errorbar(Tor, noiseTor, noiseTorStd, 'o'); 
xlabel('Max Abs Torsion'); ylabel('Max Abs Torsion, added noise'); 
grid on; title(['added noise variance = ' num2str(noisesz)]); 