% Loads N spines from spines_XYZ.mat which contains table spinesXYZ. Gets
% writhe and maximum torsion for each spine. Plots the results for
% clustering. 

load('spines_XYZ.mat');

N = 33;
writhes = zeros(N,1); torsions = zeros(N,1); 
for idx = 1:N
    % get center points 
    xyz = spinesXYZ{idx, 2:end};
    x = xyz(1:3:end); 
    y = xyz(2:3:end); 
    z = xyz(3:3:end); 
    cm = [x;y;z]';
    
    % get writhe 
    writhes(idx) = levittWrithe(cm);
    
    % get torsion 
    q = 2; % 2 vertebrae above, 2 vertebrae below, current vertebra -> 5 points to fit each cubic
    vertebrae = (1+q):(length(x)-q);
    tau = arrayfun(@(v) lewinerTorsion(cm, v, q), vertebrae);
    [~,tauvert] = max(abs(tau)); torsions(idx) = tau(tauvert); % maximum torsion
end

figure; plot(writhes, torsions, '.'); grid on;
xlabel('Writhe'); ylabel('Max Torsion');