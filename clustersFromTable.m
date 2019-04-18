load('spines_XYZ.mat');

N = 33;
writhes = zeros(N,1); torsions = zeros(N,1); 
TAU = zeros(N, 13); MM = zeros(N, 289);
for idx = 1:N
    xyz = spinesXYZ{idx, 2:end};
    x = xyz(1:3:end); 
    y = xyz(2:3:end); 
    z = xyz(3:3:end); 
    cm = [x;y;z]';
    writhes(idx) = levittWrithe(cm);
    writheM;
    MM(idx,:) = M(:)';
    
    q = 2;
    vertebrae = (1+q):(length(x)-q);
    tau = arrayfun(@(v) lewinerTorsion(cm, v, q), vertebrae);
    [~,tauvert] = max(abs(tau)); torsions(idx) = tau(tauvert);
    TAU(idx,:) = tau;
    %figure; plot(vertebrae, tau)
    %{
    concavities = zeros(size(vertebrae));
    tau = zeros(size(vertebrae));
    for vertebra = vertebrae
        [tau(vertebra-q),~,concavity] = ...
            lewinerTorsion(P, vertebra, q);
        concavities(vertebra-q) = norm(concavity);
    end
    neutral = (concavities(2:(end-1)) < concavities(1:(end-2))) & ...
        (concavities(2:(end-1)) < concavities(3:end)); 
    neutral = find([0 neutral]); neutral = vertebrae(neutral);
    [~, apex] = max(concavities); apex = vertebrae(apex);
    q = abs(apex - neutral); q = min(q, length(x)-neutral);
    torsions(idx,:) = arrayfun(@(i) lewinerTorsion(P, neutral(i), q(i)), 1:length(q));
    %}
end