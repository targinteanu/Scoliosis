load('spines_XYZ.mat');

N = 33;
Torsions = zeros(N, 1); 

for idx = 1:N
    xyz = spinesXYZ{idx, 2:end};
    x = xyz(1:3:end); 
    y = xyz(2:3:end); 
    z = xyz(3:3:end); 
    p = [x;y;z]';
    
    % estimate 2nd derivative to get neutral/apical vertebrae 
    q = 2;
    vertebrae = (1+q):(17-q);
    concavities = zeros(size(vertebrae));
    for vertebra = vertebrae
        [~, concavity] = lewinerTorsion(p, vertebra, q);
        concavities(vertebra-q) = norm(concavity);
    end
    [~, neutral] = min(concavities); neutral = vertebrae(neutral);
    [~, apex] = max(concavities); apex = vertebrae(apex);
    q = abs(apex - neutral); q = min(q, 17-neutral);
    
    % get torsion at neutral vertebra with window size determined by apical
    Torsions(idx) = lewinerTorsion(p, neutral, q);
end