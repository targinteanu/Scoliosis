load('spines_XYZ_writhe.mat');
Writhes = spinesXYZwrithe{:,2};

N = 33;
Torsions = zeros(N, 1); neutrals = Torsions; qs = neutrals;

for idx = 1:N
    %%
    for idx = [1 3 9 10 14 17]
    
    xyz = spinesXYZ{idx, 2:end};
    x = xyz(1:3:end); 
    y = xyz(2:3:end); 
    z = xyz(3:3:end); 
    p = [x;y;z]';
    
    % estimate 2nd derivative to get neutral/apical vertebrae 
    figure; hold on; grid on;
    for q = 2:5
    vertebrae = (1+q):(17-q);
    concavities = zeros(size(vertebrae));
    for vertebra = vertebrae
        [~, concavity] = lewinerTorsion(p, vertebra, q);
        concavities(vertebra-q) = norm(concavity);
    end
    [~, neutral] = min(concavities); neutral = vertebrae(neutral);
    [~, apex] = max(concavities); apex = vertebrae(apex);
    q = abs(apex - neutral); q = min(q, 17-neutral);
    plot(vertebrae, concavities);
    end
    legend('2', '3', '4', '5'); ylim([.85, 1.1]);
    title([num2str(idx), ' | ', num2str(Writhes(idx))]);
    end
    %%
    
    % get torsion at neutral vertebra with window size determined by apical
    Torsions(idx) = lewinerTorsion(p, neutral, q);
    
    qs(idx) = q; neutrals(idx) = neutral;
end
%%

figure; plot(Torsions, '.'); hold on; plot(Writhes, '.'); grid on;