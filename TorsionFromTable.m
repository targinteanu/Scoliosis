% Loads N spines from spines_XYZ.mat which contains table spinesXYZ. Gets
% torsion for each spine using method described by Kadoury: least squares
% cubics are fit to groups of vertebrae, and the "neutral" vertebra is
% estimated as the one with the smallest second derivative while the
% "apical" is estimated as the one with the largest. Then, torsion is
% estimated at the "neutral" vertebra by fitting a cubic from the "apical."  
% Also, estimate the maximum torsion of the spine. 

load('spines_XYZ.mat');

N = 33;
Torsions = zeros(N, 1); % neutral-to-apical torsion 
maxTorsions = Torsions; % max torsion on the spine 
neutrals = Torsions; qs = neutrals;

for idx = 1:N
    %%    
    % get center points 
    xyz = spinesXYZ{idx, 2:end};
    x = xyz(1:3:end); 
    y = xyz(2:3:end); 
    z = xyz(3:3:end); 
    p = [x;y;z]';
    
    % estimate 2nd derivative to get neutral/apical vertebrae 
    q = 2; % 2 vertebrae above, 2 vertebrae below, current vertebra -> 5 points to fit each cubic
    vertebrae = (1+q):(17-q); 
    tau = zeros(size(vertebrae)); % local torsion at each point on the spine 
    concavities = zeros(size(vertebrae));
    for vertebra = vertebrae
        [tau(vertebra-q), concavity] = lewinerTorsion(p, vertebra, q);
        concavities(vertebra-q) = norm(concavity);
    end
    [~, maxIdx] = max(abs(tau)); maxTorsions(idx) = tau(maxIdx);
    [~, neutral] = min(concavities); neutral = vertebrae(neutral);
    [~, apex] = max(concavities); apex = vertebrae(apex);
    q = abs(apex - neutral); % size q based on apex-to-neutral 
    q = min([q, 17-neutral, neutral-1]); % do not go out of bounds above 17 or below 1
    %%
    % get torsion at neutral vertebra with window size determined by apical
    Torsions(idx) = lewinerTorsion(p, neutral, q);
    
    qs(idx) = q; neutrals(idx) = neutral; 
    % store these to look at which vertebrae were selected as "neutral"/"apical"
end
%%

figure; plot(Torsions, '.'); grid on; hold on; plot(maxTorsions, '.');
xlabel('patient'); ylabel('Torsion');
legend('neutral-apical', 'maximum');