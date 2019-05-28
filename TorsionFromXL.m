% Loads N spines from spines_XYZ.mat which contains table spinesXYZ. Gets
% torsion for each spine using method described by Kadoury: least squares
% cubics are fit to groups of vertebrae, and the "neutral" vertebra is
% estimated as the one with the smallest second derivative while the
% "apical" is estimated as the one with the largest. Then, torsion is
% estimated at the "neutral" vertebra by fitting a cubic from the "apical."  
% Also, estimate the maximum torsion of the spine. 

[num, txt] = xlsread('Writhe-pre-post.xlsx');
XYZ = num(:, 6:end); 
N = 32;

Torsions = zeros(N, 1); % neutral-to-apical torsion 
maxTorsions = Torsions; % max torsion on the spine 
%maxTorsions2 = maxTorsions;
torsionlocs = Torsions; torsionlocs2 = torsionlocs;
neutrals = Torsions; apicals = neutrals; qs = neutrals;

for idx = 1:N
    %%    
    % get center points 
    % get center points 
    x = XYZ(idx, 1:3:51); 
    y = XYZ(idx, 2:3:51); 
    z = XYZ(idx, 3:3:51); 
    p = [x;y;z]';
    
    % estimate 2nd derivative to get neutral/apical vertebrae 
    q = 4; % 2 vertebrae above, 2 vertebrae below, current vertebra -> 5 points to fit each cubic
    vertebrae = (1+q):(17-q); 
    tau = zeros(size(vertebrae)); % local torsion at each point on the spine 
    d2 = zeros(size(vertebrae)); d1 = d2;
    for vertebra = vertebrae
        [tau(vertebra-q), d, dd] = lewinerTorsion(p, vertebra, q);
        d2(vertebra-q) = norm(dd); d1(vertebra-q) = norm(d);
    end
    
    [~, maxIdx] = max(abs(tau)); 
    maxTorsions(idx) = tau(maxIdx); torsionlocs(idx) = vertebrae(maxIdx);
%    [~, neutral] = min(d2); neutral = vertebrae(neutral);
%    [~, apex] = max(d2); apex = vertebrae(apex);
    [~, apex] = min(d1(1:(end-2))); apex = vertebrae(apex); % exclude T12, L1
    nverts = vertebrae(vertebrae >= apex);
    [~, neutral] = min( d2(vertebrae >= apex) ); neutral = nverts(neutral);

    q = abs(apex - neutral); % size q based on apex-to-neutral 
    q = min([q, 17-neutral, neutral-1]); % do not go out of bounds above 17 or below 1
    %%
    % get torsion at neutral vertebra with window size determined by apical
    Torsions(idx) = lewinerTorsion(p, neutral, q);
    
    qs(idx) = q; neutrals(idx) = neutral; apicals(idx) = apex; 
    % store these to look at which vertebrae were selected as "neutral"/"apical"
end
%%

figure; plot(Torsions, '.'); grid on; hold on; plot(maxTorsions, '.');
xlabel('patient'); ylabel('Torsion');
legend('neutral-apical', 'maximum');