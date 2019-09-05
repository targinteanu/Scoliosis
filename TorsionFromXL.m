% Loads N spines from spines_XYZ.mat which contains table spinesXYZ. Gets
% torsion for each spine using method described by Kadoury: least squares
% cubics are fit to groups of vertebrae, and the "neutral" vertebra is
% estimated as the one with the smallest second derivative while the
% "apical" is estimated as the one with the largest. Then, torsion is
% estimated at the "neutral" vertebra by fitting a cubic from the "apical."  
% Also, estimate the maximum torsion of the spine. 

%[num, txt] = xlsread('Writhe-pre-post.xlsx');
%XYZ = num(:, 6:end); 
[num, txt] = xlsread('Writhe-pre-post_new-metrics.csv');
N = 32;
num = num(1:N, :);
XYZ = num(:, 13:end); 

shapecluster = num(:,1);
writhe = num(:,6); abswrithe = num(:,7); 
tor1 = num(:,8); tor2 = num(:,9); torglob = num(:,10); 
twist = num(:,11); writhetwist = num(:,12);

Torsions = zeros(N, 1); % neutral-to-apical torsion 
maxTorsions = Torsions; % max torsion on the spine 
%maxTorsions2 = maxTorsions;
torsionlocs = Torsions; torsionlocs2 = torsionlocs;
neutrals = Torsions; apicals = neutrals; qs = neutrals;

% trying different ways to get apicals 
apicals2 = zeros(length(apicals), 2);
newTorsions = Torsions;

debugmode = false;

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

    % higher apex = lowest 1st derivative excluding T12, L1
    [~, apex] = min(d1(1:(end-2))); apex = vertebrae(apex); % exclude T12, L1
    % neutral = lowest 2nd derivative below higher apex 
    nverts = vertebrae(vertebrae >= apex);
    [~, neutral] = min( d2(vertebrae >= apex) ); neutral = nverts(neutral);

    q = abs(apex - neutral); % size q based on apex-to-neutral 
    q = min([q, 17-neutral, neutral-1]); % do not go out of bounds above 17 or below 1
    
    % get torsion at neutral vertebra with window size determined by apical
    Torsions(idx) = lewinerTorsion(p, neutral, q);
    
    qs(idx) = q; neutrals(idx) = neutral; apicals(idx) = apex; 
    % store these to look at which vertebrae were selected as "neutral"/"apical"
    
    %% new ways of getting apicals
    dist_from_z = sqrt(x.^2 + y.^2);
    %%{
    % use findpeaks
    [~, apexes, w, pp] = findpeaks(dist_from_z); 
    % if there are more than 2, get rid of the least prominent
    while length(apexes) > 2
        [~,minIdx] = min(pp); 
        apexes = apexes([1:(minIdx-1), (minIdx+1):end]);
        w = w([1:(minIdx-1), (minIdx+1):end]);
        pp = pp([1:(minIdx-1), (minIdx+1):end]);
    end
    %}
    %{
    % use max
    [~, apexes(1)] = max(dist_from_z(1:neutral));
    [~, apexes(2)] = max(dist_from_z(neutral:end)); apexes(2) = apexes(2) + neutral-1;
    %}
    apicals2(idx,:) = apexes;
    
    if (length(apexes) > 1)
        qq = abs(apexes - neutral);
        q3 = 2*lcm(qq(1),qq(2)) + 1;
        q2 = 2*qq(1)*qq(2) + 1;
        
        pwin = p(apexes(1):apexes(2),:);
        zwin = linspace(pwin(1,3), pwin(end,3), q2);
        xwin = spline(pwin(:,3), pwin(:,1), zwin);
        ywin = spline(pwin(:,3), pwin(:,2), zwin);
        %{
        pinterp = zeros(size(p,1)*q3, 3);
        for col = 1:3
            pinterp(:,col) = interp(p(:,col), q3);
        end
        %}
        
%        pcran = p(apexes(1):neutral); pcaud = p(neutral:apexes(2));
        qcran = qq(2); qcaud = qq(1);
%        pcran_interp = zeros(size(pcran,1)*qcran, 3);
%        pcaud_interp = zeros(size(pcaud,1)*qcaud, 3);
        pcran_interp = zeros(size(p,1)*qcran, 3);
        pcaud_interp = zeros(size(p,1)*qcaud, 3);
        for col = 1:3
            pcran_interp(:,col) = interp(p(:,col), qcran);
            pcaud_interp(:,col) = interp(p(:,col), qcaud);
        end
        icran = ([apexes(1), neutral]-1)*qcran + 1;
        icaud = ([neutral, apexes(2)]-1)*qcaud + 1;
        pcran = pcran_interp(icran(1):icran(2), :);
        pcaud = pcaud_interp(icaud(1):icaud(2), :);
        pinterp = [pcran; pcaud(2:end,:)];
        
        newTorsions(idx) = lewinerTorsion(pinterp);
    end
    if debugmode
        figure; plot3(x, y, z, 'o'); hold on; grid on;
        plot3(x([neutral, apexes]), y([neutral, apexes]), z([neutral, apexes]), 'x');
        %plot3(xwin, ywin, zwin, '-+');
        plot3(pinterp(:,1), pinterp(:,2), pinterp(:,3), '.');
        %plot3(pcran(:,1), pcran(:,2), pcran(:,3), '.');
        %plot3(pcaud(:,1), pcaud(:,2), pcaud(:,3), '.');
    end
    
    if debugmode
        figure; plot(dist_from_z); grid on; hold on; 
        xlabel('vertebra'); ylabel('distance from z axis'); 
        plot(apexes, dist_from_z(apexes), 'o'); 
        plot(neutral, dist_from_z(neutral), 'o'); 
        plot(apex, dist_from_z(apex), 'x');
        title(num2str(idx));
    end
end
%%

figure; plot(Torsions, '.'); grid on; hold on; plot(maxTorsions, '.'); plot(newTorsions, '.');
xlabel('patient'); ylabel('Torsion');
legend('neutral-apical', 'maximum', 'new');

%%
debugmode = true;

checkTorsions = max(abs(Torsions - torglob))
checkMaxTorsions = max(abs(maxTorsions - tor1))

newTorsions(apicals2(:,1)==apicals2(:,2)) = nan;

apicalsLow = neutrals - (apicals - neutrals);
figure; plot(neutrals, '-k'); hold on; grid on;
plot([apicals, apicalsLow], ':r'); plot(apicals2, '--b');
xlabel('patient'); ylabel('vertebra'); 
legend('neutral from d2', 'apical from d1', ...
    'apical from d1', 'apical from z-dist', 'apical from z-dist');