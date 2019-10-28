function [Torsion_from_deriv, neutral, apical_from_deriv, ...
    Torsion_from_dist, apical_from_dist, maxTorsion, maxTorsionLoc] ...
    = kadouryTorsion(p, q)
% uses Kadoury method on a sampled 3D space curve to select parameters with
% which to calculate lewinerTorsion. 

if nargin < 2
    q = 4; % 2 vertebrae above, 2 vertebrae below, current vertebra -> 5 points to fit each cubic
end

deriv_front = @(cf, n) factorial(n)*cf(:,(end-n));  

x = p(:,1); y = p(:,2); z = p(:,3);


    % lewiner on raw data ------------------------------------------------
    vertebrae = (1+q):(size(p,1)-q); 
    tau = zeros(size(vertebrae)); % local torsion at each point on the spine 
    d2 = zeros(size(vertebrae)); d1 = d2;
    for vertebra = vertebrae
        [tau(vertebra-q), d, dd] = lewinerTorsion(p, vertebra, q);
        d2(vertebra-q) = norm(dd); d1(vertebra-q) = norm(d);
    end
    [~, maxIdx] = max(abs(tau)); 
    maxTorsion = tau(maxIdx); maxTorsionLoc = vertebrae(maxIdx);
    % --------------------------------------------------------------------
    
    
    % perform spline interpolation ---------------------------------------
    pdisp = diff(p);
    L = sqrt(diag(pdisp*pdisp')); 
    L = [0; cumsum(L)]; % smaller L = closer to sacrum
    interpfactor = 1;
    Ltot = (L(end)-L(1)); interpext = 2/length(L);
    Linterp = linspace(L(1)-(interpext*Ltot), L(end)+(interpext*Ltot), ...
        length(L)*(2*interpext + 1)*interpfactor)';
    
    ppall = [spline(L, x), spline(L, y), spline(L, z)];
    cfs = {ppall.coefs};
    
    dr = cell2mat(arrayfun(@(var) deriv_front(cfs{var}, 1), 1:3, 'UniformOutput', false));
    ddr = cell2mat(arrayfun(@(var) deriv_front(cfs{var}, 2), 1:3, 'UniformOutput', false));
    dddr = cell2mat(arrayfun(@(var) deriv_front(cfs{var}, 3), 1:3, 'UniformOutput', false));
    
    tauspline = zeros(size(dr,1), 1);
    for i = 1:length(tauspline)
        cx = cross(dr(i,:), ddr(i,:));
        tauspline(i) = -(cx * dddr(i,:)')/(norm(cx)^2);
    end
    
%    df1 = arrayfun(@(i) norm(dr(i,:)), 1:size(dr,1));
%    df2 = arrayfun(@(i) norm(ddr(i,:)), 1:size(ddr,1));
%    df3 = arrayfun(@(i) norm(dddr(i,:)), 1:size(dddr,1));
    
    pinterp = cell2mat(arrayfun(@(j) ppval(ppall(j), Linterp), 1:3, 'UniformOutput', false));
    % --------------------------------------------------------------------
    
    
    % lewiner on spline interpolation ------------------------------------
    q = q*interpfactor; 
    vinterp = (1+q):(size(pinterp,1)-q); 
    tauinterp = zeros(size(vinterp)); % local torsion at each point on the spine 
    d2int = zeros(size(vinterp)); d1int = d2int; d3int = d2int;
    for vertebra = vinterp
        [tauinterp(vertebra-q), d, dd, ddd] = lewinerTorsion(pinterp, vertebra, q);
        d1int(vertebra-q) = norm(d); d2int(vertebra-q) = norm(dd); d3int(vertebra-q) = norm(ddd);
    end
        
    % higher apex = lowest 1st derivative 
    [~, apex] = min(d1(1:(end-2))); apex = vertebrae(apex); % exclude T12, L1
    
    % neutral = lowest 2nd derivative below higher apex 
    nverts = vertebrae(vertebrae >= apex);
    [~, neutral] = min( d2(vertebrae >= apex) ); neutral = nverts(neutral);
    % --------------------------------------------------------------------
    
    
    %% use derivatives to get apex and neutral ----------------------------
    
    % get lower apex using spline  
    aidx = Linterp(vinterp) > L(neutral); averts = vinterp(aidx);
    [~, ap] = min(d1int(aidx));
    ap = averts(ap); 
    
    % convert interpolated apex to real 
    [~, ap] = min(abs(L - Linterp(ap)));
    if ap == neutral
        % look for peaks instead 
        [~, ap, w, pp] = findpeaks(-d1int(aidx));
        % if there are more than 2, get rid of the farthest
        while length(ap) > 1
            %[~,minIdx] = min(pp);
            [~,maxIdx] = max(abs(ap - neutral));
            ap = ap([1:(maxIdx-1), (maxIdx+1):end]);
            w = w([1:(maxIdx-1), (maxIdx+1):end]);
            pp = pp([1:(maxIdx-1), (maxIdx+1):end]);
        end
        % if there aren't any, pick the -greatest value
        if isempty(ap)
            ap = NaN;
        end
    end
    
    % store variables 
    %neutral_from_deriv = neutral; 
    apical_from_deriv = [apex, ap];
    if ~isnan(ap)
        Torsion_from_deriv = customTorsion(p, neutral, [apex, ap]);
    end
    
    % --------------------------------------------------------------------
    
    
    %% get apicals from distance from z axis ------------------------------
    dist_from_z = sqrt(x.^2 + y.^2);
    sides = {1:neutral, neutral:length(dist_from_z)};
    apexes = [0 0];
    for s = 1:length(sides)
        [~, ap, w, pp] = findpeaks(dist_from_z(sides{s}));
        % if there are more than 2, get rid of the farthest
        while length(ap) > 1
            %[~,minIdx] = min(pp);
            [~,maxIdx] = max(abs(ap - neutral));
            ap = ap([1:(maxIdx-1), (maxIdx+1):end]);
            w = w([1:(maxIdx-1), (maxIdx+1):end]);
            pp = pp([1:(maxIdx-1), (maxIdx+1):end]);
        end
        % if there aren't any, pick the greatest value
        if isempty(ap)
            [~, ap] = max(dist_from_z(sides{s}));
        end 
        apexes(s) = ap;
    end
    apexes(2) = apexes(2) + neutral - 1;
    apexes(apexes == neutral) = nan;

    apical_from_dist = apexes;
    if ~isnan(sum(apexes))
        Torsion_from_dist = customTorsion(p, neutral, apexes);
    end

%%
function T = customTorsion(p, neutral, apexes)

    qq = abs(apexes - neutral);
    qcran = qq(2); qcaud = qq(1);
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
    pintp = [pcran; pcaud(2:end,:)];
    
    T = lewinerTorsion(pintp);

end

end