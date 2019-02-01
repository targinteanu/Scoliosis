function Torsion = getTorsion(r, range)

if nargin < 2
    range = 1:length(r(:,1));
end

dr = r(2:end,:) - r(1:(end-1),:);
ds = diag(dr*dr');
dr = [zeros(1,3); dr];
ddr = dr(2:end,:) - dr(1:(end-1),:);
ddr = [zeros(1,3); ddr];
dddr = ddr(2:end,:) - ddr(1:(end-1),:);
dddr = [zeros(1,3); dddr];

Torsion = zeros(size(range));
for t = range
    cx = cross(dr(t,:), ddr(t,:));
    Torsion(t) = -(cx*dddr(t,:)')/(norm(cx)^2);
end

end