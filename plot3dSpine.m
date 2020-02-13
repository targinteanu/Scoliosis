function plot3dSpine(cm, dir, cmLineSpec, dirLineSpec, LW)

if sum(size(dir)==1) % sp is xy-plane angle
    dir = dir * pi / 180; % deg to rad 
    U = [cos(dir), sin(dir), zeros(size(dir))];
else % sp is 3D coordinates 
    %U = dir-cm; 
    U = dir;
end
if nargin < 4
    dirLineSpec = {'-r', 'sr'};
    if nargin < 3
        cmLineSpec = {'-b', 'sb'};
    end
end
if nargin < 5
    LW = 2;
end

U = .5* (U(2:end,:) + U(1:(end-1),:));
dX = cm(2:end,:) - cm(1:end-1,:);
    % make U perp dX
    Uparl = zeros(size(U));
    for n = 1:size(U,1)
        u = U(n,:); dx = dX(n,:);
        Uparl(n,:) = ((u*dx') / (dx*dx')) * dx;
    end
    %Uparl = diag((U*dX')./(dX*dX')) .* dX;
    U = U - Uparl;
    %arrayfun(@(n) acos(U(n,:)*dX(n,:)'), 1:size(U,1))
    % make U unit
    normU = sqrt(diag(U*U'));
    U = U./normU; 
    % make U big enough to display
    U = U * 1*abs(cm(1,3)-cm(end,3))/size(cm,1);
%U = [0 0 0; U];
cm2 = .5* (cm(2:end,:) + cm(1:(end-1),:));
V = (U + cm2)*1;
%V = [0 0 0; V];

%figure; 
grid on; hold on;
for n = 1:length(cm2(:,1))
    %vector = [-cm(n,:); -V(n,:)];
    vector = [cm2(n,:); V(n,:)];
    plot3(vector(:,1), vector(:,2), vector(:,3), dirLineSpec{1}, 'LineWidth', LW*.5);
    plot3(V(:,1), V(:,2), V(:,3), dirLineSpec{2}, 'LineWidth', LW*.5);
    plot3(cm2(:,1), cm2(:,2), cm2(:,3), cmLineSpec{2}, 'LineWidth', LW*.5);
end
%plot3(-cm(:,1), -cm(:,2), -cm(:,3), 'b', LineWidth', 2);
plot3(cm(:,1), cm(:,2), cm(:,3), cmLineSpec{1}, 'LineWidth', LW);
%xlabel('coronal'); ylabel('sagittal'); zlabel('axial');
xlabel('x'); ylabel('y'); zlabel('z');

end