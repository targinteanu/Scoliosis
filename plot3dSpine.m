function plot3dSpine(cm, dir)

if sum(size(dir)==1) % sp is xy-plane angle
    dir = dir * pi / 180; % deg to rad 
    U = [cos(dir), sin(dir), zeros(size(dir))];
else % sp is 3D coordinates 
    U = dir-cm; 
end

U = U(2:end,:);
dX = cm(2:end,:) - cm(1:end-1,:);
    % make U perp dX
    Uparl = diag((U*dX')./(dX*dX')) .* dX;
    U = U - Uparl;
    % make U unit
    normU = sqrt(diag(U*U'));
    U = U./normU; 
    % make U big enough to display
    U = U * 4*abs(cm(1,3)-cm(end,3))/size(cm,1);
U = [0 0 0; U];
V = U + cm;

%figure; 
grid on; hold on;
for n = 2:length(cm(:,1))
    %vector = [-cm(n,:); -V(n,:)];
    vector = [cm(n,:); V(n,:)];
    plot3(vector(:,2), vector(:,1), vector(:,3), '-or', 'LineWidth', 1);
end
%plot3(-cm(:,2), -cm(:,1), -cm(:,3), 'b', 'LineWidth', 2);
plot3(cm(:,2), cm(:,1), cm(:,3), 'b', 'LineWidth', 2);
%xlabel('coronal'); ylabel('sagittal'); zlabel('axial');
xlabel('x'); ylabel('y'); zlabel('z');

end