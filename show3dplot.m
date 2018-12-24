% shows 3D plots of selected patients demonstrating spine centerline and
% orientation of spinous processes. Also returns twist and writhe.
% Vertebrae of selected patients must be segmented previously. 

fp = uigetdir; 
load([fp, '\vertebrae.mat']);
exist([fp, '\points.mat'], 'file')
Vertebrae = Vertebrae(20:25); % T12 to L1

vertebrae = cell(size(Vertebrae));
for n = 1:length(vertebrae)
    vertebrae{n} = unpackage(Vertebrae{n}, true);
end

info = dicominfo([fp, '\IM-0001-0001-0001.dcm']);
xdim = info.PixelSpacing(1); ydim = info.PixelSpacing(2);
zdim = info.SliceThickness; 

%%
cm = zeros(length(vertebrae), 3);
sp = cm;
for n = 1:length(vertebrae)
    cm(n, :) = CenterOfMass3(vertebrae{n});
    sp_idx = sum(sum(vertebrae{n},2),3); 
    
    sp_idx = find(sp_idx); 
    sp_idx = sp_idx(end); % if NOT posterior/anterior inverted 
    %sp_idx = sp_idx(1); % if posterior/anterior inverted
    
    sp_pln = squeeze(vertebrae{n}(sp_idx,:,:));
    sp(n,:) = [sp_idx, CenterOfMass(sp_pln)];
end

cm(:,2) = cm(:,2)*xdim; cm(:,1) = cm(:,1)*ydim; cm(:,3) = cm(:,3)*zdim;
sp(:,2) = sp(:,2)*xdim; sp(:,1) = sp(:,1)*ydim; sp(:,3) = sp(:,3)*zdim;

%%

U = sp-cm; U = U(2:end,:);
dX = cm(2:end,:) - cm(1:end-1,:);
    % make U perp dX
    Uparl = diag((U*dX')./(dX*dX')) .* dX;
    U = U - Uparl;
    % make U unit
    normU = sqrt(diag(U*U'));
    %U = U./normU; 
U = [0 0 0; U];
V = U + cm;

figure; grid on; hold on;
for n = 2:length(cm(:,1))
    vector = [-cm(n,:); -V(n,:)];
    plot3(vector(:,2), vector(:,1), vector(:,3), '-or', 'LineWidth', 1);
end
plot3(-cm(:,2), -cm(:,1), -cm(:,3), 'b', 'LineWidth', 2);
xlabel('coronal (mm)'); ylabel('sagittal (mm)'); zlabel('axial (mm)');

Writhe = getWrithe(cm)
Twist = getTwist(cm, sp)