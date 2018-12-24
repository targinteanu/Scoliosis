% obtains center-of-mass and spinous process data from specified patients'
% segmented vertebrae, plots the results, and gets twist and writhe. 

fp = uigetdir; 
load([fp, '\vertebrae.mat']);

Vertebrae = Vertebrae(2:end); % if Vertebrae{1} is "none"

tf = false(size(Vertebrae));
for n = 1:length(Vertebrae)
    sz = size(Vertebrae{n}.Volume);
    tf(n) = length(sz) > 2;
end
vertebrae = Vertebrae(tf);

for n = 1:length(vertebrae)
    vertebrae{n} = unpackage(vertebrae{n}, true);
end

info = dicominfo([fp, '\IM-0001-0001-0001.dcm']);
xdim = info.PixelSpacing(1); ydim = info.PixelSpacing(2);
zdim = info.SliceThickness; 

%%
%{
for n = 1:length(vertebrae)
    vertebrae{n} = imfill(vertebrae{n}, 'holes');
    vertebrae{n} = medfilt3(vertebrae{n});
end
%}

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

wholespine = zeros(size(vertebrae{1}));
pts = wholespine;
for n = 1:length(vertebrae)
    wholespine = wholespine + vertebrae{n};
    pts(floor(cm(n,1)), floor(cm(n,2)), floor(cm(n,3))) = 1;
    pts(floor(sp(n,1)), floor(sp(n,2)), floor(sp(n,3))) = 1;
end

dilationelement = strel('sphere', 3);
pts = imdilate(pts, dilationelement);
spinepoints = wholespine+(pts*3);

cm(:,2) = cm(:,2)*xdim; cm(:,1) = cm(:,1)*ydim; cm(:,3) = cm(:,3)*zdim;
sp(:,2) = sp(:,2)*xdim; sp(:,1) = sp(:,1)*ydim; sp(:,3) = sp(:,3)*zdim;

%%
%{
figure; plot3(-cm(:,2), -cm(:,1), -cm(:,3)); 
xlabel('x'); ylabel('y'); zlabel('z');
grid on;
hold on; plot3(-sp(:,2), -sp(:,1), -sp(:,3));
%}

%%
figure; grid on; hold on;
for n = 1:length(cm(:,1))
    vector = [-cm(n,:); -sp(n,:)];
    plot3(vector(:,2), vector(:,1), vector(:,3), '-or');
end
plot3(-cm(:,2), -cm(:,1), -cm(:,3), 'b');
xlabel('x'); ylabel('y'); zlabel('z');

%% writhe 
Writhe = getWrithe(cm)
%% twist 
Twist = getTwist(cm, sp)