N = 4; cm = cell(1,N); theta = cell(1,N);
npts = 10;

s = linspace(0, 1, npts)';

cm{1} = [zeros(size(s)), zeros(size(s)), s];
theta{1} = zeros(size(s));
cm{2} = [zeros(size(s)), zeros(size(s)), s];
theta{2} = linspace(0, 90, npts)';
cm{3} = [cos(2*pi*s), sin(2*pi*s), s];
theta{3} = zeros(size(s));
cm{4} = [cos(2*pi*s), sin(2*pi*s), s];
theta{4} = linspace(0, 90, npts)';

figure; 
for i = 1:N
    subplot(1,N,i);
    plot3dSpine(cm{i}, theta{i}); 
    Tw = getTwist(cm{i}, cm{i}+[cos(theta{i}), sin(theta{i}), zeros(size(theta{i}))]);
    title(['Twist = ' num2str(Tw)]);
    view(3);
    xlim([-2, 2]); ylim([-2, 2]); zlim([0, 1]);
end

%%
a = 1:2000; Tw = zeros(size(a)); Twact = Tw;
for ang = a
cm{2} = [zeros(size(s)), zeros(size(s)), s];
theta{2} = linspace(0, ang, npts)';
i = 2; %figure; 
%    plot3dSpine(cm{i}, theta{i}); 
    Tw(ang) = getTwist(cm{i}, cm{i}+[cos(theta{i}), sin(theta{i}), zeros(size(theta{i}))]);
    
end
figure; plot(a/360, Tw); %xticks(0:45:720);
D = cm{i}(2:end,:) - cm{i}(1:(end-1),:);
DD = D*D'; DD = diag(DD);
mean(DD)