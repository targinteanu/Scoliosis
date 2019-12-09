N = 4; cm = cell(1,N); theta = cell(1,N);
npts = 17;
amax = 180;
lbl = {'a) ', 'b) ', 'c) ', 'd) '};

s = linspace(0, 1, npts)';

%cm{1} = [zeros(size(s)), zeros(size(s)), s];
cm{1} = [cos(2*pi*s), zeros(size(s)), s];
theta{1} = zeros(size(s));
cm{2} = cm{1};
theta{2} = linspace(amax, 0, npts)';
cm{3} = [cos(2*pi*s), sin(2*pi*s), s];
theta{3} = theta{1};
cm{4} = cm{3};
theta{4} = theta{2};

figure('Position', [50 50 1400 500]); 
for i = 1:N
    subplot(1,N,i);
    plot3dSpine(cm{i}, theta{i}); 
    Tw = getTwist(cm{i}, theta{i});
    Wr = levittWrithe(cm{i});
    title([lbl{i}, 'Twist = ' num2str(Tw) ', Writhe = ' num2str(Wr)]);
    view([20, 10]);
    xlim([-2, 2]); ylim([-2, 2]); zlim([0, 1]);
end

%%
%{
a = 1:2000; Tw = zeros(size(a)); Twact = Tw;
for ang = a
  
cm{2} = [zeros(size(s)), zeros(size(s)), s];
theta{2} = linspace(0, ang * pi / 180, npts)';

i = 2; %figure; 
%    plot3dSpine(cm{i}, theta{i}); 
    Tw(ang) = getTwist(cm{i}, cm{i}+[cos(theta{i}), sin(theta{i}), zeros(size(theta{i}))]);
    
    Twact(ang) = -ang/(360);
end
figure; plot(a/360, [Tw; Twact]); %xticks(0:45:720);
D = cm{i}(2:end,:) - cm{i}(1:(end-1),:);
DD = D*D'; DD = diag(DD);
mean(DD)
%}