load('spines_XYZ.mat');

N = 33;
Writhes = zeros(N, 1); altWrithes = Writhes;
for idx = 1:N
    xyz = spinesXYZ{idx, 2:end};
    x = xyz(1:3:end); 
    y = xyz(2:3:end); 
    z = xyz(3:3:end); 
    Writhes(idx) = levittWrithe([x;y;z]');
    altWrithes(idx) = getWrithe([x;y;z]');
    %figure; plot3(x, y, z); grid on;
    %title([num2str(Writhes(idx)), ' | ', num2str(altWrithes(idx))])
    %idx
end

figure; plot(Writhes, '.'); hold on; plot(altWrithes, '.'); grid on;

%Writhes = real(Writhes);

spinesXYZwrithe = [spinesXYZ(:,1), table(Writhes)]

%save('spines_XYZ_writhe.mat', 'spinesXYZ', 'spinesXYZwrithe');