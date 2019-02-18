load('spines_XYZ.mat');

N = 33;
levittWrithes = zeros(N, 1); diffWrithes = levittWrithes; oldlevittWrithes = levittWrithes;
for idx = 1:N
    %%
    xyz = spinesXYZ{idx, 2:end};
    x = xyz(1:3:end); 
    y = xyz(2:3:end); 
    z = xyz(3:3:end); 
    %levittWrithe([x;y;z;]');
    %%
    levittWrithes(idx) = levittWrithe([x;y;z]');
    diffWrithes(idx) = getWrithe([x;y;z]');
    oldlevittWrithes(idx) = levittWritheOld([x;y;z;]');
    %figure; plot3(x, y, z); grid on;
    %title([num2str(Writhes(idx)), ' | ', num2str(altWrithes(idx))])
    %idx
end
%%
figure; plot(levittWrithes, 'o'); hold on; plot(diffWrithes, '.'); grid on;
plot(oldlevittWrithes, '.');

figure; plot(diffWrithes, levittWrithes, '.'); grid on;
xlabel('Writhe (discrete)'); ylabel('Writhe (Levitt)');

figure; plot(oldlevittWrithes, levittWrithes, '.'); grid on;
xlabel('Levitt Writhe (toren fix)'); ylabel('Levitt Writhe (DeTurck fix)');

%Writhes = real(Writhes);

spinesXYZwrithe = [spinesXYZ(:,1), table(levittWrithes)]

save('spines_XYZ_writhe.mat', 'spinesXYZ', 'spinesXYZwrithe');