patient = [6, 28, 1, 3, 12, 15, 22]; 
Range = [3:7; 1:5; 1:5; 3:7; 3:7; 3:7; 3:7]; 

writhes = zeros(size(patient)); abswrithes = zeros(size(patient)); 
torsions = zeros(size(patient)); 

for p = 1:length(patient)
    load(['points_pat' num2str(patient(p)) '.mat']);
    writhes(p) = levittWrithe(cm, Range(p,:)); 
    abswrithes(p) = levittWritheAbs(cm, Range(p,:));
    torsions(p) = lewinerTorsion(cm, 3, 2); 
end

%%
patient_ordered = 1:length(patient);

figure; 
subplot(1,2,1); bar([abs(writhes); abswrithes]'); %grid on; 
xticklabels(patient_ordered); 
legend('Writhe (positive)', 'Absolute Writhe (negative)', 'Location', 'northwest');
xlabel('Patient'); ylabel('Writhe Magnitude'); title('A) Writhe and Absolute Writhe');
subplot(1,2,2); bar(abs(torsions)); grid on; 
xticklabels(patient_ordered);
xlabel('Patient'); ylabel('Torsion'); title('B) Geometric Torsion');

%% image resolution 
fp = uigetdir;
info = dicominfo([fp, '\IM-0001-0001-0001.dcm']);
xdim = info.PixelSpacing(1); ydim = info.PixelSpacing(2);
zdim = info.SliceThickness;
[xdim ydim zdim]

%% hypothesis test 
[h1, p1] = ttest2(abs(writhes(1:2)), abs(writhes(4:end)), 'Tail', 'left', 'Vartype', 'unequal')
[h2, p2] = ttest2(abswrithes(1:2), abswrithes(3:end), 'Tail', 'left', 'Vartype', 'unequal')