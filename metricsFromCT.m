patients = {'1', '2', '3', '4', '6', '12', '28'};
patients = {'1', '2', '3', '4'}

for i = 1:length(patients)
    patient = patients{i}
    load(['points_pat' patient '.mat'])
    Writhe = levittWrithe(cm)
end

%%

load(['points_pat' patient '.mat'])

Writhe = levittWrithe(cm)

