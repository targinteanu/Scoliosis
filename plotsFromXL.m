% gets data from XL spreadsheet and draws specified plots 

%% load data

vertnames = [arrayfun(@(i) ['T' num2str(i)], 1:12, 'UniformOutput', 0), ...
    arrayfun(@(i) ['L' num2str(i)], 1:5, 'UniformOutput', 0)];

[num, txt] = xlsread('Writhe-pre-post_new-metrics.xlsx');

shapecluster = num(:,3); 
writhe_preop = num(:,4); 
writhe_postop = num(:,5);

abswrithe = num(:,6);
maxtorsion1 = num(:,7);
maxtor1vert = num(:,8);
maxtorsion2 = num(:,9); 
maxtor2vert = num(:,10); 
globtorsion = num(:,11);
vertNeutral = num(:,12); 
vertApical = num(:,13); 
twist = num(:,14);

N = length(shapecluster); 
lbl = arrayfun(@(n) num2str(n), 1:N, 'UniformOutput', false);

tbl = table(shapecluster, writhe_postop, writhe_preop, abswrithe, twist, ...
    globtorsion, vertNeutral, vertApical, ...
    maxtorsion1, maxtor1vert, maxtorsion2, maxtor2vert, ...
    'VariableNames', {...
    'ShapeCluster', 'WrithePostOp', 'Writhe', 'AbsoluteWrithe', 'Twist', ...
    'KadouryGlobalTorsion', 'NeutralVertebra', 'ApicalVertebra', ...
    'MaxTorsion', 'MaxTorsionVertebra', ...
    'Max2Torsion', 'Max2TorsionVertebra', ...
    });

%% feature plots 

feature = [3, 5; 3, 6];

for f = 1:length(feature(:,1))
    f1name = tbl(:, feature(f,1)).Properties.VariableNames{:};
    f2name = tbl(:, feature(f,2)).Properties.VariableNames{:};
    f1 = tbl{:, feature(f,1)};
    f2 = tbl{:, feature(f,2)};
    
    figure; plot(f1, f2, '.'); grid on; 
    xlabel(f1name); ylabel(f2name); 
    text(f1, f2, lbl)
end