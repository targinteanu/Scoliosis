% gets data from XL spreadsheet and draws specified plots 

%% load data

vertnames = [arrayfun(@(i) ['T' num2str(i)], 1:12, 'UniformOutput', 0), ...
    arrayfun(@(i) ['L' num2str(i)], 1:5, 'UniformOutput', 0)];

[num, txt] = xlsread('Writhe-pre-post_new-metrics.xlsx');

XYZ = num(:, 15:65); 
ROT = num(:, 66:end);

X = XYZ(:, 1:3:51);
Y = XYZ(:, 2:3:51);
Z = XYZ(:, 3:3:51); 
uX = cos(ROT); uY = sin(ROT); 

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

% put variables in matrix so they can be compared/plotted in a loop

vars = [shapecluster, writhe_postop, writhe_preop, abswrithe, twist, ...
    globtorsion, vertNeutral, vertApical, ...
    maxtorsion1, maxtor1vert, maxtorsion2, maxtor2vert, ...
    abs(writhe_postop), abs(writhe_preop), abs(twist), ...
    abs(globtorsion), abs(maxtorsion1), abs(maxtorsion2), ...
    ];

varnames = {...
    'Shape Cluster', 'Writhe Post-Op', 'Writhe', 'Absolute Writhe', 'Twist', ...
    'Kadoury Global Torsion', 'Neutral Vertebra', 'Apical Vertebra', ...
    'Max Torsion', 'Max Torsion Vertebra', ...
    '2nd-Max Torsion', '2nd-Max Torsion Vertebra', ...
    '|Writhe Post-Op|', '|Writhe|', '|Twist|', ...
    '|Kadoury Global Torsion|', '|Max Torsion|', '|2nd-Max Torsion|', ...
    };

%% plots 

feature = [3, 5; 14, 15; 3, 6; 3, 9; 3, 16; 3, 17; 14, 16; 14, 17];
% Each row is a pair of features to be compared. Feature numbers correspond
% to viariables in the matrix vars with descriptions given by varnames. 

for f = 1:length(feature(:,1))
    
    % set up features for plots
    
    f1name = varnames{feature(f,1)};
    f2name = varnames{feature(f,2)};
    f1 = vars(:, feature(f,1));
    f2 = vars(:, feature(f,2));
    
    %{
    % plot feature 1 vs feature 2 with patient numbers 
    figure; plot(f1, f2, '.'); grid on; 
    xlabel(f1name); ylabel(f2name); 
    text(f1, f2, lbl)
    %}
    
    % patient curve plots 
    
    % separate based on low feature 1 / high feature 2 and low feature 2 /
    % high feature 1. Low/High is defined as lower/higher than the mean of
    % that feature. 
    lo1 = f1 < mean(f1); hi1 = f1 > mean(f1); 
    lo2 = f2 < mean(f2); hi2 = f2 > mean(f2); 
    lo1hi2 = lo1&hi2; hi1lo2 = hi1&lo2; 
    % group a = low feat 1 / high feat 2; group b = vice versa 
    Xa = (X(lo1hi2,:)); Ya = (Y(lo1hi2,:)); Za = (Z(lo1hi2,:));
    uXa = (uX(lo1hi2,:)); uYa = (uY(lo1hi2,:)); 
    Xb = (X(hi1lo2,:)); Yb = (Y(hi1lo2,:)); Zb = (Z(hi1lo2,:));
    uXb = (uX(hi1lo2,:)); uYb = (uY(hi1lo2,:)); 
    
    a = ['low ' f1name ' high ' f2name]; b = ['low ' f2name ' high ' f1name]; 
    figure('Position', [50 50 1400 850]);
    
    % plot the average spines (averaged x, y, and z points) for group a in
    % 3 views with 1 standard deviation error bars. 
    subplot(231);  
    errorbar(mean(Ya), mean(Za), std(Za), std(Za), std(Ya), std(Ya)); 
    xlabel('y'); ylabel('z'); grid on;
    subplot(232);  
    errorbar(mean(Xa), mean(Za), std(Za), std(Za), std(Xa), std(Xa)); 
    xlabel('x'); ylabel('z'); grid on;
    title(a); 
    subplot(233);  
    plot3(mean(Xa), mean(Ya), mean(Za), '-o'); 
    xlabel('x'); ylabel('y'); zlabel('z'); grid on;
    view([0,90]);
    
    % repeat for group b 
    subplot(234);  
    errorbar(mean(Yb), mean(Zb), std(Zb), std(Zb), std(Yb), std(Yb)); 
    xlabel('y'); ylabel('z'); grid on;
    subplot(235);  
    errorbar(mean(Xb), mean(Zb), std(Zb), std(Zb), std(Xb), std(Xb)); 
    xlabel('x'); ylabel('z'); grid on;
    title(b); 
    subplot(236);  
    plot3(mean(Xb), mean(Yb), mean(Zb), '-o'); 
    xlabel('x'); ylabel('y'); zlabel('z'); grid on;
    view([0,90]);
end