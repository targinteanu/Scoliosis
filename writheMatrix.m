load('spines_XYZ.mat')
idx = 14;

    xyz = spinesXYZ{idx, 2:end};
    x = xyz(1:3:end); 
    y = xyz(2:3:end); 
    z = xyz(3:3:end); 
    
cm = [x; y; z]';

figure; plot3(x,y,z,'-o'); grid on;

%%
range = 2:length(cm);
M = zeros(length(cm));

for i = range
    for j = range
        
        if abs(i-j) > 2
        
        p1 = cm(i-1,:); p2 = cm(i,:); 
        p3 = cm(j-1,:); p4 = cm(j,:);
        
        r34 = p4-p3; r12 = p2-p1;
        r13 = p3-p1; r14 = p4-p1; r23 = p3-p2; r24 = p4-p2;
        
        n = zeros(4, 3);
        n(1,:) = cross(r13, r14);
        n(2,:) = cross(r14, r24);
        n(3,:) = cross(r24, r23);
        n(4,:) = cross(r23, r13);
        
        for k = 1:4
            if norm(n(k,:))
                % normalize 
                n(k,:) = n(k,:)/norm(n(k,:));
            end
        end
        
        N = n([2:end 1],:);
        as = diag(n*N'); angles = asin(as) + pi/2;

        Omega = sum(angles) - 2*pi;
        Omega = Omega*sign( (cross(r34,r12)) *r13');
        
        M(i,j) = Omega;
        end
    end
end

figure; heatmap(M);

%%
writheTop = zeros(1, 17); writheBot = writheTop; writheCro = writheTop;
for p = 2:(length(x)-1)
    writheTop(p) = levittWrithe(cm, 1:p); 
    writheBot(p) = levittWrithe(cm, p:17);
    writheCro(p) = (levittWrithe(cm) - writheTop(p) - writheBot(p))/2;
end
figure; plot([writheTop; writheBot; writheCro]')