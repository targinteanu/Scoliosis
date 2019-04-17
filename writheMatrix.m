load('spines_XYZ.mat')
idx = 20;

    xyz = spinesXYZ{idx, 2:end};
    x = xyz(1:3:end); 
    y = xyz(2:3:end); 
    z = xyz(3:3:end); 
    
cm = [x; y; z]';

figure; plot3(x,y,z,'-o'); grid on; view([90 0])
title(num2str(idx));

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

zerolim = 5e-4; M = M.*(abs(M) > zerolim);
figure; heatmap(sign(M));
%figure; heatmap(M);
title(num2str(idx));

%%
writheTop = zeros(1, 17); writheBot = writheTop; writheCro = writheTop;
for p = 2:(length(x)-1)
    writheTop(p) = levittWrithe(cm, 1:p); 
    writheBot(p) = levittWrithe(cm, p:17);
    writheCro(p) = (levittWrithe(cm) - writheTop(p) - writheBot(p))/2;
end
%figure; plot([writheTop; writheBot; writheCro]')
figure; plot(writheTop + writheBot);
title(num2str(idx));

%%
writheTop2 = zeros(17, 17); writheBot2 = writheTop2; writheMid2 = writheTop2; writheCro2 = writheTop2;
for p = 2:(length(x)-1)
    for q = p:(length(x)-1)
        writheTop2(p,q) = levittWrithe(cm, 1:p); 
        writheMid2(p,q) = levittWrithe(cm, p:q);
        writheBot2(p,q) = levittWrithe(cm, q:17);
        writheCro2(p,q) = (levittWrithe(cm) - writheTop2(p,q) - writheBot2(p,q) - writheMid2(p,q))/2;
    end
end
%figure; plot([writheTop2; writheBot2; writheMid2]')
figure; heatmap(abs(writheTop2 + writheMid2 + writheBot2));
title(num2str(idx));
