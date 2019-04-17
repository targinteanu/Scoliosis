function M = writheSphere(patient_number)

load('spines_XYZ.mat')

    xyz = spinesXYZ{patient_number, 2:end};
    x = xyz(1:3:end); 
    y = xyz(2:3:end); 
    z = xyz(3:3:end); 
    
cm = [x; y; z]';
[sx,sy,sz] = sphere;

views = [0 0; 90 0; 0 90];
cursori = cell(1,3); cursorj = cell(1,3);
figure('Position', [50 500 1400 350]); 
for idx = 1:3
    subplot(1,3,idx); 
    plot3(x/500, y/500, z/500 - 1, 'LineWidth', 3);
    hold on; mesh(sx, sy, sz, zeros(size(sx))); hidden off; grid on;
    view(views(idx,:));
    cursori{idx} = plot3(x(2)/500, y(2)/500, z(2)/500 - 1, '*k', 'LineWidth', 1.1);
    cursorj{idx} = plot3(x(2)/500, y(2)/500, z(2)/500 - 1, 'ok', 'LineWidth', 1.1);
end

%%
range = 2:length(cm);
M = cell(1, (length(cm)-1)^2);
for i = range
    for j = range
        
        for idx = 1:3
            cursori{idx}.XData = x(i)/500; cursori{idx}.YData = y(i)/500; 
            cursori{idx}.ZData = z(i)/500 - 1;
            cursorj{idx}.XData = x(j)/500; cursorj{idx}.YData = y(j)/500; 
            cursorj{idx}.ZData = z(j)/500 - 1;
        end
        
        if abs(i-j) > 2
        
        p1 = cm(i-1,:); p2 = cm(i,:); 
        p3 = cm(j-1,:); p4 = cm(j,:);
        
        r34 = p4-p3; r12 = p2-p1;
        r13 = p3-p1; r14 = p4-p1; r23 = p3-p2; r24 = p4-p2;
        sgn = sign( (cross(r34,r12)) *r13');
        
        r13 = r13/norm(r13); r14 = r14/norm(r14); r23 = r23/norm(r23); r24 = r24/norm(r24);
        
        R = [r13; r14; r23; r24];
        px = R(:,1); py = R(:,2); pz = R(:,3);
        colr = {'red', 'blue', 'green'};
        for idx = 1:3
            subplot(1,3,idx);
            patch(px, py, pz, colr{sgn + 2});
        end
        
        F = getframe(gcf);
        M{(i-1)*(length(cm)-1) + j} = F.cdata;
        
        pause(.1);
        
        end
    end
end

end