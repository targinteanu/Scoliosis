%function M = writheSphere2(patient_number)

[num, txt] = xlsread('Writhe-pre-post_new-metrics.csv');
N = 32;
num = num(1:N, :);
XYZ = num(:, 13:end); 

    x = XYZ(patient_number, 1:3:51); 
    y = XYZ(patient_number, 2:3:51); 
    z = XYZ(patient_number, 3:3:51); 
    
cm = [x; y; z]';
[sx,sy,sz] = sphere;

views = [0 0; 0 90];
numplots = 2;
cursori = cell(1,numplots); cursorj = cell(1,numplots);
Rvector = cell(4, numplots);
figure('Position', [50 300 1400 525]); 
for idx = fliplr(1:numplots)
    subplot(1,numplots,idx); 
    plot3(x/500, y/500, z/500 - 1, 'LineWidth', 3);
    hold on; mesh(sx, sy, sz, zeros(size(sx))); hidden off; grid on;
    view(views(idx,:));
    cursori{idx} = plot3(x([1,2])/500, y([1,2])/500, z([1,2])/500 - 1, '*k', 'LineWidth', 1.1);
    cursorj{idx} = plot3(x([1,2])/500, y([1,2])/500, z([1,2])/500 - 1, 'ok', 'LineWidth', 1.1);
    for iidx = 1:4
        Rvector{iidx, idx} = plot3([1 0], [0 1], [0 0], 'k', 'LineWidth', 1.1);
    end
end
grid off;

dirold = [0 0 -1];

%%
range = 2:length(cm);
        frames = 20;
        T = .1;
M = cell(1, frames*((length(cm)-1)^2));
for i = range
    for j = range
        
        for idx = 1:numplots
            cursori{idx}.XData = x([i-1,i])/500; 
            cursori{idx}.YData = y([i-1,i])/500; 
            cursori{idx}.ZData = z([i-1,i])/500 - 1;
            cursorj{idx}.XData = x([j-1,j])/500; 
            cursorj{idx}.YData = y([j-1,j])/500; 
            cursorj{idx}.ZData = z([j-1,j])/500 - 1;
        end
        
        if abs(i-j) > 2
        
        p1 = cm(i-1,:); p2 = cm(i,:); 
        p3 = cm(j-1,:); p4 = cm(j,:);
        
        r34 = p4-p3; r12 = p2-p1;
        r13 = p3-p1; r14 = p4-p1; r23 = p3-p2; r24 = p4-p2;
        sgn = sign( (cross(r34,r12)) *r13');
                
        r13 = r13/norm(r13); r14 = r14/norm(r14); r23 = r23/norm(r23); r24 = r24/norm(r24);
        
        %view(r13);
        bigAngle = acos(r13*dirold');
        smallAngle = bigAngle/frames;
        for f = 1:frames
            newview = changeviews(dirold, r13, f*smallAngle, bigAngle);
            if norm(newview)
                view(newview);
            end
            pause(T/frames);
            
        
            F = getframe(gcf);
            M{ ((i-1)*(length(cm)-1) + j)*(f-1)*frames + f } = F.cdata;
        end
        
        R = [r13; r14; r24; r23];
        px = R(:,1); py = R(:,2); pz = R(:,3);
        colr = {'red', 'blue', 'green'};
        for idx = fliplr(1:numplots)
            subplot(1,numplots,idx);
            patch(px, py, pz, colr{sgn + 2});
            for iidx = 1:2
                Rvector{iidx,idx}.XData = [0, px(iidx)] + x(i-1)/500;
                Rvector{iidx,idx}.YData = [0, py(iidx)] + y(i-1)/500;
                Rvector{iidx,idx}.ZData = [0, pz(iidx)] + z(i-1)/500 - 1;
            end
            for iidx = 3:4
                Rvector{iidx,idx}.XData = [0, px(iidx)] + x(i)/500;
                Rvector{iidx,idx}.YData = [0, py(iidx)] + y(i)/500;
                Rvector{iidx,idx}.ZData = [0, pz(iidx)] + z(i)/500 - 1;
            end
        end
        
        dirold = r13;
        
        pause(.25);
        
        end
    end
end

function view3 = changeviews(view1, view2, angle13, angle12)
    %{
    A = rref([view1; view2]);
    a = A(1,3); b = A(2,3); 
    syms x y z;
    eqn1 = a*x + b*y == z; % view3 in plane spanned by view1,2
    eqn2 = view1*[x;y;z] == cos(angle); % angle is between view1,3
    eqn3 = acos(view1*[x;y;z]) + acos(view2*[x;y;z]) == acos(view1*view2'); % v3 is between v1,v2
    eqn4 = norm([x;y;z]) == 1; % view3 is unit
    [X, Y, Z] = solve([eqn1, eqn2, eqn3, eqn4], [x, y, z]);
    view3 = [X, Y, Z];
    %}
    %ANG = acos(view1*view2');
    A = rref([view1*view1', view1*view2', cos(angle13); view1*view2', view2*view2', cos(angle12-angle13)]);
    view3 = A(1,3)*view1 + A(2,3)*view2;
end

%end