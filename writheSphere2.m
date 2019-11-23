%function M = writheSphere2(patient_number)

[num, txt] = xlsread('Writhe-pre-post_new-metrics.csv');
N = 32;
num = num(1:N, :);
XYZ = num(:, 13:end); 

    x = XYZ(patient_number, 1:3:51); x0=x;
    y = XYZ(patient_number, 2:3:51); y0=y;
    z = XYZ(patient_number, 3:3:51); z0=z;
    
cm = [x; y; z]';
[sx,sy,sz] = sphere;

spscl = 750;
views = [0 0; 0 90];
numplots = 2;
cursori = cell(1,numplots); cursorj = cell(1,numplots);
spn = cell(1, numplots);
Rvector = cell(4, numplots);
figure('Position', [50 300 1400 525]); 
for idx = fliplr(1:numplots)
    subplot(1,numplots,idx); 
    spn{idx} = plot3(x/spscl, y/spscl, z/spscl - 1, 'LineWidth', 3);
    hold on; mesh(sx, sy, sz, zeros(size(sx))); hidden off; grid on;
    view(views(idx,:));
    cursori{idx} = plot3(x([1,2])/spscl, y([1,2])/spscl, z([1,2])/spscl - 1, '*k', 'LineWidth', 1.1);
    cursorj{idx} = plot3(x([1,2])/spscl, y([1,2])/spscl, z([1,2])/spscl - 1, 'ok', 'LineWidth', 1.1);
    for iidx = 1:4
        Rvector{iidx, idx} = plot3([1 0], [0 1], [0 0], '--k', 'LineWidth', 1.1);
    end
end
grid off;
subplot(1,numplots,1); xlim([-1 1]); ylim([-1 1]); zlim([-2 2]);

dirold = [0 0 -1];

%%
range = 2:length(cm);
        frames = 20;
        T = .001; % pause between frames
        T2 = .15; % pause between plots 
M = cell(1, frames*((length(cm)-1)^2));
for i = range
    for j = range
        
        if abs(i-j) > 2
        %% get writhe stuff 
        p1 = cm(i-1,:); p2 = cm(i,:); 
        p3 = cm(j-1,:); p4 = cm(j,:);
        
        r34 = p4-p3; r12 = p2-p1;
        r13 = p3-p1; r14 = p4-p1; r23 = p3-p2; r24 = p4-p2;
        sgn = sign( (cross(r34,r12)) *r13');
                
        r13 = r13/norm(r13); r14 = r14/norm(r14); r23 = r23/norm(r23); r24 = r24/norm(r24);
        
        %% update view 
        %view(r13);
        %dir = r13*[0,0,1;0,1,0;-1,0,0]; 
        dir = r13*[cos(pi/4), 0, sin(pi/4); 0,1,0; -sin(pi/4), 0, cos(pi/4)]*...
            [cos(pi/4), -sin(pi/4), 0; sin(pi/4), cos(pi/4), 0; 0,0,1];
        dir = dir/norm(dir);
        bigAngle = acos(dir*dirold');
        smallAngle = bigAngle/frames;
        %%{
        for f = 1:frames
            newview = changeviews(dirold, dir, f*smallAngle, bigAngle);
            if norm(newview)
                view(newview);
            end
            pause(T/frames);
            
        
            F = getframe(gcf);
            M{ ((i-1)*(length(cm)-1) + j)*(f-1)*frames + f } = F.cdata;
        end
        %}
        
        dirold = dir;
        
        %% update plots 
        x = x0-p1(1); y = y0-p1(2); z = z0-p1(3); % center on p1
        for idx = 1:numplots
            spn{idx}.XData = x/spscl; spn{idx}.YData = y/spscl; spn{idx}.ZData = z/spscl;
            cursori{idx}.XData = x([i-1,i])/spscl; 
            cursori{idx}.YData = y([i-1,i])/spscl; 
            cursori{idx}.ZData = z([i-1,i])/spscl;
            cursorj{idx}.XData = x([j-1,j])/spscl; 
            cursorj{idx}.YData = y([j-1,j])/spscl; 
            cursorj{idx}.ZData = z([j-1,j])/spscl;
        end
        
        %R = [r13; r14; r24; r23] + [p3;p4;p4;p3]/spscl - [zeros(4,2), ones(4,1)];
        R = [r13; r14; r24; r23];
        px = R(:,1); py = R(:,2); pz = R(:,3);
        p1=p1-p1; p2=p2-p1; p3=p3-p1; p4=p4-p1;% center on p1
        colr = {'red', 'blue', 'green'};
        for idx = fliplr(1:numplots)
            subplot(1,numplots,idx);
            patch(px, py, pz, colr{sgn + 2});
            
            Rvector{1,idx}.XData = [0, r13(1)]; 
            Rvector{1,idx}.YData = [0, r13(2)];
            Rvector{1,idx}.ZData = [0, r13(3)];
            
            Rvector{2,idx}.XData = [0, r14(1)]; 
            Rvector{2,idx}.YData = [0, r14(2)];
            Rvector{2,idx}.ZData = [0, r14(3)];
            
            Rvector{3,idx}.XData = [0, r23(1)]; 
            Rvector{3,idx}.YData = [0, r23(2)];
            Rvector{3,idx}.ZData = [0, r23(3)];
            
            Rvector{4,idx}.XData = [0, r24(1)]; 
            Rvector{4,idx}.YData = [0, r24(2)];
            Rvector{4,idx}.ZData = [0, r24(3)];
            %{
            Rvector{1,idx}.XData = [p1(1), p4(1), p4(1)]/spscl + [0,0,r14(1)];
            Rvector{1,idx}.YData = [p1(2), p4(2), p4(2)]/spscl + [0,0,r14(2)];
            Rvector{1,idx}.ZData = [p1(3), p4(3), p4(3)]/spscl - 1 + [0,0,r14(3)];
            
            Rvector{2,idx}.XData = [p1(1), p3(1), p3(1)]/spscl + [0,0,r13(1)];
            Rvector{2,idx}.YData = [p1(2), p3(2), p3(2)]/spscl + [0,0,r13(2)];
            Rvector{2,idx}.ZData = [p1(3), p3(3), p3(3)]/spscl - 1 + [0,0,r13(3)];
            
            Rvector{3,idx}.XData = [p2(1), p3(1), p3(1)]/spscl + [0,0,r23(1)];
            Rvector{3,idx}.YData = [p2(2), p3(2), p3(2)]/spscl + [0,0,r23(2)];
            Rvector{3,idx}.ZData = [p2(3), p3(3), p3(3)]/spscl - 1 + [0,0,r23(3)];
            
            Rvector{4,idx}.XData = [p2(1), p4(1), p4(1)]/spscl + [0,0,r24(1)];
            Rvector{4,idx}.YData = [p2(2), p4(2), p4(2)]/spscl + [0,0,r24(2)];
            Rvector{4,idx}.ZData = [p2(3), p4(3), p4(3)]/spscl - 1 + [0,0,r24(3)];
            %}
            %{
            for iidx = 1:2
                Rvector{iidx,idx}.XData = [x(i-1)/spscl, px(iidx)];
                Rvector{iidx,idx}.YData = [y(i-1)/spscl, py(iidx)];
                Rvector{iidx,idx}.ZData = [z(i-1)/spscl - 1, pz(iidx)];
            end
            for iidx = 3:4
                Rvector{iidx,idx}.XData = [x(i)/spscl, px(iidx)];
                Rvector{iidx,idx}.YData = [y(i)/spscl, py(iidx)];
                Rvector{iidx,idx}.ZData = [z(i)/spscl - 1, pz(iidx)];
            end
            %}
        end
                
        pause(T2);
        
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