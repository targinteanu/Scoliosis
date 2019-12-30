N = 5; cm = cell(1,N); theta = cell(1,N);
npts = 17;
amax = 180;
lbl = {'a) ', 'b) ', 'c) ', 'd) ', 'e) '};

s = linspace(0, 2*pi*sqrt(2), npts)';

%cm{1} = [zeros(size(s)), zeros(size(s)), s];
cm{1} = [zeros(size(s)), cos(s/sqrt(2)), s/sqrt(2)];
theta{1} = 90 + zeros(size(s));
cm{2} = cm{1};
theta{2} = 90 + linspace(0, amax, npts)';
cm{3} = [sin(s/sqrt(2)), cos(s/sqrt(2)), s/sqrt(2)];
theta{3} = theta{1};
cm{4} = cm{3};
theta{4} = theta{2};
cm{5} = [zeros(size(s)), zeros(size(s)), s/sqrt(2)];
theta{5} = theta{2};

%figure('Position', [50 50 1400 500]); 
figure('Position', [0 50 1500 500]);
for i = 1:N
    subplot(1,N,i);
    plot3dSpine(cm{i}, theta{i}); 
    Tw = getTwist(cm{i}, theta{i});
    Wr = levittWrithe(cm{i});
    title([lbl{i}, 'Twist = ' num2str(Tw) ', Writhe = ' num2str(Wr)]);
    view([-30, 10]);
    xlim([-2, 2]); ylim([-2, 2]); zlim([0, 2*pi]);
end

%%
%%{
a = 1:2000; Tw = zeros(size(a)); Twact = Tw;
for ang = a
  
%CM = [zeros(size(s)), zeros(size(s)), s];
CM = cm{2};
th = 90 + linspace(0, ang, npts)';

%i = 5; %figure; 
%    plot3dSpine(cm{i}, theta{i}); 
    %Tw(ang) = getTwist(cm{i}, cm{i}+[cos(theta{i}), sin(theta{i}), zeros(size(theta{i}))]);
    Tw(ang) = getTwist(CM, th);
    
    Twact(ang) = -ang/(360);
end
figure; plot(a/360, [Tw; Twact]); %xticks(0:45:720);
%D = cm{i}(2:end,:) - cm{i}(1:(end-1),:);
%DD = D*D'; DD = diag(DD);
%mean(DD)
%}

%% writhe by hand 
r = @(t) [cos(t), sin(t), t]; dr = @(t) [-sin(t), cos(t), 1];
q = @(t) [cos(t), -sin(t), t]; dq = @(t) [-sin(t), -cos(t), 1];

ti = 0; tf = 10; dt = .2;
T = ti:dt:tf; T = T';
R = r(T); Q = q(T); 

wr_r = integral2( @(t1,t2) ddWr(t1,t2, r, r, dr, dr), ...
    ti, tf, ti, tf, 'Method', 'iterated');
wr_q = integral2( @(t1,t2) ddWr(t1,t2, q, q, dq, dq), ...
    ti, tf, ti, tf, 'Method', 'iterated');
wr_r = wr_r/(4*pi); wr_q = wr_q/(4*pi);
wr_R = levittWrithe(R); wr_Q = levittWrithe(Q);

figure; 
subplot(1,2,1); plot3(R(:,1), R(:,2), R(:,3), '-o'); grid on; 
title([num2str(wr_r), ' | ', num2str(wr_R)]);
xlabel('x'); ylabel('y'); zlabel('z'); view([20, 10]);
subplot(1,2,2); plot3(Q(:,1), Q(:,2), Q(:,3), '-o'); grid on; 
title([num2str(wr_q), ' | ', num2str(wr_Q)]);
xlabel('x'); ylabel('y'); zlabel('z'); view([20, 10]);

%% twist by hand
proj = @(u, x) ((u*x')/(x*x'))*x;
r = @(s) [sin(s/sqrt(2)), cos(s/sqrt(2)), s/sqrt(2)];
dr = @(s) (1/sqrt(2))*[cos(s/sqrt(2)), -sin(s/sqrt(2)), 1];
q = @(s) [sin(s/sqrt(2)), -cos(s/sqrt(2)), s/sqrt(2)];
dq = @(s) (1/sqrt(2))*[cos(s/sqrt(2)), sin(s/sqrt(2)), 1];
v = @(s) [0, 1, 0];
vr = @(s) v(s) - proj(v(s), r(s)); ur = @(s) vr(s)/norm(vr(s));
vq = @(s) v(s) - proj(v(s), q(s)); uq = @(s) vq(s)/norm(vq(s));

R = r(s); Q = q(s); 
%VR = vr(s); VQ = vq(s);
VR = theta{3}; VQ = VR;

tw_r = integral( @(t) dTw(t, dr, ur), s(1), s(end));%, 'Method', 'iterated');
tw_q = integral( @(t) dTw(t, dq, uq), s(1), s(end));%, 'Method', 'iterated');
tw_r = tw_r/(2*pi); tw_q = tw_q/(2*pi);
tw_R = getTwist(R, VR); tw_Q = getTwist(Q, VQ);

figure; 
subplot(1,2,1); plot3dSpine(R, VR); view([-30,10]);
xlim([-2, 2]); ylim([-2, 2]); zlim([0, 2*pi]);
title([num2str(tw_r), ' | ', num2str(tw_R)]);
subplot(1,2,2); plot3dSpine(Q, VQ); view([-30,10]);
xlim([-2, 2]); ylim([-2, 2]); zlim([0, 2*pi]);
title([num2str(tw_q), ' | ', num2str(tw_Q)]);

%% twist and writhe functions 
    function dTw_ = dTw(t, dr, u)
    % u, du must be entered s.t. u perp dr
        t = t';
        dTw_ = zeros(size(t)); imax = length(t);
        for i = 1:imax
            tt = t(i);
            if (i == 1) 
                if (i+1 <= imax)
                    dtt = t(i+1)-t(i);
                else
                    dtt = .0001;
                end
            else
                dtt = t(i)-t(i-1);
            end
            du = (u(tt + dtt) - u(tt))/dtt;
            dTw_(i) = cross(du, u(tt)) * dr(tt)';
        end
        dTw_ = dTw_';
    end

    function ddWr_ = ddWr(t1, t2, r1, r2, dr1, dr2)
        t1 = t1'; t2 = t2';
        r12 = r1(t1) - r2(t2);
        if length(t1) <= length(t2)
            sz = size(t1); imax = length(t1);
        else
            sz = size(t2); imax = length(t2);
        end
        ddWr_ = zeros(sz);
        for i = 1:imax
            tt1 = t1(i); tt2 = t2(i); rr12 = r12(i,:);
            if norm(rr12)
                ddWr_(i) = (cross(dr1(tt1), dr2(tt2)) * rr12')/( norm(rr12)^3 );
            else
                ddWr_(i) = 0;
            end
        end
        ddWr_ = ddWr_';
        %figure; plot3(t1, t2, ddWr_, '.'); grid on; xlabel('t1'); ylabel('t2');
    end