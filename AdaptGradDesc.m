function args = AdaptGradDesc(gradfun, yfun, Xtrue, Ytrue, b, mu)
% gradfun: of the form d_err/d_b = gradfun(b_old)
% yfun: of the form y = yfun(b, x)
% Xtrue, Ytrue: paired actual vaules of Y, X used to calculate error 
% b: initial coefficients 
% mu: initial learning rate

% init vars
nsteps = 10000;
errs = zeros(1,nsteps); derrs = zeros(size(errs)); gradnorms = zeros(size(errs));
mus = zeros(1,nsteps); mus(1) = mu;
bs = zeros(length(b),nsteps); bs(:,1)=b;
derr = 0;  
KP = 1e-10;

% initial gradient, error calculations
err = sum( arrayfun(@(j) Ytrue(j) - yfun(b, Xtrue(j)), 1:length(Xtrue)) );
derr_db = gradfun(b); gradnorm = norm(derr_db); 
errs(1) = err; gradnorms(1) = gradnorm;

% init plots
colr = {'k',  'b',  'r',  'c',  'm',  'g',  'y'};
styl = {'ok', 'ok', 'or', 'oc', 'om', 'og', 'oy'};
figure('Position', [50 100 1300 700]);
plts(1).ax = subplot(1,3,1); plts(1).t = ['Error = ' num2str(err)];
plts(1).plt{1} = plot(errs, colr{1}); hold on; plts(1).plt{2} = plot(derrs, colr{2}); 
    plts(1).plt{3} = plot(gradnorms, colr{3});
plts(1).nextval{1} = plot(2,err,styl{1}); plts(1).nextval{2} = plot(2,derr,styl{2}); 
    plts(1).nextval{3} = plot(2,err,styl{3});
plts(2).ax = subplot(1,3,2); plts(1).t = ['\mu = ' num2str(mu)];
plts(2).plt{1} = plot(mus, colr{1});
plts(2).nextval{1} = plot(2,mu,styl{1});
plts(2).ax = subplot(1,3,2); plts(1).t = ['b0 = '];
for j = 1:length(b)
    plts(3).plt{j} = plot(bs(j,:), colr{j}); hold on; 
    plts(3).nextval{j} = plot(2,b(j),styl{j});
end
xlim(plts.ax, [0, 10]);
pause(.01);


% iterate
for n = 2:nsteps
    retry = true;
    while retry
        retry = false;
    
        % try performing a gradient step and getting the new error
        b_new = b + mu*derr_db;
        try
            err = sum( arrayfun(@(j) Ytrue(j) - yfun(b_new, Xtrue(j)), 1:length(Xtrue)) );
        catch
            retry = true;
        end
        
        if ~retry
            derr = err - errs(n-1);
            
            if n > 10
                % ADJUST LEARNING RATE
                mu = mu - KP*derr;
                retry = retry|(derr>0);
            else
                retry = false;
            end
        end
        
        if ~retry
            b = b_new;
            try
                derr_db = gradfun(b); 
            catch
                retry = true;
            end
            
            if ~retry
                gradnorm = norm(derr_db);
                
                bs(:,n) = b;
                errs(n) = err; derrs(n) = derr; gradnorms(n) = gradnorm;
                mus(n) = mu;
            end
        end
        
        % update plots ------------------------------------------------
        
        % update next val indicator
        for hp = plts
            for j = 1:length(hp.nextval)
                hp.nextval{j}.XData = n;
            end
        end
        plts(1).nextval{1}.YData = err; plts(1).nextval{2}.YData = derr; plts(1).nextval{3}.YData = gradnorm;
        plts(2).nextval{1}.YData = mu;
        for j = 1:length(b_new)
            plts(3).nextval{j}.YData = b_new(j);
        end
        
        % update plots
        if ~retry
            plts(1).plt{1}.YData = errs; plts(1).plt{2}.YData = derrs; plts(1).plt{3}.YData = gradnorms;
            plts(2).plt{1}.YData = mus;
            for j = 1:length(b_new)
                plts(3).plt{j}.YData = bs(j);
            end
        end
        
        % update title, view range
        for hp = plts
            hp.t.String = num2str(hp.nextval{1}.YData(1));
            win = [max(0, n-100), max(n, 10)];
            xlim(hp.ax, win);
            mn = zeros(length(hp.plt)); mx = zeros(size(mn));
            for j = 1:length(hp.plt)
                mn(j) = min(hp.nextval{j}.YData, min(hp.plt{j}.YData(win(1):win(2))));
                mx(j) = max(hp.nextval{j}.YData, max(hp.plt{j}.YData(win(1):win(2))));
            end
            mn = min(mn); mx = max(mx);
            if mx > mn
                ylim(hp.ax, [mn,mx] + [-1,1]*.25*(mx-mn));
            end
        end
    
    end
    
end


end