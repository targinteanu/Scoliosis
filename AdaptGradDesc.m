function beta = AdaptGradDesc(gradfun, yfun, Xtrue, Ytrue, b, mu, PID, p2s)
% gradfun: of the form d_err/d_b = gradfun(b_old, solvedParams)
% yfun: of the form y = yfun(b, x, solvedParams)
% Xtrue, Ytrue: paired actual vaules of Y, X used to calculate error 
% b: initial coefficients 
% mu: initial learning rate
% PID: [KP, KI, KD]
% p2s (paramsToSolve): function of b that returns 2xm cell array / solvedVars = fsolve(eqn, Vars)
%                Row 1: sys of equations of the form 0 = funs(Vars, b, solvedParams)
%                Row 2: starting Vars0 points corresponding to equations 

% init vars
nsteps = 3000;
errs = nan(nsteps,1); derrs = nan(size(errs)); gradnorms = nan(size(errs));
mus = nan(nsteps,1); mus(1) = mu;
bs = nan(nsteps,length(b)); bs(1,:)=b;
derr = 0; 
KP = PID(1); KI = PID(2);

% initial gradient, error calculations
options = optimoptions('fsolve','Display','none');
paramsToSolve = p2s(b);
solvedParams = cell(size(paramsToSolve,2),1);
for j = 1:length(paramsToSolve)
    SP = fsolve(@(Vars) paramsToSolve{1,j}(Vars, b, solvedParams), ...
        paramsToSolve{2,j}, options);
    solvedParams{j} = SP;
end
err = sum( arrayfun(@(j) (Ytrue(j) - yfun(b, Xtrue(j), solvedParams)).^2, 1:length(Xtrue)) );
derr_db = gradfun(b, solvedParams); gradnorm = norm(derr_db); 
errs(1) = err; gradnorms(1) = gradnorm;

% init plots
colr = {'k',  'b',  'r',  'c',  'm',  'g',  'y'};
styl = {'ok', 'ok', 'or', 'oc', 'om', 'og', 'oy'};
figure('Position', [50 100 1300 700]);
plts(1).ax = subplot(1,3,1); 
plts(1).plt{1} = plot(errs, colr{1}); hold on; plts(1).plt{2} = plot(derrs, colr{2}); 
%    plts(1).plt{3} = plot(gradnorms, colr{3});
plts(1).nextval{1} = plot(2,err,styl{1}); plts(1).nextval{2} = plot(2,derr,styl{2}); 
%    plts(1).nextval{3} = plot(2,err,styl{3});
plts(1).t = title(['Error = ' num2str(err)]);
plts(2).ax = subplot(1,3,2); 
plts(2).plt{1} = plot(mus, colr{1}); hold on;
plts(2).nextval{1} = plot(2,mu,styl{1});
plts(2).t = title(['\mu = ' num2str(mu)]);
plts(3).ax = subplot(1,3,3); 
for j = 1:length(b)
    jj = mod(j,length(colr)) + 1;
    plts(3).plt{j} = plot(bs(:,j), colr{jj}); hold on; 
    plts(3).nextval{j} = plot(2,b(j),styl{jj});
end
plts(3).t = title(['b0 = ']);
%{
for hp = plts(1:2)
    xlim(hp.ax, [0, 10]);
end
%}
pause(1);


% iterate
n = 1;
notyetfound = true;
mulimit = 1e-100; enablemulimit = false;
while notyetfound
    n = n+1;
    retry = true;
    while retry&notyetfound
        retry = false;
        enablemulimit = enablemulimit|(abs(mu)>mulimit);
    
        % try performing a gradient step and getting the new error
        b_new = b - mu*derr_db;
        
        try
            paramsToSolve = p2s(b_new);
            for j = 1:length(paramsToSolve)
                solvedParams{j} = fsolve(@(Vars) paramsToSolve{1,j}(Vars, b_new, solvedParams), ...
                    paramsToSolve{2,j}, options);
            end
        catch
            retry = true;
            mu = mu/2;
        end
        
        if ~retry
            err = sum( arrayfun(@(j) (Ytrue(j) - yfun(b_new, Xtrue(j), solvedParams)).^2, 1:length(Xtrue)) );
            
            derr = err - errs(n-1); 
            db = norm(b_new-b);
            if db
                derr=derr/db;
            end
            
            if n > 5
                if n > 7
                    notyetfound = (n<nsteps);%&((abs(derr)>1e-6)|retry);
                    if enablemulimit&(abs(mu)<mulimit)
                        notyetfound = false;
                    end
                end
                
                % ADJUST LEARNING RATE
                PIDvar = -KP*derr +KI*err;
                if (derr>(0*err))|(sum([b_new(:);err])==Inf)
                    retry = true;
                    mu = .5*mu;
                else
                    mu = mu + PIDvar;
                end
            else
                retry = false;
            end
        end
        
        if ~retry
            b = b_new;
            try
                derr_db = gradfun(b, solvedParams); 
            catch
                retry = true;
                mu = mu/2;
            end
            
            if ~retry
                gradnorm = norm(derr_db);
                
                bs(n,:) = b;
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
        plts(1).nextval{1}.YData = err; plts(1).nextval{2}.YData = derr; 
        %    plts(1).nextval{3}.YData = gradnorm;
        plts(2).nextval{1}.YData = mu;
        for j = 1:length(b_new)
            plts(3).nextval{j}.YData = b_new(j);
        end
        
        % update plots
        if ~retry
            plts(1).plt{1}.YData = errs; plts(1).plt{2}.YData = derrs; 
            %    plts(1).plt{3}.YData = gradnorms;
            plts(2).plt{1}.YData = mus;
            for j = 1:length(b_new)
                plts(3).plt{j}.YData = bs(:,j);
            end
        end
        
        % update title, view range
        for hp = plts(1:2)
            hp.t.String = num2str(hp.nextval{1}.YData(1));
            %{
            win = [max(1, n-50), max(n, 10)];
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
            %}
        end
        
        pause(.00001);
    
    end
    
end

beta = b';

end