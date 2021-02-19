function [A_MLE,I] = MLE_A(t,u,T,eta,mu,A)
% MLE of A. Use true value as initialization. No constraint.
    d = size(eta,1);
    N = size(eta,2);
    lin = zeros(d,1);
    I = cell(1,d);
    lower_bound = 0;
    for i=1:N
        lin(u(i)) = lin(u(i)) + (1-exp(-T+t(i)));
    end
    A_MLE = zeros(d,d);
    for i=1:d
        alpha = A(i,:);
        t = zeros(1,d);
        H = (eta'./([eta',ones(N,1)]*[alpha,mu(i)]'))'...
            *diag(u==i)*(eta'./([eta',ones(N,1)]*[alpha,mu(i)]'));
        if rank(H,0.1)<d
            a(i,:) = NaN;
            I{i} = NaN*ones(d,d);
            continue;
        end
        while max(abs(t-alpha))>0.0001
            %newton
            t = alpha;
            gd = (u==i)*(eta'./([eta',ones(N,1)]*[t,mu(i)]')) -lin';
            H = (eta'./([eta',ones(N,1)]*[alpha,mu(i)]'))'...
            *diag(u==i)*(eta'./([eta',ones(N,1)]*[alpha,mu(i)]'));
            idx = 1:d;
            while ~prod(alpha(idx)>lower_bound | gd(idx)/H(idx,idx)>=0)
                idx = idx(alpha(idx)>lower_bound | gd(idx)/H(idx,idx)>=0);
            end
            step = gd(idx)/H(idx,idx);
            step = step/max([-step./(alpha(idx)-lower_bound),1]);
            alpha(idx) = t(idx) + step;
        end
        A_MLE(i,:) = alpha;
        I{i} = (eta'./([eta',ones(N,1)]*[alpha,mu(i)]'))'...
            *diag(u==i)*(eta'./([eta',ones(N,1)]*[alpha,mu(i)]'));
    end
end

