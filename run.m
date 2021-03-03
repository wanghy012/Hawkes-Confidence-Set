function [t,u,eta] = run(A,mu,T)
Mu = sum(mu);
R = size(A,1);
t = 0;
u = 0;
dst = zeros(R,1);
eta = zeros(R,1);

while 1
    gap = get_next(Mu,sum(dst));
    if t(end) + gap > T
        break;
    end
    t(end+1) = t(end) + gap;
    dst = dst * exp(-gap);
    r = rand()*(Mu+sum(dst));
    v=1;
    s=0;
    for j=1:R
        s = s+mu(j)+dst(j);
        if r>s
            v = v+1;
        else
            break;
        end
    end
    u(end+1) = v;
    dst = dst + A(:,u(end));
    eta(:,end) = eta(:,end)*exp(-gap);
    eta(:,end+1) = eta(:,end);
    eta(u(end),end) = eta(u(end),end)+1;
end
t=t(2:end);
u=u(2:end);
eta = eta(:,1:end-1);
end
