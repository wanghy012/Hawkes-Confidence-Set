function [t,u,Udst] = run(A,mu,T)
Mu = sum(mu);
R = size(A,1);
t = 0;
u = 0;
dst = zeros(R,1);
Udst = zeros(R,1);

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
    Udst(:,end) = Udst(:,end)*exp(-gap);
    Udst(:,end+1) = Udst(:,end);
    Udst(u(end),end) = Udst(u(end),end)+1;
end
t=t(2:end);
u=u(2:end);
Udst = Udst(:,1:end-1);
end
