function [A_MLE,I] = MLE_nonnegative(t,u,T,eta,mu,A)
%estimate A with constraint (A > = 0) (use A as initialization)
d = size(eta,1);
N = size(eta,2);
lin = zeros(d,1);
I = cell(1,d);
for i=1:N
    lin(u(i)) = lin(u(i)) + (1-exp(-T+t(i)));
end
A_MLE = zeros(d,d);
for i=1:d
    alpha = ones(1,d);
    t = A(i,:);
    %check if Hessian is rank deficient
    H = (eta'./([eta',ones(N,1)]*[alpha,mu(i)]'))'...
        *diag(u==i)*(eta'./([eta',ones(N,1)]*[alpha,mu(i)]'));
%     if rank(H,0.1)<d
%         A_MLE(i,:) = NaN;
%         I{i} = NaN*ones(d,d);
%         continue;
%     end
    %if Hessian is good
    while max(abs(t-alpha))>0.0001
        %coordinate descent
        alpha=t;
        for j=1:d
            t1 = 2;
            while(abs(t1-t(j))>0.0001)
                t(j)=t1;
                gd = (u==i)*(eta(j,:)'./([eta',ones(N,1)]*[t,mu(i)]')) -lin(j);
                H = (u==i)*((eta(j,:)'./([eta',ones(N,1)]*[t,mu(i)]')).^2);
                if H == 0
                    t1 = 0;
                else
                    t1 = max(0,t1 + gd/H);
                end
            end
        end
    end
    A_MLE(i,:) = alpha;
    I{i} = (eta'./([eta',ones(N,1)]*[alpha,mu(i)]'))'...
        *diag(u==i)*(eta'./([eta',ones(N,1)]*[alpha,mu(i)]'));
end
end


