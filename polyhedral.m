function [lb,ub] = polyhedral(A_MLE,mu,t,u,T,eta,epsilon)
%Compute the confidence set C_{i,epsilon}^p. Requires A_MLE be
%non-negative.
%implemented for phi_ij = exp(-t)
%Input: A_MLE, the MLE
%mu, the true parameter
%t, an array of event time
%u, index of nodes where events happen
%T, total time
%eta, each column is eta at time t(i). Since we assume a universal phi, eta(t) is the same for each node
%epsilon, overall confidence level. example: epsilon = 0.05.
%Output: the upper and lower bound for each entry alpha_ij
D = length(mu);
lb = zeros(D,D);
ub = zeros(D,D);
epsilon = epsilon/2/D/D;
for i=1:D % for node i
    g_alpha = zeros(2*D,D); %g'(alpha_i hat)
    g_alpha2 = zeros(2*D,D);
    g_alpha3 = zeros(2*D,D);
    g_MLE = zeros(2*D,1);%g(alpha_i hat)
    I_hat = zeros(D,D); %t times the estimated Fisher Information at MLE.
    z = zeros(D,2*D);
    I_alpha = zeros(D,D,D);%symmetric. 3d tensor.
    V = zeros(2*D,1);
    V_cum = zeros(2*D,length(t)-1);
    V_2_approx = zeros(2*D,1);
    %update the values above each time an event happens
    for k=1:length(t)-1
        % compute z for [t(k),t(k+1)).
        r = rank(I_hat);
        if r<D
            z = kron([1,-1],eye(D)*sqrt(2*log(1/epsilon)))/sqrt(T);
            z_alpha = zeros(D,2*D,D);
        else
            inverse = inv(I_hat);
            z = kron([1,-1],inverse.*sqrt(2*t(k)*log(1/epsilon)./diag(inverse)'/T));
            for j = 1:D % get the partial derivative of z over each entry alpha_ij
                inverse_alpha = -inverse*I_alpha(:,:,j)*inverse;
                z_alpha(:,:,j) = kron([1,-1],...
                inverse_alpha.*sqrt(2*t(k)*log(1/epsilon)./diag(inverse)'/T)...
                -sqrt(2*t(k)*log(1/epsilon)/T)*inverse.*diag(inverse_alpha)'./(diag(inverse)'.^1.5)/2....
                );
            end
        end
        %eta2: eta at t(k)
        eta2 = eta(:,k);
        eta2(u(k)) = eta2(u(k))+1;
        %update g_MLE till t(k+1)
        g_MLE = g_MLE - z'*eta2*(1-exp(t(k)-t(k+1)));
        if u(k+1) == i
            g_MLE = g_MLE + z'*eta(:,k+1)/(A_MLE(i,:)*eta(:,k+1)+mu(i));
            
        end
        f = @(x) (mu(i) + A_MLE(i,:)*eta2*exp(-x))...
            *exp(z'*eta2*exp(-x)/(mu(i) + A_MLE(i,:)*eta2*exp(-x)))...
            -z'*eta2*exp(-x) - (mu(i) + A_MLE(i,:)*eta2*exp(-x));%integrand for V
        temp = integral(f,0,t(k+1)-t(k),'ArrayValued',1);
        g_MLE = g_MLE - temp;
        V = V + temp;

        %%update g_alpha till t(k+1)
        %%partial dS/partial alpha * z
        if u(k+1) == i
            g_alpha = g_alpha ...
                - z'*eta(:,k+1)*eta(:,k+1)'/(A_MLE(i,:)*eta(:,k+1)+mu(i))^2;
        end
        %%dS*partial z/partial alpha
        g_alpha = g_alpha - reshape(sum(z_alpha.*eta2,1),2*D,D)*(1-exp(t(k)-t(k+1)));
        if u(k+1) == i
            g_alpha = g_alpha + reshape(sum(z_alpha.*eta(:,k+1),1),2*D,D)/(mu(i)+A_MLE(i,:)*eta(:,k+1));
        end
        
        
        
        %%partial V/partial alpha
        f = @(x) eta2'*exp(-x).*exp(z'*eta2*exp(-x)/(mu(i) + A_MLE(i,:)*eta2*exp(-x)))...
            +(mu(i) + A_MLE(i,:)*eta2*exp(-x))...
            *exp(z'*eta2*exp(-x)/(mu(i) + A_MLE(i,:)*eta2*exp(-x)))...
            .*(reshape(sum(z_alpha.*eta2,1),2*D,D)*exp(-x)/(mu(i) + A_MLE(i,:)*eta2*exp(-x))...
            -z'*eta2*exp(-x).*eta2'*exp(-x)/(mu(i) + A_MLE(i,:)*eta2*exp(-x))^2)...
            -reshape(sum(z_alpha.*eta2,1),2*D,D)*exp(-x) - eta2'*exp(-x);
        temp = integral(f,0,t(k+1)-t(k),'ArrayValued',1);
        g_alpha = g_alpha - temp;
        g_alpha3  = g_alpha3 -  temp;
        f = @(x) eta2'*exp(-x).*exp(z'*eta2*exp(-x)/(mu(i) + A_MLE(i,:)*eta2*exp(-x)))...
            -exp(z'*eta2*exp(-x)/(mu(i) + A_MLE(i,:)*eta2*exp(-x)))...
            .*z'*eta2*exp(-x).*eta2'*exp(-x)/(mu(i) + A_MLE(i,:)*eta2*exp(-x))...
            - eta2'*exp(-x);
        g_alpha2 = g_alpha2 - integral(f,0,t(k+1)-t(k),'ArrayValued',1);
        %%update I_hat and I_alpha
        f = @(x) exp(-2*x)/(mu(i) + (A_MLE(i,:)*eta2)*exp(-x));
        I_hat = I_hat + eta2*eta2'*integral(f,0,t(k+1)-t(k),'ArrayValued',1);
       
        f = @(x) exp(-3*x)/(mu(i) + A_MLE(i,:)*eta2*exp(-x))^2;
        I_alpha = I_alpha ...
            - reshape(kron(eta2*eta2',eta2),D,D,D)*integral(f,0,t(k+1)-t(k),'ArrayValued',1);
    end
    
    for j=1:D %get confidence interval for each alpha_ij by linear programming.
        ej = zeros(D,1);
        ej(j) = 1;
        options = optimoptions('linprog','Display','none');
        [~,lb(i,j)]  = linprog(ej,g_alpha,...
            log(1/epsilon) - g_MLE,[],[],[],[],options);

        [~,ub(i,j)]  = linprog(-ej,g_alpha,...
            log(1/epsilon)-g_MLE,[],[],[],[],options);
        clear ej;
    end
end
ub = A_MLE-ub;
lb = A_MLE+lb;
end
    
    
