%generate model. influence functions are exp(-t).
% D = 5;
% A = rand(D,D).*(rand(D,D)<2/D)*0.2;
% T = 2000;
% mu = rand(D,1);
% figure(1)
% plot(digraph(A),'EdgeLabel',A(find(A>0)),'LineWidth',2,...
%     'EdgeColor','k','NodeColor','k','ArrowSize',10,'MarkerSize',10,...
%     'EdgeFontSize',20,'NodeFontSize',20);
% %generate data
[t,u,eta] = run(A,mu,T);


[A_MLE,I] = MLE_nonnegative(t,u,T,eta,mu,A);
epsilon = 0.05*2*D*D;
% asymptotic bound
lb_a = zeros(D,D);
ub_a = zeros(D,D);
for i=1:D
    lb_a(i,:) = A_MLE(i,:) + norminv(epsilon/2/D/D)*sqrt(diag(inv(I{i})))';
    ub_a(i,:) = A_MLE(i,:) - norminv(epsilon/2/D/D)*sqrt(diag(inv(I{i})))';
end

%non-asymptotic bound
[lb_p,ub_p,g_MLE,g_alpha] = polyhedral2(A_MLE,mu,t,u,T,eta,epsilon);
