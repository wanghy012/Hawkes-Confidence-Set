function [ub_a,lb_a,ub_p,lb_p,A_MLE,h] = CI_comparison_over_time(t,u,eta,T,epsilon,mu,A) 
D = length(mu);
% decide the ts for CI computation 
K = 10;
t_tic = exp((-K+0.5:-0.5)/K*2)*T;
ub_a = zeros(length(t_tic),D,D);
lb_a = zeros(length(t_tic),D,D);
ub_p = zeros(length(t_tic),D,D);
lb_p = zeros(length(t_tic),D,D);
A_MLE = zeros(length(t_tic),D,D);
% for each t_tic, compute the CI
for k=1:length(t_tic)
    fprintf('iteration = %d\n',k);
    idx = max(find(t<t_tic(k)));
    [A_MLE(k,:,:),I] = MLE_nonnegative(t(1:idx),u(1:idx),t_tic(k),eta(:,1:idx),mu,A);
    
    for i=1:D
        Iinv = inv(I{i});
        ub_a(k,i,:) = reshape(A_MLE(k,i,:),1,D) - norminv(epsilon/2)*sqrt(diag(Iinv))';
        lb_a(k,i,:) = reshape(A_MLE(k,i,:),1,D) + norminv(epsilon/2)*sqrt(diag(Iinv))';
    end
    [lb_p(k,:,:),ub_p(k,:,:)] = polyhedral22(max(0,reshape(A_MLE(k,:,:),D,D)),mu,t(1:idx),u(1:idx),t_tic(k),eta(:,1:idx),epsilon*D*D);
end 
%plot CI over time.
h = figure(1);
for i=1:D
    for j=1:D
        subplot(D,D,D*(i-1)+j);
        hold on;
        plot(1:length(t_tic),ub_a(:,i,j),'r-');
        plot(1:length(t_tic),lb_a(:,i,j),'r-');
        plot(1:length(t_tic),ub_p(:,i,j),'b-');
        plot(1:length(t_tic),lb_p(:,i,j),'b-');
        plot(1:length(t_tic),A(i,j)*ones(length(t_tic),1),'k-');
        plot(1:length(t_tic),A_MLE(:,i,j),'g-');
%         if A(i,j)>0
%             plot(xaxis,zeros(length(xaxis),1),'k--');
%         end
        axis([1,length(t_tic),0,0.6])
    end
end

end