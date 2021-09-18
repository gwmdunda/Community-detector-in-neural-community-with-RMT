%Suppose we have matrix ts

diffs = diff(ts);     % to returns
N = size(diffs,2);      % N is the number of time series
T = size(diffs,1);      % T is the lenght of each series

C = cov(diffs);    % Create a covariance matrix and ensure
C = .5 * (C+C');        % it's symmetric

[V,D] = eig(C);
[s_eigvals, ind]=sort(diag(D),'ascend');
V = V(:,ind);
D=diag(sort(diag(D),'ascend'));

syms x
LHS = 0;
for m=1:N
    LHS = LHS + s_eigvals(m)/((s_eigvals(m)-x)*N);
end
S = vpasolve(LHS == T/N,x,'Random',true);
mu = double(S);
est_eigs = T*(s_eigvals - mu);
est_eigs = sort(est_eigs,'ascend');


est_cov = zeros(N,N);
for m=1:N
    est_eigv = zeros(N,N);
    for k=1:N
        if m ~= k
            est_eigv = est_eigv + (-1)*(s_eigvals(m)/(s_eigvals(k)-s_eigvals(m))-mu(m)/(s_eigvals(k)-mu(m)))*V(:,k)*V(:,k).';
        else
            coeff = 0;
            for r=1:N
                if r == m
                    continue
                end
                coeff = coeff + s_eigvals(r)/(s_eigvals(k)-s_eigvals(r))-mu(r)/(s_eigvals(k)-mu(r));
            end
            est_eigv = est_eigv + est_eigs(m)*(1+coeff)*V(:,k)*V(:,k).';
        end
    end
    est_cov = est_cov + est_eigs(m)*est_eigv; 
end
est_cov = .5*(est_cov + est_cov.');
est_corr = corrcov(est_cov);

eigv = sort(eig(est_corr),'ascend');

figure(1) %plot with global mode
histn(eigv,0,0.2,20);
xlabel('eigenvalue');
ylabel('PDF');
hold on;

Q = T/N;
variance = 1 - max(eig(est_corr))/N;
d_pos = variance*(1+1/Q+2*sqrt(1/Q));
d_neg = variance*(1+1/Q-2*sqrt(1/Q));
d = d_neg:0.01:d_pos;
probs = Q/(2*pi)*sqrt((d_pos-d).*(d-d_neg))./(d*variance);
plot(d, probs,'r','LineWidth',1);

figure(2) %plot without global mode
histn(eigv(1:end-1),0,0.2,20);
xlabel('eigenvalue');
ylabel('PDF');
hold on;

Q = T/N;
variance = 0.6;
d_pos = variance*(1+1/Q+2*sqrt(1/Q));
d_neg = variance*(1+1/Q-2*sqrt(1/Q));
d = d_neg:0.01:d_pos;
probs = Q/(2*pi)*sqrt((d_pos-d).*(d-d_neg))./(d*variance);
plot(d, probs,'r','LineWidth',1);
 
