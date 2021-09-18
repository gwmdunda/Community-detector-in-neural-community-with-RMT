diffs = diff(ts);     % to returns
p = size(diffs,2);      % N is the number of time series
n = size(diffs,1);      % T is the lenght of each series

C = corr(diffs);    % Create a covariance matrix and ensure
C = .5 * (C+C');

eigv = eig(C);
figure(1);
histn(eigv,0,0.2,90);
xlabel('eigenvalue');
ylabel('PDF');
hold on;

variance = 1 - max(eigv)/p;
Q = n/p;
d_pos = (1+1/Q+2*sqrt(1/Q))*variance;
d_neg = (1+1/Q-2*sqrt(1/Q))*variance;
d = d_neg:0.01:d_pos;
probs = Q/(2*pi)*sqrt((d_pos-d).*(d-d_neg))./(d*variance);
plot(d, probs,'r','LineWidth',1);