function [ sigma, Q, M] = FuncSignature(TS, nullModel, cdAlg, estimator )
%% FuncSignature
% source [1] with modification
% FuncSignature takes an array of timeseries and returns a set of
% communities (based on modularity optimization) of the time series and a
% quantitative measure of the modularity.
%% Syntax
% [sigma, Q,M]=FuncSignature(TS,nullModel, cdAlg)
%
%% Description
% This function takes a set of time series (see FinRMT for more details on
% the composition of time series) FuncSignatureand a numerical index which represents a
% specific community detection algorithm to use. It will generate a null
% model using Random Matrix Theory to filter the correlation matrix of the
% timeseries data, by removing the market mode and noise. It will then use
% the chosed community detection algorithm in conjunction with the filtered
% correlation matrix to ascertain a set of communities which are internally
% correlated and externally anti correlated. The results of this
% calcualtion will differ slightly depending on the choice of community
% detection algorithm used and even between different runs using the same
% algorithm.
%
%
%% Inputs arguments:
% TS : an mxn matrix containing timeseries' of stock prices. Each column
% should be a time series for one financial instrument and each row should
% correspond to the value of each instrument at a point in time. For example
% 32.00   9.43   127.25   ...
% 32.07   9.48   126.98   ...
% 32.08   9.53   126.99   ...
%  ...    ...     ....    ...
% No header columns or timestamp columns should be included
%
% nullModel is an integer value specifying the type of null model to use
% when constructing the modularity matrix
% 1 = infinite time series with no global mode
% 2 = finite time series with no global mode
% 3 = finite time series with a global mode
% for more details see the corresponding paper
%
% cdAlg : is an integer value (1, 2 or 3) which representes the type of
% community detection algorithm which will be used
% 1 for the generalized Louvain method
% 2 for Newman's spectral decomposition
% 3 for a Simulated Annealing approach
%
% estimator: is an integer value (1, 2 or 3) which represents preprocessing
% of the estimated eigenvalues
% 1 for no preprocessing
% 2 for preprocessing according to Meistre paper
% 
%
%% Outputs:
% sigma: is a vector of community membership where each integer in the
% vector specifies the community to which the corresponding time
% series from the input matrix is a member. For example sigma = [1, 2, 2,
% 1, 2, 1] implies that timeseries 1,4 & 6 are all in community 1 and time
% series 2,3 & 5 are in community 2.
% M : the filered modulatiry matrix
%
%% Example:
% [communities, Q, Modu]=FuncSignature(myData,3,1,1)
%
%% Disclaimer
% This software is provided as is and by using it the user acknowledges
% that the authors accept no liability whatsoever for the results of it's
% use
%
%% Acknowledgements
% The community detection code provide here are provided as examples of
% community detection algorithms which can be used in conjunction with the
% RMT null model formation. The code for the generalized louvain method is
% provide unalterd from it's original form, developed by :
%        Inderjit S. Jutla, Lucas G. S. Jeub, and Peter J. Mucha,
%        "A generalized Louvain method for community detection implemented
%        in MATLAB," http://netwiki.amath.unc.edu/GenLouvain (2011-2012).
%
% The Newman method has been slightly modified from it's original version
% developed by :
%        Mika Rubinov, UNSW
%        Jonathan Power, WUSTL
%        Dani Bassett, UCSB
% The modifications are to accomodate for the fact that the original
% version depends on certain properties of the modularity matrix that do
% not exist when dealing with correlation matrices. It has also been
% modified to take a modularity matrix directly, instead of constructing
% one from an adjacency matrix.
%
% The Simmulated Annealing method was written by us as a proof of concept
% for the use of Simulated Annealing. It is fairly rudementery in it's
% applicatoin of Simulated Annealing and does not expose many of the
% traditional parameters of the algorithm
%%


% Construct the filtered correlation matrix by using RMT to define the null
% model and remove it from the original correlation matrix
M = FinRMT(TS,nullModel,estimator);

% Using the filterd correlation matrix as the modularity matrix, we can now
% use conventional, network based community detection algorithms to
% determine the community structure. Note that some minor modifictions were
% necessary to the Newman method (details are provided in the corresponding
% paper)

if (cdAlg == 1)
    sigma = genlouvain(M,2,0,0);
    sigma = sigma';
elseif (cdAlg == 2)
    sigma = Newman(M);
    sigma = sigma';
elseif (cdAlg == 3)
    sigma = Potts(M);
end


% Define a normalization for the value of Q, depending on the type of time
% series being used a different type of normalization might be appropriate.
% So QNorm can be redefined as necessary.
% e.g. QNorm = sum(abs(M(:)));
QNorm = sum(M(:));

% Calculate the Modularity
sMat = repmat(sigma,size(sigma,2),1);
mod=~(sMat-sMat').*M;
mod=sum(mod(:));
Q = mod/QNorm;

end


function M = FinRMT(TS, nullModel, estimator)
%% FinRMT
% FinRMT uses Random Matrix Theory (RMT) to create a filtered correlation
% matrix from a set of financial time series price data, for example the
% daily closing prices of the stocks in the S&P
%% Syntax
% M=FinRMT(TS, nullModel, estimator)
%
%% Description
% This function eigendecomposes a correlation matrix of time series
% and splits it into three components, Crandom, Cgroup and Cmarket,
% according to techniques from the literature (See, "Systematic Identification
% of Group Identification in Stock Markets, Kim & Jeong, (2008).") and
% returns a filtered correlation matrix containging only the Cgroup
% components.
% The function is intended to be used in conjunction with a community
% detection algorithm (such as the Louvain method) to allow for community
% detecion on time series based networks.
%
%
%% Inputs arguments:
% TS : an mxn matrix containing timeseries of signals. Each column
% should be a time series for one ROI and each row should
% correspond to the value of each ROI at a point in time. For example
% 32.00   9.43   127.25   ...
% 32.07   9.48   126.98   ...
% 32.08   9.53   126.99   ...
%  ...    ...     ....    ...
% No header columns or timestamp columns should be included
%
% nullModel is an integer value specifying the type of null model to use
% when constructing the modularity matrix
% 1 = infinite time series with no global mode
% 2 = finite time series with no global mode
% 3 = finite time series with a global mode
% for more details see the corresponding paper
%% Outputs:
% M : The filtered correlation matrix. This matrix can be passed directly to
% a community detection algorithm in place of the modularity matrix
%
%% Example:
% ModularityMatrix = FinRMT(myData, 3)
%  ...
% Communities = myCommunityDectionAlg(ModularityMatrix)
%
%% Issues & Comments
% Note that the output of this function can serve as the Modularity
% Matrix (Not the Adjacency matrix) for a generalized Community Detection
% Algorithm. Specifically, one which does not rely on properties of the
% Adjaceny Matrix to create the Modularity Matrix. The Louvain Method
% and methods based on simulated annealing are examples of such.
%
%%

    
%diffs = diff(log(TS));     % to log returns
diffs = diff(TS);     % to returns
N = size(diffs,2);      % N is the number of time series
T = size(diffs,1);      % T is the lenght of each series


%Create the initial correlation matrix and ensure it's symmetric
%It should be symmetric but sometimes matlab introduces small roundoff
%errors that would prevent an IsSymmetric call from returning true.
%This folding of the matrix will suffice to solve that.

C = cov(diffs);    % Create a correlation matrix and ensure
C = .5 * (C+C');        % it's symmetric

%preprocess correlation matrix
if estimator == 1
    C = corrcov(C);
elseif estimator == 2
    [V,D] = eig(C);
    [s_eigvals, ind]=sort(diag(D),'ascend');
    V = V(:,ind);
    
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
    C = est_corr;
end
    

% Decompose the correlation matrix into its eigenvalues and eigenvectors,
% store the indices of which columns the sorted eigenvalues come from
% and arrange the columns in this order

% if the user request the null model for infinite time series then we
% simply subtract the identity from the correlation matrix and return the
% result. Otherwise we proceed to use RMT to filter out the noise expected
% from a non infinite time series.
if nullModel == 1
    M = C-eye(N);
else
    [V,D] = eig(C);
    [eigvals, ind]=sort(diag(D),'ascend');
    V = V(:,ind);
    D=diag(sort(diag(D),'ascend'));
    
    
    % Find the index of the predicted lambda_max, ensuring to check boundary
    % conditions
    Q=T/N;
    sigma =  1 - max(eigvals)/N;
    RMTmaxEig = sigma*(1 + (1.0/Q) + 2*sqrt(1/Q));
    RMTmaxIndex = find(eigvals > RMTmaxEig);
    if isempty(RMTmaxIndex)
        RMTmaxIndex = N;
    else
        RMTmaxIndex = RMTmaxIndex(1);
    end
    
    % Find the index of the predicted lambda_min, ensuring the check boundary
    % conditions
    RMTminEig = sigma*(1 + (1.0/Q) - 2*sqrt(1/Q));
    RMTminIndex = find(eigvals < RMTminEig);
    if isempty(RMTminIndex)
        RMTminIndex = 1;
    else
        RMTminIndex = RMTminIndex(end);
    end
    
    
    % Determine the average Eigenvalue to rebalance the matrix after removing
    % Any of the noise and/or market mode components
    avgEigenValue = mean(eigvals(1:(RMTmaxIndex-1)));
    
    
    
    % Build a new diagonal matrix consisting of the group eigenvalues
    Dg = zeros(N,N);
    
    % Replace the random component with average values.
    Dg(1 : (N+1) : (RMTmaxIndex-1)*(N+1)) = avgEigenValue;
    
    % Add the group component. The N+1 here is just used to increment to the
    % next diagonal element in the matrix
    Dg(1+(N+1)*(RMTmaxIndex-1) : (N+1) : end-(N+1)) = D(1+(N+1)*(RMTmaxIndex-1) : (N+1) : end-(N+1));
    
    
    % If the user requested a nullModel that includes the market mode component
    % that we add it in here. If instead they wanted to remove the market mode
    % then we do nothing, as it's already missing from the new matrix we're
    % constructing, Dg.
    if nullModel == 2
        Dg(N,N) = D(N,N);
    end
    % Build the component correlation matrix from the new diagonal eigenvalue
    % matrix and eigenvector matrix. The eigenvectors corresponding to zero
    % valued eigenvalue entries in Dg will not contribute to M
    
    M = V * Dg * V.';
    % Replace the diagonals with 1s
    M = M - diag(diag(M)) + eye(N);
    
end

end

function [S,Q] = genlouvain(B,limit,verbose,randord)
%GENLOUVAIN  Louvain-like community detection, specified quality function.
%   Version 1.2 (July 2012)
%
%   [S,Q] = GENLOUVAIN(B) with matrix B implements a Louvain-like greedy
%   community detection method using the modularity/quality matrix B that
%   encodes the quality function Q, defined by summing over all elements
%   B(i,j) such that nodes i and j are placed in the same community.
%   Following Blondel et al. 2008, the algorithm proceeds in two phases
%   repeated iteratively: quality is optimized by moving one node at a time
%   until no such moves improve quality; the communities found to that
%   point are then aggregated to build a new network where each node
%   represents a community.  The output vector S encodes the obtained
%   community assignments, with S(i) identifying the community to which
%   node i has been assigned.  The output Q gives the quality of the
%   resulting partition of the network.
%
%   [S,Q] = GENLOUVAIN(B) with function handle B such that B(i) returns
%   the ith column of the modularity/quality matrix uses this function
%   handle (to reduce the memory footprint for large networks) until the
%   number of groups is less than 10000 and then builds the B matrix
%   corresponding to the new aggregated network in subsequent passes.  Use
%   [S,Q] = GENLOUVAIN(B,limit) to change this default=10000 limit.
%
%   [S,Q] = GENLOUVAIN(B,limit,0) suppresses displayed text output.
%
%   [S,Q] = GENLOUVAIN(B,limit,verbose,0) forces index-ordered (cf.
%   randperm-ordered) consideration of nodes, for deterministic results.
%
%   Example (using adjacency matrix A)
%         k = full(sum(A));
%         twom = sum(k);
%         B = @(v) A(:,v) - k'*k(v)/twom;
%         [S,Q] = genlouvain(B);
%         Q = Q/twom;
%     finds community assignments for the undirected network encoded by the
%     symmetric adjacency matrix A.  For small networks, one may obtain
%     reasonably efficient results even more simply by handling the full
%     modularity/quality matrix
%         B = A - k'*k/twom;
%     instead of the function handle.  Intended use also includes the
%     "multislice" network quality function of Mucha et al. 2010, where B
%     encodes the interactions as an equivalent matrix (see examples posted
%     online at http://netwiki.amath.unc.edu/GenLouvain).
%
%   Notes:
%     The matrix represented by B must be both symmetric and square.  This
%     condition is not checked thoroughly if B is a function handle, but is
%     essential to the proper use of this routine.
%
%     Under default options, this routine can return different results from
%     run to run because it considers nodes in pseudorandom (randperm)
%     order.  Because of the potentially large number of nearly-optimal
%     partitions (Good et al. 2010), one is encouraged to investigate
%     results of repeated applications of this code (and, if possible, of
%     other computational heuristics).  To force deterministic behavior,
%     ordering nodes by their index, pass zero as the fourth input:
%     GENLOUVAIN(B,limit,verbose,0).
%
%     This algorithm is only "Louvain-like" in the sense that the two
%     phases are used iteratively in the same manner as in the Louvain
%     algorithm (Blondel et al. 2008).  Because it operates on a general
%     quality/modularity matrix B, it does not include any analytical
%     formulas for quickly identifying the change in modularity from a
%     proposed move nor any improved efficiency obtained by their use.  If
%     your problem uses one of the well-used null models included in other
%     codes, those codes should be much faster for your task.
%
%     Past versions had a problem where accumulated subtraction error might
%     lead to an infinite loop with each pass oscillating between two or
%     more partitions yet incorrectly identifying increases in quality.  We
%     believe this problem has been corrected by the relative change checks
%     in lines 178 and 273.  If you encounter a similar problem, notify
%     Peter Mucha (<a href="mailto:mucha@unc.edu">mucha@unc.edu</a>).
%
%     The output Q provides the sum over the appropriate elements of B
%     without any rescaling.  As such, we have rescaled Q in the example
%     above by 2m = sum(k) so that Q <= 1.
%
%     The '~' for ignoring function returns (used for "max" below) are not
%     supported prior to R2009b.  Replace (e.g. 'dummy') for pre-2009b.
%
%     By using this code, the user implicitly acknowledges that the authors
%     accept no liability associated with that use.  (What are you doing
%     with it anyway that might cause there to be a potential liability?!?)
%
%   References:
%     Blondel, Vincent D., Jean-Loup Guillaume, Renaud Lambiotte, and
%     Etienne Lefebvre, "Fast unfolding of communities in large networks,"
%     Journal of Statistical Mechanics: Theory and Experiment, P10008
%     (2008).
%
%     Fortunato, Santo, "Community detection in graphs," Physics Reports
%     486, 75-174 (2010).
%
%     Mucha, Peter J., Thomas Richardson, Kevin Macon, Mason A. Porter, and
%     Jukka-Pekka Onnela. "Community Structure in Time-Dependent,
%     Multiscale, and Multiplex Networks," Science 328, 876-878 (2010).
%
%     Porter, M. A., J. P. Onnela, and P. J. Mucha, "Communities in
%     networks," Notices of the American Mathematical Society 56, 1082-1097
%     & 1164-1166 (2009).
%
%   Acknowledgments:
%     A special thank you to Stephen Reid, whose greedy.m code was the
%     original version that has over time developed into the present code.
%     Thank you also to Dani Bassett, Jesse Blocher, Mason Porter and Simi
%     Wang for inspiring improvements to the code.
%
%   Citation: If you use this code, please cite as
%        Inderjit S. Jutla, Lucas G. S. Jeub, and Peter J. Mucha,
%        "A generalized Louvain method for community detection implemented
%        in MATLAB," http://netwiki.amath.unc.edu/GenLouvain (2011-2012).

%set default for maximum size of modularity matrix
if nargin<2
    limit = 10000;
end

%set level of reported/displayed text output
if nargin<3
    verbose = 1;
end
if verbose
    mydisp = @(s) disp(s);
else
    mydisp = @(s) disp('');
end

%set randperm- v. index-ordered
if nargin<4
    randord = 1;
end
if randord
    myord = @(n) randperm(n);
else
    myord = @(n) 1:n;
end

%initialise variables and do symmetry check
if isa(B,'function_handle')
    n=length(B(1));
    S=(1:n)';
    M=B;
    it(:,1)=M(1);
    ii=find(it(2:end)>0,3)+1;
    ii=[1,ii'];
    for i=2:length(ii),
        it(:,i)=M(ii(i));
    end
    it=it(ii,:);
    if nnz(it-it'),
        disp('WARNING: Function handle does not correspond to a symmetric matrix')
    end
else
    n = length(B);
    S = (1:n)';
    M=B;
    if nnz(M-M'),
        B=(B+B')/2; disp('WARNING: Forced symmetric B matrix')
    end
end

dtot=0; %keeps track of total change in modularity

%Run using function handle, if provided
while (isa(M,'function_handle')) %loop around each "pass" (in language of Blondel et al) with B function handle
    
    y = unique(S); %unique also puts elements in ascending order
    Sb=S;
    yb = [];
    
    clocktime=clock;
    mydisp(['Merging ',num2str(length(y)),' communities  ',num2str(clocktime(4:6))]);
    
    dstep=1;	%keeps track of change in modularity in pass
    
    while (~isequal(yb,y))&&(dstep/dtot>2*eps) %This is the loop around Blondel et al's "first phase"
        %        Q = 0;
        %        %improves performance considerably if one doesn't compute modularity
        %        %for the first pass (for display purposes only)
        %         P = sparse(y,1:length(y),1); %Modularity Calculation
        %         for i = 1:length(M(1))
        %             Q = Q + (P*M(i))'*P(:,i);
        %         end
        %         mydisp([num2str(length(unique(y))),' ',num2str(Q)])
        yb = y;
        G=sparse(1:length(y),y,1);      %no-mex version
        dstep=0;
        
        for i = myord(length(M(1)))     %loop over nodes in pseudorandom order
            
            Mi = M(i);
            
            u = unique([y(i);y(Mi>0)]);
            
            dH=Mi'*G(:,u);              %no-mex version
            %dH=modchange_y(Mi,y,u);
            
            yi=find(u==y(i));
            dH(yi) = dH(yi) - Mi(i);
            
            [~, k] = max(dH);
            
            %only move to different group if it is more optimized than
            %staying in same group (up to error with double precision)
            if(dH(k)>(dH(yi)))
                dtot=dtot+dH(k)-dH(yi);
                dstep=dstep+dH(k)-dH(yi);
                G(i,y(i))=0;            %no-mex version
                G(i,u(k))=1;            %no-mex version
                y(i) = u(k);
            end
            
        end
        
        mydisp([num2str(length(unique(y))),' change: ',num2str(dstep),...
            ' total: ',num2str(dtot),' relative: ',num2str(dstep/dtot)]);
    end
    
    %[S,y] = tidyconfig_c(S,y);  %note tidyconfig reorders along node numbers
    y = tidyconfig(y);                  %no-mex version
    for i = 1:length(y)                 %no-mex version
        S(S==i) = y(i);                 %no-mex version
    end                                 %no-mex version
    
    %calculate modularity and return if converged
    if isequal(Sb,S)
        Q=0;
        P=sparse(y,1:length(y),1);
        for i=1:length(M(1))
            Q=Q+(P*M(i))'*P(:,i);
        end
        return
    end
    
    %check wether #groups < limit
    t = length(unique(S));
    if (t>limit)
        M=@(i) metanetwork_i(B,S,t,i);   %use function handle if #groups>limit
    else
        J = zeros(t);   %convert to matrix if #groups small enough
        for c=1:t
            J(:,c)=metanetwork_i(B,S,t,c);
        end
        B = J;
        M=B;
    end
    
end


S2 = (1:length(B))';
Sb = [];
while ~isequal(Sb,S2) %loop around each "pass" (in language of Blondel et al) with B matrix
    
    y = unique(S2);  %unique also puts elements in ascending order
    Sb = S2;
    
    clocktime=clock;
    mydisp(['Merging ',num2str(length(y)),' communities  ',num2str(clocktime(4:6))]);
    
    yb = [];
    
    G=sparse(1:length(y),y,1);
    
    dstep=1;
    
    % P = G';
    % Q = sum(sum((P*M).*(P)));
    % Qb = -inf;
    
    while (~isequal(yb,y)) && (dstep/dtot>2*eps) %This is the loop around Blondel et al's "first phase"
        
        % mydisp([num2str(length(unique(y))),' ',num2str(Q)])
        yb = y;
        % Qb=Q;
        
        dstep=0;
        
        for i = myord(length(M))
            u = unique([y(i);y(M(:,i)>0)]);
            % dH = modchange_y(M(:,i),y,u); %relative changes in modularities
            dH = (M(:,i)'*G(:,u));
            
            yi=find(u==y(i));
            dH(yi) = dH(yi) - M(i,i);
            [~, k] = max(dH);
            %%only move to different group if it is more optimized than
            %%staying in same group (up to error with double precision)
            if(dH(k)>(dH(yi)))
                dtot=dtot+dH(k)-dH(yi);
                dstep=dstep+dH(k)-dH(yi);
                G(i,y(i))=0;
                G(i,u(k))=1;
                y(i) = u(k);
            end
        end
        
        % P=sparse(y,1:length(y),1);
        % Q = sum(sum((P*M).*(P)));
        
    end
    
    y = tidyconfig(y);  %note tidyconfig reorders along node numbers
    for i = 1:length(y)
        S(S==i) = y(i);
        S2(S2==i) = y(i);
    end
    
    if isequal(Sb,S2)
        P=G';
        Q=sum(sum((P*M).*P));
        return
    end
    
    M = metanetwork(B,S2);
end
end
%-----%
function M = metanetwork(J,S)
%Computes new aggregated network (communities --> nodes)
if(issparse(J))
    m=max(S);
    [i,j,v]=find(J);
    M = sparse(S(i),S(j),v,m,m);
else
    PP = sparse(1:length(S),S,1);
    M = PP'*J*PP;
end
end
%-----%
function Mi = metanetwork_i(J,S,t,i)
%ith column of metanetwork (used to create function handle)
%J is a function handle

Mi=sparse([],[],[],t,1);
for j=find(S==i)'
    Jj=J(j);
    [ii,k,v]=find(Jj);
    Mi=Mi+sparse(S(ii),k,v,t,1);
end
end
%-----%
function S = tidyconfig(S)
%This function remains almost identical to that originally written by
%Stephen Reid for his greedy.m code.
%   tidy up S i.e.  S = [2 4 2 6] -> S = [1 2 1 3]
T = zeros(length(S),1);
for i = 1:length(S)
    if T(i) == 0
        T(S==S(i)) = max(T) + 1;
    end
end
S = T;
end

function Ci=Newman(C)
%MODULARITY_DIR     Optimal community structure and modularity
%
%   Ci = modularity_dir(W);
%   [Ci Q] = modularity_dir(W);
%
%   The optimal community structure is a subdivision of the network into
%   nonoverlapping groups of nodes in a way that maximizes the number of
%   within-group edges, and minimizes the number of between-group edges.
%   The modularity is a statistic that quantifies the degree to which the
%   network may be subdivided into such clearly delineated groups.
%
%   Input:      C,      Covariance Matrix
%
%   Outputs:    Ci,     optimal community structure
%               Q,      maximized modularity
%
%   Note: Ci and Q may vary from run to run, due to heuristics in the
%   algorithm. Consequently, it may be worth to compare multiple runs.
%   Also see Good et al. (2010) Phys. Rev. E 81:046106.
%
%   Reference: Leicht and Newman (2008) Phys Rev Lett 100:118703.
%
%
%   2008-2010
%   Mika Rubinov, UNSW
%   Jonathan Power, WUSTL
%   Dani Bassett, UCSB


%   Modification History:
%   Jul 2008: Original (Mika Rubinov)
%   Oct 2008: Positive eigenvalues are now insufficient for division (Jonathan Power)
%   Dec 2008: Fine-tuning is now consistent with Newman's description (Jonathan Power)
%   Dec 2008: Fine-tuning is now vectorized (Mika Rubinov)
%   Sep 2010: Node identities are now permuted (Dani Bassett)

N=length(C);                            %number of vertices

n_perm = randperm(N);                   %DB: randomly permute order of nodes
C = C(n_perm,n_perm);                   %DB: use permuted matrix for subsequent analysis
%Ki=sum(A,1);                            %in-degree
%Ko=sum(A,2);                            %out-degree
%m=sum(Ki);                           	%number of edges
%b=A-(Ko*Ki).'/m;
%B=b+b.';                                %directed modularity matrix
Ci=ones(N,1);                           %community indices
cn=1;                                   %number of communities
U=[1 0];                                %array of unexamined communites

ind=1:N;
Cg=C;
Ng=N;
Ctot = sum(Cg(:));
while U(1)                              %examine community U(1)
    [V D]=eig(Cg);
    [d1 i1]=max(diag(D));               %most positive eigenvalue of Bg
    v1=V(:,i1);                         %corresponding eigenvector
    %Determine sum of elements
    S=ones(Ng,1);
    S(v1<0)=-1;                         %set s(i) to -1 when corresponding constiuents in v are negative
    q=S.'*Cg*S + Ctot;                  %contribution to modularity
    
    if q>1e-10                       	%contribution positive: U(1) is divisible
        qmax=q;                          %maximal contribution to modularity
        Cg(logical(eye(Ng)))=0;      	%Cg is modified, to enable fine-tuning
        indg=ones(Ng,1);                %array of unmoved indices
        Sit=S;
        while any(indg);                %iterative fine-tuning
            Qit=qmax-Sit.*(Cg*Sit); 	%this line is equivalent to:
            qmax=max(Qit.*indg);        %for i=1:Ng
            imax=(Qit==qmax);           %	Sit(i)=-Sit(i);
            Sit(imax)=-Sit(imax);       %	Qit(i)=Sit.'*Bg*Sit;
            indg(imax)=nan;             %	Sit(i)=-Sit(i);
            if qmax>q;                  %end
                q=qmax;
                S=Sit;
            end
        end
        
        if abs(sum(S))==Ng              %unsuccessful splitting of U(1)
            U(1)=[];
        else
            cn=cn+1;
            Ci(ind(S==1))=U(1);         %split old U(1) into new U(1) and into cn
            Ci(ind(S==-1))=cn;
            U=[cn U];                   %Grow U and add the next community to the fron
        end                             %U behaves like a stack of communities
    else                                %contribution nonpositive: U(1) is indivisible
        U(1)=[];                        %Pop the most recent community of the stack
    end
    Ctot=0;
    ind=find(Ci==U(1));                 %indices of unexamined community U(1)
    cg=C(ind,ind);
    Cg=cg-diag(sum(cg));                %modularity matrix for U(1)
    Ng=length(ind);                     %number of vertices in U(1)
end


Ci_corrected = zeros(N,1);              % DB: initialize Ci_corrected
Ci_corrected(n_perm) = Ci;              % DB: return order of nodes to the order used at the input stage.
Ci = Ci_corrected;                      % DB: output corrected community assignments
end


function sigma = Potts(M)
% Simulated Annealing algorithm for solving Modularity Maximization based
% on the Potts Model for Community Detection
%
% Input:
% - M                   - The modularity matrix
%
% Output:
% - sigma               - a vector such that node i is in community sigma[i]
%
% Author: Mel MacMahon
% Last modified: Sept 15st, 2012


% eval_budget represents the speed of annealing. The larger this number the
% longer the algorithm will take to run but the more accurate the results
% will be. Programatically, it represents the number of function
% evaluations available to the algorithm.
eval_budget = 1000;

% N is the number of nodes
N = size(M,2);

% Establish initial constructs
% sigma[i] stores the community (index into comNodeCount) to which
% each stock, i belongs to. We star off with each stock in its own
% community.
sigma = (1:1:N);

% comNodeCount[i] is the number of stocks in community i. As the annealing
% process continues, more and more communities will become empty and some
% rescaling will need to occur whereby we will reduce the size of
% comNodeCount to match the active number of communities plus one empty
% community. To begin with each community contains one node.
comNodeCount = ones(1,N);

% acl[i] stores the index of the ith active community (community with nodes
% in it) in comNodesCount. acl will also rescale as the number of active
% communities diminishes (or grows, if that were to ever happen). To begin
% with, there are N active communities.
acl = (1:1:N);

% q is the number of active communities, so it starts off as N and will
% decrease as the number of active communities diminishes.
q = size (acl, 2);

% CN is the number of communites that were active at the time of the last
% rescaling. It starts off as N since we start with N communities and will
% rescale with comNodeCount.
CN = size (comNodeCount, 2);


% K is the number of iterations over which the temperature is kept constant
k = ceil(eval_budget/100);

% Set temperature and cooling conditions
T =100;
alpha = 0.001;


while T > 1
    
    for j = 1:k
        
        % Choose a node at random
        node = randi(N,1);
        
        % Asses the strength of its correlation with its community
        ch = cohesion (node, sigma, M, N);
        
        % Choose a new community at random but such that it is not the
        % community to which node currently belongs
        cNew = sigma(node);
        while (cNew == sigma(node))
            cNew = randi(CN,1);
        end
        
        % Calculate the Hamiltonian of the system where node is moved from
        % its current community to the new one, cNew
        Hnew = ch - adhesion(node, cNew, sigma, M, N);
        
        % If this potential change results in a reduction in energy (Hnew
        % is negative) then we want to accept the change outright.
        % Otherwise we calculate a probability of changing, based on the
        % new energies of the system, should the change be accepted.
        if Hnew < 0
            pSwitch = 1;
        else
            pSwitch = P(Hnew, ch, cNew, T, q, N, node);
        end
        
        % Determine whether or not to make the switch and if so update the
        % community arrays appropriately
        if pSwitch > rand()
            if Switch(node, cNew) == 1
                UpdateActiveComList();
            end
        end
    end
    
    % Temperature update
    T = T/(1 + alpha * T);
    
    % Remove excess empty communities
    Rescale();
    
end

% Perform a post processing stage where by we iterate over all communities
% and try to place their nodes in every other community, taking all
% positive contributions to the modularity
bestNode = 0;
bestCom = 0;
bestEnergy = 0;
switched = 1;

while (switched == 1)
    switched = 0;
    for i = 1:N
        for j = 1:q
            if (sigma(i) ~= j)
                ch = cohesion (i, sigma, M, N);
                ah = adhesion (i, j, sigma, M, N);
                if (ch-ah < bestEnergy)
                    bestEnergy = ch-ah;
                    bestNode = i;
                    bestCom = j;
                end
            end
        end
        if (bestNode ~= 0 && bestCom ~= 0)
            if Switch(bestNode, bestCom) == 1
                UpdateActiveComList();
            end
            bestNode = 0;
            bestCom = 0;
            bestEnergy = 0;
            switched = 1;
        end
    end
    Rescale();
end



    function Rescale()
        
        % Determine the number of empty communities
        empty = size(find(comNodeCount==0),2);
        
        % Create a new community node list with one extra, empty community
        newSize = CN-empty+1;
        newComNC = zeros (1,newSize);
        
        % Resize the active community list appropriately
        acl = (1:1:newSize);
        q = newSize;
        
        % Populate the new community list
        n = 0;
        for i = 1:CN
            if comNodeCount(i)>0
                n = n + 1;
                newComNC(n) = comNodeCount(i);
                comNodeCount(i) = n;
            end
        end
        
        for i = 1:N
            sigma(i) = comNodeCount(sigma(i));
        end
        
        % reassign the new array and size
        CN = newSize;
        comNodeCount = newComNC;
        
        
    end




    function UpdateActiveComList()
        % If a community has been reduced to zero members or an empty community
        % has had a member added to it we need to update the list of active
        % communities appropriately
        n = 0;
        for i = 1:CN
            if (comNodeCount(i) > 0)
                n = n+1;
                acl(n) = i;
            end
        end
        
        
    end

    function c = Switch(node, cNew)
        % Move node from its current community into cNew and then determine whether
        % doing so has reduced it's current community to zero or started a new
        % community. If so, return 1, if not return 0.
        
        cOld = sigma(node);
        sigma(node) = cNew;
        
        % Update the node counts in the affected communities
        comNodeCount(cOld) = comNodeCount(cOld) - 1;
        comNodeCount(cNew) = comNodeCount(cNew) + 1;
        
        if comNodeCount(cOld) == 0 || comNodeCount(cNew) == 1
            c = 1;
        else
            c = 0;
        end
        
    end


    function p = P(Hnew, ch, cNew, T, q, N, node)
        % This basicaly calculates the energy of the system in the proposed
        % transition divided by the partition function, which in this case is
        % the sum of the energies of all the possible transition involving the
        % same node.
        
        % get the community to which the node currently belongs
        cCur = sigma(node);
        beta = (1.0/T);
        
        % The energy of the proposed transition
        Enew = exp(-beta*Hnew);
        
        % The total energy
        Etot = Enew;
        for i = 1:q
            if (acl(i) ~= cCur) && (acl(i) ~= cNew)
                Etot = Etot + exp(-beta * (ch - adhesion(node, acl(i), sigma, M, N)));
            end
        end
        
        p = Enew/Etot;
        
    end




end


function ch = cohesion (node, sigma, C, N)

ch = 0;

% get the community to which node belongs
curCom = sigma(node);

% calculate the correlation between node and every other node
% in the same community, curCom
for i = 1:N
    if sigma(i) == curCom
        ch = ch + C(node, i);
    end
end
% we don't want to count the correlation with itself, C(node, node) which
% will be added in the above sum when i = node.
ch = ch - C(node, node);

end


function ah = adhesion(node, cNew, sigma, C, N)

ah = 0;

% If node is in community cNew then there is no adhesion, only cohesion
if sigma(node) == cNew
    ah = 0;
else
    for i = 1:N
        if sigma(i) == cNew
            ah = ah + C(node, i);
        end
    end
end


% Since node will never be in cNew we don't have to worry about the
% addition of C(node, node) as we do for cohesion, since it can't happen.

end




