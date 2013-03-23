
% +++ CONSTRUCTING EFFICIENT PORTFOLIO USING MEAN VARIANCE OPTIMIZATION +++
% -------------------------------------------------------------------------
%        %This program will use the Harry Markowitz MVO
% technique to % find the efficient portfolio for the user defined number 
% of rate of % returns that are equidistant between the lower and upper 
% bound, which is % again defined by the user at run time.
%
% NOTE:- Rmin and Rmax SHOULD be entered as decimal values. If you want to
% enter 8%, you SHOULD enter this value as 0.08
%
function [r,v,Z] = Markowitz(D, Rmin, Rmax, k, dotable, doplot)

% INPUT
% -----
% D == Matrix of size mxn
% Rmin, Rmax == user defined scalar bounds for rate of return
% k == number of values equidistantly placed between Rmin and Rmax
% dotable, doplotf == chart variables

n=length(D(:,1));
m=length(D(1,:));
S=ones(n,m);

% STEP 1: CONVERT THE RAW DATA INTO RATE OF RETURN MATRIX
% -------------------------------------------------------
% Convert the raw data matrix into rate of return for each asset i at time
% period t = 1,2,3,....,T

for i =2:n
    for t=1:m
        S(i,t) = (D(i,t)-D(i-1,t))/D(i-1,t);
    end
end

% STEP 2: COMPUTE THE COVARIANCE FROM RESULTANT MATRIX OF LAST ITERATION
% ----------------------------------------------------------------------
% Compute Covariance from the rate of return matrix

S(1,:) = [];
Q = cov(S);
n1=length(Q(:,1)); 
e = ones(n1,1);
f = linspace(Rmin,Rmax,k);
g = length(f);

% STEP 3: COMPUTE THE GEOMETRIC MEAN
% ----------------------------------
% Computing the Geometric mean

a1 = geo_mean(S(:,1)+1)-1;
a2 = geo_mean(S(:,2)+1)-1;
a3 = geo_mean(S(:,3)+1)-1;
mu = [a1;a2;a3];


% STEP 4: USE CVX SOLVER TO SOLVE THIS QUADRATIC PROGRAMMING
% ----------------------------------------------------------
% Repeatedly solve MVO for different mean return intervals

k=0;
for u=1:g
    %R = 0.065:0.005:0.105 -- If you use this value above, you will get the
    %exact results from the page 147 from the text book
    k=k+1;
    R=f(u);
    cvx_begin
    cvx_quiet(true)
    variable x(n1)
    minimize(x'*Q*x)
    subject to 
    mu'*x>=R; e'*x == 1; x >= 0;
    cvx_end
    r(k,1) = R; % A vector which can remember the return value 
    v(k,1) = cvx_optval; % remember the optimal variance
    Z(k,:) = x'; % Remember the optimal portfolio
end

% STEP 5: DISPLAY THE RESULTS IN TABLE AND PLOT IF GIVEN IN INPUT
% ---------------------------------------------------------------
% 	

if (dotable ==1)
    display('R      Variance   Stocks   Bonds   MM');
    display('-------------------------------------');
    for j=1:k
    display(sprintf('%1.3f   %1.4f   %1.2f   %1.2f   %1.2f',r(j),v(j),Z(j,1),Z(j,2),Z(j,3)));
    end
end

% Plot the efficient frontier
if (doplot ==1)
    hold on
    subplot(1,2,1)
    plot(sqrt(v)*100, r, '-ro')
    axis('square');
    title('Efficient Frontier')
    xlabel('Standard deviation (\%)')
    ylabel('Expected return (\%)')
    % Plot the composition of efficient portfolios
    subplot(1,2,2)
    area(r'*100,Z)
    axis('square')
    axis([r(1)*100, r(k)*100, 0, 1])
    title('Composition of efficient portfolios')
    xlabel('Expected return of efficient portfolios (\%)')
    ylabel('Percent invested in different asset classes')
    legend('Stocks','Bonds','Money Market')
    colormap summer
    hold off
end

% ... and that's all she wrote :-)
