function [w ag dg] = csi(U, V)

% R. J. G. B. Campello, "Generalized external indexes for comparing data
% partitions with overlapping categories," Pattern Recogn. Lett., vol. 31,
% no. 9, pp. 966?975, 2010.

n = size(U,2);

clustersU = U2clusters(U);
clustersV = U2clusters(V);

alphaU = 0;

if matlabpool('size') >= 2
    
    spmd
        if labindex==1
            [alphaU betaU] = alphabeta(clustersU, n);
        elseif labindex==2
            [alphaV betaV] = alphabeta(clustersV, n);
        end
    end
    
    alphaU = alphaU{1};
    betaU = betaU{1};
    alphaV = alphaV{2};
    betaV = betaV{2};
else
    [alphaU betaU] = alphabeta(clustersU, n);
    [alphaV betaV] = alphabeta(clustersV, n);
end

Ag = min(alphaU, alphaV);
Dg = abs(alphaU - alphaV);

minbeta = min(betaU,betaV);
absbeta = abs(betaU-betaV);

for t = 1:n-1
    
    Ag(t,t+1:n) = Ag(t,t+1:n) + minbeta(t) + minbeta(t+1:n);
    Ag(t+1:n,t) = Ag(t,t+1:n);
    Dg(t,t+1:n) = Dg(t,t+1:n) + absbeta(t) + absbeta(t+1:n);
    Dg(t+1:n,t) = Dg(t,t+1:n);
end

ag=sum(sum(triu(Ag)));
dg=sum(sum(triu(Dg)));

if ag + dg == 0
    w = 1;
else
    w = ag / (ag+dg);
end

end

function [alpha beta] = alphabeta(clusters, n)

alpha = zeros(n,n);
beta = zeros(1,n);

k = length(clusters);

for i = 1:k
    
    for t = 1:length(clusters{i})-1
        
%         for j = t+1:length(clusters{i})
%             
%             alpha(clusters{i}(t), clusters{i}(j)) = alpha(clusters{i}(t), clusters{i}(j)) + 1;
%             alpha(clusters{i}(j), clusters{i}(t)) = alpha(clusters{i}(t), clusters{i}(j));
%             
%         end
        
        alpha(clusters{i}(t), clusters{i}(t+1:length(clusters{i}))) = alpha(clusters{i}(t), clusters{i}(t+1:length(clusters{i}))) + 1;
        alpha(clusters{i}(t+1:length(clusters{i})), clusters{i}(t)) = alpha(clusters{i}(t), clusters{i}(t+1:length(clusters{i})));
    end
    
    beta(clusters{i}) = beta(clusters{i}) + 1;
end

beta = beta - 1;

end
