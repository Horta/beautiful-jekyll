function w=csi_fast(ba, bb, n, p)

pa = biclusters2pclusters(ba, n, p);
pa = sort_clusters(pa);

pb = biclusters2pclusters(bb, n, p);
pb = sort_clusters(pb);

cc = collide_clusters([pa pb]);

cpa = cell(1,length(pa));
cpb = cell(1,length(pb));
mult = zeros(1,length(cc));
for i = 1:length(cc)
    
    obj = min(cc{i});
    mult(i) = length(cc{i});
    cpa = insert_obj(cpa, pa, i, obj);
    cpb = insert_obj(cpb, pb, i, obj);
end

[w ag dg] = csi(cpa, cpb, mult, n*p);
end

function clusters = sort_clusters(clusters)

for i = 1:length(clusters)
    clusters{i} = sort(clusters{i});
end

end

function cpp = insert_obj(cpp, pp, index, obj)

for i = 1:length(pp)
    if ismember(obj, pp{i})
        
        cpp{i} = [cpp{i} index];
    end
end

end

function qout = collide_clusters(pp)

qin = pp;
qout = {};

while ~isempty(qin)
    
    spot = qin{end};
    qin = qin(1:end-1);
    
    qaux = {};
    for i = 1:length(qout)
        
        res = intersect_sorted(spot, qout{i});
        if isempty(res)
            qaux = [qaux qout(i)];
        else
            qaux = [qaux res];
            spot = setdiff_sorted(spot, res);
            qoutres = setdiff_sorted(qout{i}, res);
            if ~isempty(qoutres)
                qaux = [qaux qoutres];
            end
        end
        if isempty(spot)
            break
        end
    end
    if ~isempty(spot)
        qaux = [qaux spot];
    end
    qout = qaux;
end

end

function [w ag dg] = csi(cpa, cpb, mult, ntotal)

n = length(mult);

[alphaU, betaU] = alphabeta(cpa, n);
[alphaV, betaV] = alphabeta(cpb, n);

Ag = min(alphaU, alphaV);
Dg = abs(alphaU - alphaV);

minbeta = min(betaU,betaV);
absbeta = abs(betaU-betaV);

for t = 1:n
    
    Ag(t,t:n) = Ag(t,t:n) + minbeta(t) + minbeta(t:n);
    Ag(t:n,t) = Ag(t,t:n);
    Dg(t,t:n) = Dg(t,t:n) + absbeta(t) + absbeta(t:n);
    Dg(t:n,t) = Dg(t,t:n);
end

for i = 1:n-1
    for j = i+1:n
        Ag(i,j) = Ag(i,j)*mult(i)*mult(j);
        Dg(i,j) = Dg(i,j)*mult(i)*mult(j);
    end
end

for i = 1:n
    if mult(i)==1
        Ag(i,i) = 0;
        Dg(i,i) = 0;
    else
        Ag(i,i) = Ag(i,i)*nchoosek(mult(i),2);
        Dg(i,i) = Dg(i,i)*nchoosek(mult(i),2);
    end
end

ag=sum(sum(triu(Ag)));
dg=sum(sum(triu(Dg)));

ag = ag + (ntotal-sum(mult))*sum(mult.*minbeta);
dg = dg + (ntotal-sum(mult))*sum(mult.*absbeta);


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
    
    for t = 1:length(clusters{i})
        
        alpha(clusters{i}(t), clusters{i}(t:length(clusters{i}))) = alpha(clusters{i}(t), clusters{i}(t:length(clusters{i}))) + 1;
        alpha(clusters{i}(t:length(clusters{i})), clusters{i}(t)) = alpha(clusters{i}(t), clusters{i}(t:length(clusters{i})));
    end
    
    beta(clusters{i}) = beta(clusters{i}) + 1;
end

beta = beta - 1;
beta(beta==-1) = 0; % to handle singletons not in |clusters|

end

function iv = intersect_sorted(veca, vecb)

iv = veca( ismembc(veca, vecb) );

end

function dv = setdiff_sorted(veca, vecb)

dv = veca( ~ismembc(veca, vecb) );

end
