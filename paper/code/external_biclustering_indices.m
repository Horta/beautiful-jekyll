function r = external_biclustering_indices(which, ba, bb, nr, nc)

if isempty(ba) && isempty(bb)
    r = 1;
    return
end

if isempty(ba) || isempty(bb)
    r = 0;
    return
end

%indices = {'prel' 'prec' 'rnia' 'ce' 'l&w' 'stm' 'wjac' 'wdic' 'fabia'...
%           'u' 'e' 'ay' 'erel' 'erec' 'csi' 'ebc'};

if strcmp(which,'prel')
    
    r = applyPrelicRelevance(ba, bb);
    
elseif strcmp(which,'prec')
    
    r = applyPrelicRecovery(ba, bb);

elseif strcmp(which,'rnia')
    
    r = anne_rnia(ba, bb, nr, nc);

elseif strcmp(which,'ce')
   
    r = biclusteringError(ba, bb, nr, nc);
    
elseif strcmp(which,'l&w')
    
    r = apply_lew(ba, bb);
    
elseif strcmp(which,'stm')
    
    [r assoc] = apply_stmaria(ba, bb, nr, nc);
    
elseif strcmp(which,'wjac')
    
    [r assoc] = apply_wjac(ba, bb, nr, nc);
    
elseif strcmp(which,'wdic')
    
    [r assoc] = apply_wdic(ba, bb, nr, nc);
    
elseif strcmp(which,'fabia')
    
    [r assoc] = apply_fabia(ba, bb, nr, nc);

elseif strcmp(which,'u')
    
    r = applyBozdagUncovered(ba, bb, nr, nc);
    
elseif strcmp(which,'e')
    
    r = applyBozdagExtra(ba, bb, nr, nc);
    
elseif strcmp(which,'ay')
    
    r = applyAyadi(ba, bb);
        
elseif strcmp(which,'erel')
    
    r = apply_erenRelev(ba, bb, nr, nc);
    
elseif strcmp(which,'erec')
    
    r = apply_erenRecov(ba, bb, nr, nc);
    
elseif strcmp(which,'naive-csi')
    
    r = applyNaiveCSI(ba, bb, nr, nc);

elseif strcmp(which,'fast-csi')
    
    r = applyFastCSI(ba, bb, nr, nc);
    
elseif strcmp(which,'ebc')
    
    r = applyEBC(ba, bb, nr, nc);
    
else

    assert(false, 'Unrecognized measure');
end

end

function r = applyEBC(ba, bb, nr, nc)

U = biclusters2UBackground(ba, nr, nc);
V = biclusters2UBackground(bb, nr, nc);

r = exbcubed(U, V);

end

function r = applyAyadi(ba, bb)

us = inline('length(union(a,b))');
is = inline('length(intersect(a,b))');

r = 0;
for i = 1:length(ba)
    
    maior = -inf;
    for j = 1:length(bb)
        
        top    = is(ba(i).rows,bb(j).rows)*is(ba(i).cols,bb(j).cols);
        bottom = us(ba(i).rows,bb(j).rows)*us(ba(i).cols,bb(j).cols);
        
        maior = max([maior top/bottom]);
    end
    r = r + maior;
end

r = r/length(ba);

end

function r = applyBozdagUncovered(ba, bb, nr, nc)

pba = biclusters2pclusters(ba, nr, nc);
pbb = biclusters2pclusters(bb, nr, nc);

uniona = unique( [pba{:}] );
unionb = unique( [pbb{:}] );

r = length(intersect(uniona,unionb))/length(unionb);

end

function r = applyBozdagExtra(ba, bb, nr, nc)

r = applyBozdagUncovered(bb, ba, nr, nc);

end


function r = applyPrelicRelevance(ba, bb)

r = sqrt(applyPrelicGenes(ba, bb) * applyPrelicChips(ba, bb));

end

function r = applyPrelicRecovery(ba, bb)

r = applyPrelicRelevance(bb, ba);

end

function r = applyPrelicGenes(ba, bb)

r = 0;
for i = 1:length(ba)
    
    rbest = -inf;
    for j = 1:length(bb)
        
        rbest = max([rbest length(intersect( ba(i).rows, bb(j).rows )) / length(union( ba(i).rows, bb(j).rows ))]);
    end
    r = r + rbest;
end

r = r/length(ba);

end

function r = applyPrelicChips(ba, bb)

r = 0;
for i = 1:length(ba)
    
    rbest = -inf;
    for j = 1:length(bb)
        
        rbest = max([rbest length(intersect( ba(i).cols, bb(j).cols )) / length(union( ba(i).cols, bb(j).cols ))]);
    end
    r = r + rbest;
end

r = r/length(ba);

end

function r = applyNaiveCSI(ba, bb, nr, nc)

Ua = biclusters2UBackground(ba, nr, nc);
Ub = biclusters2UBackground(bb, nr, nc);

[r,~,~] = csi(Ua, Ub);
end

function r = applyFastCSI(ba, bb, nr, nc)

r = fast_csi(ba, bb, nr, nc);
end

function r = apply_lew(ba, bb)

k = length(ba);
q = length(bb);

S = zeros(k,q);

for i = 1:k
    
    for j = 1:q
        
        rrint = intersect(ba(i).rows, bb(j).rows);
        ccint = intersect(ba(i).cols, bb(j).cols);
        
        rruni = union(ba(i).rows, bb(j).rows);
        ccuni = union(ba(i).cols, bb(j).cols);
        
        S(i,j) = (length(rrint) + length(ccint)) / (length(rruni) + length(ccuni));
    end
end

[vals ~] = max(S,[],2);
r = sum(vals) / k;

end

function [r assoc] = apply_stmaria(ba, bb, nr, nc)

D = calculate_dice(ba, bb, nr, nc);

[vals indices] = max(D,[],2);
r = sum(vals) / length(ba);

assoc = [(1:length(ba))' indices];

end

function [r assoc] = apply_wjac(ba, bb, nr, nc)

J = calculate_jaccard(ba, bb, nr, nc);
pclustersa = biclusters2pclusters(ba, nr, nc);

k = length(ba);

up = 0;
down = 0;

assoc = zeros(k,2);

for i = 1:k
    
    [val index] = max(J(i,:));
    
    up = up + length(pclustersa{i}) * val;
    
    assoc(i,1) = i; assoc(i,2) = index;
    
    down = down + length(pclustersa{i});
end

r = up / down;

end

function [r assoc] = apply_wdic(ba, bb, nr, nc)

D = calculate_dice(ba, bb, nr, nc);
pclustersa = biclusters2pclusters(ba, nr, nc);

k = length(ba);

up = 0;
down = 0;

assoc = zeros(k,2);

for i = 1:k
    
    [val index] = max(D(i,:));
    up = up + length(pclustersa{i}) * val;
    
    down = down + length(pclustersa{i});
    
    assoc(i,1) = i; assoc(i,2) = index;
end

r = up / down;

end

function [r assoc] = apply_fabia(ba, bb, nr, nc)

J = calculate_jaccard(ba, bb, nr, nc);

k = length(ba);
q = length(bb);

assignment = munkres(-1*J);
[rows,cols,~] = find(assignment);

r = 0;
for i = 1:min(k,q)
    r = r + J(rows(i),cols(i));
end

r = r / max(k,q);

assoc = [rows cols];

end

function D = calculate_dice(ba, bb, nr, nc)

pclustersa = biclusters2pclusters(ba, nr, nc);
pclustersb = biclusters2pclusters(bb, nr, nc);

k = length(ba);
q = length(bb);

D = zeros(k,q);

for i = 1:k
    
    for j = 1:q
        
        s = 2*length(intersect(pclustersa{i}, pclustersb{j}));
        D(i,j) = s / (length(pclustersa{i}) + length(pclustersb{j}));
    end
end

end

function J = calculate_jaccard(ba, bb, nr, nc)

pclustersa = biclusters2pclusters(ba, nr, nc);
pclustersb = biclusters2pclusters(bb, nr, nc);

k = length(ba);
q = length(bb);

J = zeros(k,q);

for i = 1:k
    
    for j = 1:q
        
        int = length(intersect(pclustersa{i}, pclustersb{j}));
        uni = length(union(pclustersa{i}, pclustersb{j}));
        
        J(i,j) = int / uni;
    end
end

end

function r = apply_erenRelev(ba, bb, nr, nc)

pclustersa = biclusters2pclusters(ba, nr, nc);
pclustersb = biclusters2pclusters(bb, nr, nc);

r = 0;
for i = 1:length(pclustersa)
    
    best = -inf;
    for j = 1:length(pclustersb)
        
        int = length(intersect(pclustersa{i}, pclustersb{j}));
        uni = length(union(pclustersa{i}, pclustersb{j}));
        t = int / uni;
        
        if t > best
            best = t;
        end
    end
    r = r + best;
end

r = r / length(pclustersa);

end

function r = apply_erenRecov(ba, bb, nr, nc)

r = apply_erenRelev(bb, ba, nr, nc);

end
