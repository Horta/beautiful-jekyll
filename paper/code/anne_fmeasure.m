function [fmeasure recall precision] = anne_fmeasure(biclustersa, biclustersb, n, p)
% O biclustersb eh considerado biclustering de referencia (golden truth).

Is = anne_intersection_size(biclustersa, biclustersb, n, p);

Usb = anne_union_size(biclustersb, biclustersb, n, p);
if Usb == 0
    recall = 0;
else
    recall = Is / Usb;
end

Usa = anne_union_size(biclustersa, biclustersa, n, p);
if Usa == 0
    precision = 0;
else
    precision = Is / Usa;
end

if precision+recall == 0
    fmeasure = 0;
else
    fmeasure = 2 * (precision*recall) / (precision+recall);
end

end