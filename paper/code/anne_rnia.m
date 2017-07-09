function r = anne_rnia(biclustersa, biclustersb, n, p)

Is = anne_intersection_size(biclustersa, biclustersb, n, p);
Us = anne_union_size(biclustersa, biclustersb, n, p);

r = 1 - ((Us - Is) / Us);

end
