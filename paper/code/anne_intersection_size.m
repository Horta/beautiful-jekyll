function Is = anne_intersection_size(biclustersa, biclustersb, n, p)

Na = anne_num_cluster_belonging_to(biclustersa, n, p);
Nb = anne_num_cluster_belonging_to(biclustersb, n, p);

Is = sum(sum(min(Na,Nb)));

end
