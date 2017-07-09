function Us = anne_union_size(biclustersa, biclustersb, n, p)

Na = anne_num_cluster_belonging_to(biclustersa, n, p);
Nb = anne_num_cluster_belonging_to(biclustersb, n, p);

Us = sum(sum(max(Na,Nb)));

end
