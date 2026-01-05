# from motherbase to divided by cluster

motherbase <- read.csv2("/Users/elviracastilloalmansa/Desktop/serriola_phenotypes/motherbase.csv")

## Morphocluster
mc_nr_clusters <- unique(motherbase$morpho_cluster)
for (i in mc_nr_clusters){
  sub <- motherbase[motherbase$morpho_cluster == i, c("LATITUDE","LONGITUDE")]
  sub <- na.omit(sub)
  write.csv(sub, paste0("/Users/elviracastilloalmansa/Desktop/serriola_phenotypes/coords_morphocluster_", i, ".csv"), row.names = FALSE)
  
  sub1 <- motherbase[motherbase$morpho_cluster == i, c("TKIID", "LKID")]
  write.csv(sub1, paste0("/Users/elviracastilloalmansa/Desktop/serriola_phenotypes/id_morphocluster_", i, ".csv"), row.names = FALSE)
}

## Morphotaxocluster 
mtc_nr_clusters <- unique(motherbase$morphotaxo_cluster)
for (i in mtc_nr_clusters){
  sub2 <- motherbase[motherbase$morphotaxo_cluster == i, c("LATITUDE","LONGITUDE")]
  sub2 <- na.omit(sub2)
  write.csv(sub2, paste0("/Users/elviracastilloalmansa/Desktop/serriola_phenotypes/coords_morphotaxocluster_", i, ".csv"), row.names = FALSE)
  
  sub3 <- motherbase[motherbase$morphotaxo_cluster == i, c("TKIID", "LKID")]
  write.csv(sub3, paste0("/Users/elviracastilloalmansa/Desktop/serriola_phenotypes/id_morphotaxocluster_", i, ".csv"), row.names = FALSE)
}
