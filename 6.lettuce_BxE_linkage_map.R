setwd("/Users/xt25/Desktop/BI/onemap/projects/lettuce/BxE_RIL_DArt_data/code")


# onemap 3.1.2 fix redundant marker bug
# devtools::install_github("Cristianetaniguti/onemap")
#install.packages("onemap")
library(onemap)

# read mapmaker file
lettuce_ril_onemap <- read_onemap(inputfile="../input/lettuce_onemap_input.raw")

# to produce a graphic with information about the raw data
plot(lettuce_ril_onemap)

# removes 4 individuals with high crossover/double crossover numbers
# Lac_0232, Lac_0243, Lac_0247, and Lac_0250
lettuce_ril <- remove_inds(lettuce_ril_onemap, c('Lac_0232', 'Lac_0243', 'Lac_0247', 'Lac_0250'))
lettuce_ril

# filter missing data
dat <- filter_missing(lettuce_ril, threshold = 0.1, by = "markers")
plot(dat)
dat1 <- filter_missing(dat, threshold = 0.1, by = "individuals")
plot(dat1)

# filter redundant markers
bins <- find_bins(dat1, exact = T)

(lettuce_bins <- create_data_bins(dat1, bins))

# exporting .raw file from OneMap object
write_onemap_raw(lettuce_bins, file.name = "../output/record/lettuce_new_dataset.raw")

# segregation tests
ril_test <- test_segregation(lettuce_bins)

(alpha <- Bonferroni_alpha(ril_test))

plot(ril_test)

# to show the markers numbers without segregation distortion
no_dist <- select_segreg(ril_test, distorted = FALSE, numbers = TRUE, threshold = alpha)
#length(no_dist)

# find an initial value to use for their linkage test
(LOD_sug <- suggest_lod(lettuce_bins))

twopts_ril <- rf_2pts(input.obj = lettuce_bins, LOD=LOD_sug, rm_mks = TRUE)

mark_no_dist_ril <- make_seq(twopts_ril, no_dist)
#LG1 <- make_seq(twopts_ril, no_dist[which(no_dist %in% which(twopts_ril$CHROM == "1"))])

LGs_ril <- group(mark_no_dist_ril, LOD = 5, max.rf = 0.4)


LGs_ril
set_map_fun(type = "kosambi")

groups <- LGs_ril$n.groups
LG_ril_record <- heatmaps_record <- vector("list", groups)

for(i in 1:groups){
  LG_temp_ril <- make_seq(LGs_ril, i)
  LG_ril_record[[i]] <- record(input.seq = LG_temp_ril, hmm = F)

  heatmaps_record[[i]] <- rf_graph_table(LG_ril_record[[i]], mrk.axis = "none", main = paste("LG",i," record map", sep=""))
}


heatmaps_record

# estimate the genetic distances
library(parallel)
(n.cores <-parallel::detectCores()-1) 
clust <- makeCluster(n.cores)
maps_record <- parLapply(clust, LG_ril_record, map, global_error=0.05)
stopCluster(clust)

# add redundant markers
maps_record_redundant <- vector("list", groups)
for(j in 1:groups){
  maps_record_redundant[[j]] <- add_redundants(maps_record[[j]], lettuce_ril, bins)
}

draw_map2(maps_record, output = "../output/record/map_record.png")
write_map(maps_record, "../output/record/Lettuce_BxE_RIL_map_record.map")
write_map(maps_record_redundant, "../output/record/Lettuce_BxE_RIL_map_record_redundant.map")


knitr::include_graphics("../output/record/map_record.png")
knitr::include_graphics("../output/record/lettuce_vs_phys_record.png")

(mbXcm <- plot_genome_vs_cm(maps_record, "kosambi", groups))