library(patchwork)
library(sf)
library(terra)
library(tidyterra)
library(tidyverse)
library(vegan)
theme_set(theme_classic())

cluster_gap <- function(mtx, max_k, seed = 139) {
  cluster_fn <- function(x, k) {
    list(cluster = stats::cutree(hclust(vegan::vegdist(x, method = "bray"),
                                        method = "ward.D2"),
                                 k = k))
  }
  set.seed(seed)
  gap <- suppressWarnings(
    cluster::clusGap(mtx,
                     FUN = cluster_fn,
                     K.max = 10,
                     B = 500,
                     verbose = FALSE)
  )
  as_tibble(gap$Tab) %>%
    mutate(
      k = seq_along(logW),
      is_optimum = k == cluster::maxSE(gap, SE.sim, method = "Tibs2001SEmax")
    )
}

stress <- function(dist, k) {
  tibble(
    k = k,
    stress = map_dbl(k, \(k) {
      vegan::metaMDS(
        dist,
        k = k,
        trymax = 100,
        trace = FALSE
      )$stress
    })
  )
}

nmds_to_df <- function(nmds, clust) {
  nmds$points %>%
    as_tibble(rownames = "station") %>%
    rename_with(\(axis) paste0("N", axis), starts_with("MDS")) %>%
    mutate(clust = factor(clust[station]))
}

zoop <- read_csv("data/Zoop_EsupTotal.csv") %>%
  sample_n(500)
zoop_mtx <- zoop %>%
  select(`Euphausia crystallorophias`:Gymnosomata) %>%
  as.matrix()
rownames(zoop_mtx) <- zoop$TowID
# pseudo_log
zoop_mtx <- log(zoop_mtx + 1)

zoop_dist <- vegdist(zoop_mtx, method = "bray")
zoop_clust <- hclust(zoop_dist, method = "ward.D2")
zoop_gap <- cluster_gap(zoop_mtx, max_k = 10, seed = 321)
ggplot(zoop_gap, aes(k, gap)) +
  geom_line() +
  geom_linerange(aes(ymin = gap - SE.sim, ymax = gap + SE.sim)) +
  geom_point(aes(color = is_optimum)) +
  expand_limits(x = 1, y = 0) +
  scale_color_manual(values = c("black", "red")) +
  scale_x_continuous(breaks = c(1, 4, 7, 10)) +
  labs(x = "Number of clusters", y = "Gap statistic") +
  theme(legend.position = "none")
zoop_cut <- cutree(zoop_clust, k = 2)
zoop_stress <- stress(zoop_mtx, k = 1:6)
ggplot(zoop_stress, aes(k, stress)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = c(0.2, 0.1, 0.05), linetype = "dotted") +
  scale_x_continuous(breaks = zoop_stress$k) +
  expand_limits(y = 0)
nmds_k <- 4
zoop_nmds <- vegan::metaMDS(zoop_dist, k = nmds_k, trymax = 100, trace = FALSE)
zoop_nmds_df <- nmds_to_df(zoop_nmds, zoop_cut)
zoop_nmds_plots <- map(seq(nmds_k - 1), \(i) {
  map(seq(i + 1, nmds_k), \(j) {
    axes <- paste0("NMDS", seq(nmds_k))
    ggplot(zoop_nmds_df, aes(.data[[axes[i]]], .data[[axes[j]]])) +
      geom_point(aes(color = clust), shape = 21, size = 0.5) +
      scale_color_brewer(palette = "Dark2") +
      coord_fixed() +
      labs(color = "Cluster") +
      theme(legend.position = "bottom")
  })
}) %>%
  unlist(recursive = FALSE)
reduce(zoop_nmds_plots, `+`) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")
# TODO: add envfit

# Make map
zoop_sp <- vect(mutate(zoop, cluster = factor(zoop_cut[TowID])),
                crs = "EPSG:4326",
                geom = c("LongitudeStart", "LatitudeStart"))
# Bathymetry contours
ibsco <- rast("data/nogithub/IBCSO_v2_bed.tif")
zoop_bbox_ibsco <- zoop_sp %>%
  project(crs(ibsco)) %>%
  ext() * 5
ibsco <- crop(ibsco, zoop_bbox_ibsco)
bathy_cont <- as.contour(ibsco, levels = c(-500, -1000, -2500)) %>%
  mutate(depth = factor(-level))
land <- as.polygons(ibsco >= 0) %>% filter(elevation == 1)
map_proj <- '+proj=aea +lon_0=-66 +lat_1=-68.25 +lat_2=-61.65 +lat_0=-64.95 +datum=WGS84 +units=m +no_defs'

map_bbox <- project(zoop_sp, map_proj) %>% ext() * 1.5
graticule <- st_graticule(lon = seq(-80, -40, by = 10),
                          lat = seq(-70, -58, by = 2)) %>%
  st_transform(map_proj)
ggplot() +
  geom_sf(data = graticule, color = "grey") +
  geom_spatvector(data = land, fill = "grey20", color = NA) +
  geom_spatvector(aes(linetype = depth), bathy_cont,
                  color = "navy", linewidth = 0.25, alpha = 0.75) +
  geom_spatvector(aes(color = cluster), zoop_sp) +
  scale_color_brewer(palette = "Dark2") +
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) +
  coord_sf(crs = sf::st_crs(map_proj),
           xlim = map_bbox[1:2], ylim = map_bbox[3:4],
           label_axes = "ENEN") +
  theme(legend.position = "bottom")
ggsave("scratch/lteramlrclusters.png")
