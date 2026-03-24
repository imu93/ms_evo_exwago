setwd("/home/isaac/backup/mac/Data/ms_exwago_evo/analysis/repeat_comparison/")
pacman::p_load(dplyr, Biostrings,  ggplot2, ggpubr, ggdist)

wbpFiles = list.files(path = "./wbp_libs/",pattern = ".*.fa$") 
curatedFiles = list.files(path = "./curated_libs/", pattern = ".*.fa$")

wbp_libs = lapply(paste0("./wbp_libs/", wbpFiles), readDNAStringSet)
cur_libs = lapply(paste0("./curated_libs/", curatedFiles), readDNAStringSet)

names(wbp_libs) = sub("\\..*","", wbpFiles)
names(cur_libs) = sub("\\..*","", curatedFiles)



info_lib = function(x){
  nat_lib_fa = x[grepl("^rnd|^ltr", names(x)),]
  names(nat_lib_fa) = sub(" .*", "", names(nat_lib_fa))
  wdt_nat = nat_lib_fa %>% width() 
  df_nat = wdt_nat %>% data.frame()
  rownames(df_nat) = names(nat_lib_fa)
  df_nat$class = sub(".*#", "", rownames(df_nat)) %>%  sub("\\/.*", "", .)
  colnames(df_nat) = sub("\\.", "length", colnames(df_nat))
  df_nat = df_nat[grepl("DNA|RC|LINE|LTR|PLE|SINE|Unk", df_nat$class),]
  df_nat$class = ifelse(grepl("RC", rownames(df_nat)), "DNA", df_nat$class)
  df_nat$class = ifelse(grepl("PLE", rownames(df_nat)), "LINE", df_nat$class)
  df_nat$class = sub("\\?", "", df_nat$class)
  return(df_nat)
}


species_list = list()
for (i in names(cur_libs)) {
  print(i)
  tmp.wbp = wbp_libs[[i]]
  tmp.cur = cur_libs[[i]]
  wbp_df = info_lib(tmp.wbp)
  cur_df = info_lib(tmp.cur)
  wbp_df$set = "WBP"
  wbp_df$species = i
  wbp_df$class = sub(" .*", "", wbp_df$class)
  cur_df$set = "Curated"
  cur_df$species = i
  ndf = rbind(cur_df, wbp_df)
  species_list[[i]] = ndf
}
names(species_list) = NULL
df = do.call(rbind, species_list)
df$class %>% table()
df$class = factor(df$class, levels = c("DNA", "LINE", "LTR", "Unknown", "SINE"))
df$species = factor(df$species, levels = c("heligmosomoides_bakeri", "nippostrongylus_brasiliensis",
                                           "teladorsagia_circumcincta", "ancylostoma_ceylanicum"))
df$set = factor(df$set, levels = c("WBP", "Curated"))



# split violin 
GeomSplitViolin <- ggproto(
  "GeomSplitViolin", 
  GeomViolin, 
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data, 
                      xminv = x - violinwidth * (x - xmin), 
                      xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1,'group']
    newdata <- plyr::arrange(
      transform(data, x = if(grp%%2==1) xminv else xmaxv), 
      if(grp%%2==1) y else -y
    )
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin", 
                       grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    } else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

geom_split_violin <- function (mapping = NULL, 
                               data = NULL, 
                               stat = "ydensity", 
                               position = "identity", ..., 
                               draw_quantiles = NULL, 
                               trim = TRUE, 
                               scale = "area", 
                               na.rm = FALSE, 
                               show.legend = NA, 
                               inherit.aes = TRUE) {
  layer(data = data, 
        mapping = mapping, 
        stat = stat, 
        geom = GeomSplitViolin, 
        position = position, 
        show.legend = show.legend, 
        inherit.aes = inherit.aes, 
        params = list(trim = trim, 
                      scale = scale, 
                      draw_quantiles = draw_quantiles, 
                      na.rm = na.rm, ...)
  )
}


te_data = df %>% split(df$set) %>% lapply(function(cls){split(cls,cls$species) %>% lapply(function(sp){split(sp, sp$class)}) }) %>%  lapply(., function(lvl1) {
  lapply(lvl1, function(lvl2) {
    sapply(lvl2, nrow)
  })
}) 


# Step 1: Convert your nested list into a data frame
extract_counts <- function(data, source) {
  do.call(rbind, lapply(names(data), function(species) {
    df <- as.data.frame(t(data[[species]]))
    df$species <- species
    df$set <- source
    df$class <- rownames(df)
    colnames(df)[1] <- "Count"
    return(df)
  }))
}

# Extract counts from both WBP and Curated datasets
df_wbp = extract_counts(te_data$WBP, "WBP")
df_curated = extract_counts(te_data$Curated, "Curated")
counts_df = rbind(df_wbp, df_curated)




rs = reshape2::melt(counts_df)
rs = rs[,c(1,2,4,5)]
rs$variable = sub("Count", "DNA", rs$variable)
colnames(rs) = sub("variable","class", colnames(rs))
rs$class = factor(rs$class, levels = c("DNA","LINE", "LTR", "Unknown", "SINE"))
rs$species = sub("heligmosomoides_bakeri", "H. bakeri", rs$species)
rs$species = sub("nippostrongylus_brasiliensis", "N. brasiliensis", rs$species)
rs$species = sub("teladorsagia_circumcincta", "T. circumcincta", rs$species)
rs$species = sub("ancylostoma_ceylanicum", "A. ceylanicum", rs$species)


rs$species = factor(rs$species,levels=c("H. bakeri",
                                        "N. brasiliensis",
                                        "T. circumcincta",
                                        "A. ceylanicum"))
rs_wbp = rs[rs$set == "WBP",]
rs_cur = rs[rs$set == "Curated",]
tbl_wbp =  rs_wbp %>% tibble() %>% group_by(species, set)
tbl_cur =  rs_cur %>% tibble() %>% group_by(species, set)

df$species = as.character(df$species)
df$species =  sub("heligmosomoides_bakeri", "H. bakeri", df$species)
df$species =  sub("nippostrongylus_brasiliensis", "N. brasiliensis", df$species)
df$species =  sub("teladorsagia_circumcincta", "T. circumcincta", df$species)
df$species =  sub("ancylostoma_ceylanicum", "A. ceylanicum", df$species)

df$species = factor(df$species,levels=c("H. bakeri",
                                        "N. brasiliensis",
                                        "T. circumcincta",
                                        "A. ceylanicum"))




p1 = df %>% ggplot(aes(x = class, y = log10(length), fill = set)) +
  geom_split_violin(linewidth = .001,  scale = "width") +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
               geom = "point", 
               position = position_dodge(width = .25)
  ) +
  theme_test() +
  xlab("TE class") + ylab(bquote(~log[10]~ 'Consensus length')) +
  scale_fill_manual(values = c("#16CAF5","#E3014F") %>% 
                      adjustcolor(alpha.f = .8), name=NULL, label=c("Uncurated", "Curated")) +
  theme(axis.title.x = element_text(size = 18), 
        axis.text.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(face = "bold", size = 18),
        axis.text.y = element_text(size = 16)) +
  facet_wrap(~species, ncol=1) + 
  theme( strip.text = element_text(size = 14, face = "italic", colour = "black")) +
  theme(legend.position = "bottom") +
  theme(legend.text = element_text(size = 14)) +
  guides(color = guide_legend(override.aes = list(size=10))) 
  

df_tmp = df[!df$class == "SINE",]
results <- df_tmp %>%
  filter(set %in% c("Curated", "WBP")) %>%  # Filter for the sets of interest
  group_by(species, class) %>%  # Group by species and class
  summarise(p_value = wilcox.test(
    log10(length[set == "Curated"]), 
    log10(length[set == "WBP"]), 
    alternative = "greater"
  )$p.value)

p1 +  geom_text(
  data = results,
  aes(x = class, y = 5.4, label = ifelse( p_value < 0.05, "***", "ns")),
  inherit.aes = FALSE,
  position = position_dodge(width = 0.25),
  size = 5
) + 
geom_text(data = tbl_wbp, 
              aes(x = class, y=4.8, label = value),
            size = 4) +
  geom_text(data = tbl_cur, 
            aes(x = class, y=4.4, label = value),
            size = 4)


ggsave("Fig1A_GBE.pdf", device = "pdf", 
       width = 5, height = 10, dpi = 300, plot = last_plot(),
       path = "~/storage/Data/ms_evo_exwago/analysis/figures/")
