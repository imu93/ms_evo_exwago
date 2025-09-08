setwd("~/storage/Data/ms_evo_exwago/analysis/de_tables/rrna_norm/cele_wagos/")
pacman::p_load(vegan, stringr, dplyr, edgeR, tibble, ggplot2, tidyr, kableExtra, ggpubr, patchwork)


# reda simplified topTable
files = list.files(pattern = "deTbl.txt")
files = files[!grepl("alg5", files)]
names(files) = c("A_ceylanicum", "alg1", "alg2", "alg3", "alg4", "csr1",
                 "ergo1", "hrde1", "nrde3", "ppw1", "ppw2", "prg1", "rde1",
                 "sago1", "sago2", "vsra1", "wago1", "wago10","wago4",
                 "H_bakeri","H_polygyrus", "N_brasiliensis", "T_circumcincta")
de_tabs = lapply(files, read.delim)

# Now Rds containing dge objects after IP vs input dea
dges = list.files("../../dge/cele_wagos/")
dges = dges[!grepl("alg5", dges)]
names(dges) = names(files)
dge_lst = lapply(paste0("../../dge/cele_wagos/",dges), readRDS)
names(dge_lst) = names(dges)
# I will estimate average CPM
ave_cpm = lapply(dge_lst, function(x){cpmByGroup(x, prior.count = 0.01)})
ave_ip_cpm = lapply(ave_cpm, function(x){x[,grepl("Ip|IP", colnames(x))]})

de_ip = list()
for (i in names(de_tabs)) {
  x = de_tabs[[i]]
  up = x[x$cl == 1,]
  tmp.cpm = ave_ip_cpm[[i]]
  up$cls = sub(":.*","", rownames(up))
  up$cls = up$cls %>%  str_replace("_rnd.*(_S|_As)", "\\1")
  up$cls = up$cls %>% str_replace("_rnd.*", "")
  up$cls = up$cls %>% str_replace("pseudo", "")
  up$cls = up$cls %>% str_replace("18s_|28s_|8s_", "")
  up$cls = up$cls %>% str_replace("^_", "")
  up$cls = up$cls %>% str_replace("_.*(_S|_As)", "\\1")
  up$cls = up$cls %>% str_replace("LSU_", "rRNA_")
  up$cls = ifelse(grepl("/", up$cls), sub("\\/.*_", "_", up$cls), up$cls)
  up$cls = ifelse(grepl("/", up$cls), sub("\\/.*", "", up$cls), up$cls)
  up$cls = ifelse(grepl("Unknown_", up$cls), sub("_.*", "", up$cls), up$cls)
  up$cls = ifelse(grepl("\\?", up$cls), sub("\\?", "", up$cls), up$cls)
  up$IP_cpm = tmp.cpm[match(rownames(up), names(tmp.cpm))]

  
  de_ip[[i]] = up
}



ids = lapply(de_ip, function(x){
  tmp.x = x
  tmp.x$cls %>% table %>% names
}) %>% Reduce(union,.)

counts= list()
for (i in names(de_ip)) {
  tmp.x = de_ip[[i]]
  tmp.x = split(tmp.x, tmp.x$cls) %>% lapply(function(x){x[,ncol(x)] %>% sum()}) %>% unlist()
  counts[[i]] = tmp.x
}

lst1 = list()
for (i in names(counts)) {
  tmp.x = counts[[i]]
  ms_ids = setdiff(ids, names(tmp.x))
  vls = rep(0, length(ms_ids))
  names(vls) = ms_ids
  tmp.x = c(tmp.x, vls)
  tmp.x = tmp.x[ids]
  lst1[[i]] = tmp.x 
}
df = do.call(rbind, lst1)

write.table(df, "agos_count_diversity.txt", sep = "\t", quote = F)


shannon = diversity(df, index = "shannon", MARGIN = 1) 
simpson = diversity(df, index = "simpson", MARGIN = 1) 

shannon[parasites] %>% mean()

index = rbind(shannon, simpson)
index = round(index,3)
index = index[,c(2:23, 1)]
index = index[,c(1:8, 10:12, 15:18, 9, 13,14, 19:23)]
index = as.data.frame(index)

index_clean <- index %>%
  rownames_to_column("Index")

index_long <- index_clean %>%
  pivot_longer(-Index, names_to = "AGO", values_to = "Diversity")

index_long$AGO <- recode(index_long$AGO,
                         "H_bakeri" = "italic('H. bakeri')",
                         "H_polygyrus" = "italic('H. polygyrus')",
                         "N_brasiliensis" = "italic('N. brasiliensis')",
                         "T_circumcincta" = "italic('T. circumcincta')",
                         "A_ceylanicum" = "italic('A. ceylanicum')")

index_long$exp = "Small RNA diversity"
p0 = ggplot(index_long, aes(x = reorder(AGO, Diversity), y = Diversity, fill = Index)) +
  geom_bar(stat = "identity", position = "dodge") + theme_test() +
  coord_flip() +
  scale_fill_manual( values = c("#FA3937", "#00BFFF") %>% adjustcolor(alpha.f = .6))+
  labs(x = "AGO", y = "Diversity Index", fill = "Index") +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  annotate("text", x = Inf, y = 1, label = "Simpson max", hjust = 1.1,
           vjust = -0.5, size = 3.5, color = "grey40") +
  facet_wrap(~exp) + theme(legend.position = "bottom") +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13),
    legend.position = "bottom",
    strip.text = element_text(size = 15))

ggsave("ago_diversity_no_rarefly.pdf", device = "pdf", width = 6, height = 7, dpi = 300, 
       path = "~/storage/Data/ms_evo_exwago/analysis/figures/", plot = p0)
  
# Define parasite names
parasites <- c("H_bakeri", "H_polygyrus", "N_brasiliensis", "T_circumcincta", "A_ceylanicum")

# Clean back AGO labels for group assignment
index_long <- index_clean %>%
  pivot_longer(-Index, names_to = "AGO", values_to = "Diversity") %>%
  mutate(
    Group = ifelse(AGO %in% parasites, "exWAGO (parasites)", "C. elegans AGOs")
  )


p1 <- index_long %>%
  filter(Index == "shannon") %>% mutate(exp = "Shannon Diversity") %>% 
  ggplot(aes(x = Group, y = Diversity, fill = Group)) +
  geom_boxplot(alpha = 0.8) +
  geom_point(position = position_jitter(width = 0.05, height = 0), 
             shape = 21, size = 2, stroke = 0.2) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  labs(x = "", y = "Index Value") +
  scale_fill_brewer(palette = "Set2") +
  theme_test(base_size = 14) +
  theme(legend.position = "none", 
        plot.title = element_text(size = 16, face = "bold")) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13),
    legend.position = "bottom",
    strip.text = element_text(size = 15)) +
  facet_wrap(~exp)


p2 <- index_long %>%
  filter(Index == "simpson") %>% mutate(exp = "Simpson Diversity")  %>%  
  ggplot(aes(x = Group, y = Diversity, fill = Group)) +
  geom_boxplot(alpha = 0.8) +
  geom_point(position = position_jitter(width = 0.05, height = 0), 
             shape = 21, size = 2, stroke = 0.2) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  labs(x = "", y = "Index Value") +
  scale_fill_brewer(palette = "Set2") +
  theme_test(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(size = 16, face = "bold")) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13),
    legend.position = "bottom",
    strip.text = element_text(size = 15)) +
  facet_wrap(~exp)


panel_A <- stak_prop +
  scale_x_discrete(labels = c(
    "H_bakeri" = expression(italic("H. bakeri")),
    "H_polygyrus" = expression(italic("H. polygyrus")),
    "N_brasiliensis" = expression(italic("N. brasiliensis")),
    "T_circumcincta" = expression(italic("T. circumcincta")),
    "A_ceylanicum" = expression(italic("A. ceylanicum")),
    "PPW-1" = "PPW-1",
    "SAGO-1" = "SAGO-1",
    "SAGO-2" = "SAGO-2"
  )) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))

top_row = ggarrange(panel_A, p0, labels = c("A", "B"), ncol = 2, widths = c(1.4, 1))

bottom_row = ggarrange(p1, p2, labels = c("C", "D"), ncol = 2, widths = c(.6, .6))

final_plot = ggarrange(top_row, bottom_row, nrow = 2, heights = c(1.2, 0.8), widths = c(1, .5), font.label = list(size = 18))



ggsave("ago_diversity.pdf", device = "pdf", width = 11, height = 13, dpi = 300, 
       path = "~/storage/Data/ms_evo_exwago/analysis/figures/", plot = final_plot)

kable(index_clean, format = "html", escape = FALSE, caption = "Diversity Indices for Argonaute Proteins") %>%
  kable_classic("striped", full_width = FALSE) %>%
  add_header_above(c(" " = 1, "C. elegans AGOs" = 18, "exWAGO (parasites)" = 5)) %>%
  column_spec(1, bold = TRUE)



cele_index = index[,!grepl("A_cey|H_ba|H_poly|N_bra|T_circ|alg|^rde", colnames(index))]
cele_index[2,] %>% range()
cele_index[2,] %>%  rowMeans()
cele_index[2,] %>% sd()

st_index = index[,grepl("A_cey|H_ba|H_poly|N_bra|T_circ", colnames(index))]
st_index[2,] %>% range()
st_index[2,] %>% rowMeans()
st_index[2,] %>% sd()



de_ip_wagos = de_ip[!grepl("alg|^rde|prg", names(de_ip))]
names(de_ip_wagos)
samp_size = lapply(de_ip_wagos, nrow) %>% unlist() %>% min()

set.seed(123)
counts= list()
for (i in names(de_ip_wagos)) {
  tmp.x = de_ip[[i]]
  tmp.x = tmp.x[sample(1:nrow(tmp.x) ,replace = F, size = samp_size),]
  tmp.x = split(tmp.x, tmp.x$cls) %>% lapply(function(x){x[,ncol(x)] %>% sum()}) %>% unlist()
  counts[[i]] = tmp.x
}


lst2 = list()
for (i in names(counts)) {
  tmp.x = counts[[i]]
  ms_ids = setdiff(ids, names(tmp.x))
  vls = rep(0, length(ms_ids))
  names(vls) = ms_ids
  tmp.x = c(tmp.x, vls)
  tmp.x = tmp.x[ids]
  lst2[[i]] = tmp.x 
}
df = do.call(rbind, lst2)

shannon = diversity(df, index = "shannon", MARGIN = 1) 
simpson = diversity(df, index = "simpson", MARGIN = 1) 

index2 = rbind(shannon, simpson)


rarefied_counts_list = replicate(100, {
  tmp = list()
  for (i in names(de_ip_wagos)) {
    tmp.x = de_ip[[i]]
    tmp.x = tmp.x[sample(nrow(tmp.x), size = samp_size), ]
    tmp.x = split(tmp.x, tmp.x$cls) %>% 
      lapply(\(x) sum(x$IP_cpm)) %>% 
      unlist()
    tmp[[i]] = tmp.x
  }
  # Standardize like lst2
  out = list()
  for (i in names(tmp)) {
    tx = tmp[[i]]
    vls = rep(0, length(ids)); names(vls) = ids
    out[[i]] = c(tx, vls[!names(vls) %in% names(tx)])[ids]
  }
  do.call(rbind, out)
}, simplify = FALSE)

# Calculate average diversity across iterations
shannon_vals = sapply(rarefied_counts_list, \(mat) diversity(mat, index = "shannon", MARGIN = 1))
simpson_vals = sapply(rarefied_counts_list, \(mat) diversity(mat, index = "simpson", MARGIN = 1))

index_shannon_mean = rowMeans(shannon_vals)
index_simpson_mean = rowMeans(simpson_vals)
index2_avg = rbind(shannon = index_shannon_mean, simpson = index_simpson_mean)


p02 = index2_avg %>% 
  as.data.frame() %>%
  rownames_to_column("Index") %>%
  pivot_longer(-Index, names_to = "AGO", values_to = "Diversity") %>% 
  mutate(exp="Small RNA diversity") %>% 
  ggplot( aes(x = reorder(AGO, Diversity), y = Diversity, fill = Index)) +
  geom_bar(stat = "identity", position = "dodge") + theme_test() +
  coord_flip() +
  labs(x = "AGO", y = "Diversity Index", fill = "Index") +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  annotate("text", x = Inf, y = 1, label = "Simpson max", hjust = 1.1,
           vjust = -0.5, size = 3.5, color = "grey40") +
  facet_wrap(~exp) + theme(legend.position = "bottom") +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13),
    legend.position = "bottom",
    strip.text = element_text(size = 15))


index_long2 = index2_avg %>% as.data.frame() %>% 
  rownames_to_column("Index") %>%
  pivot_longer(-Index, names_to = "AGO", values_to = "Diversity") %>%
  mutate(
    Group = ifelse(AGO %in% parasites, "exWAGO (parasites)", "C. elegans AGOs")
  )


p21 = index_long2 %>%
  filter(Index == "shannon") %>% mutate(exp = "Shannon Diversity") %>% 
  ggplot(aes(x = Group, y = Diversity, fill = Group)) +
  geom_boxplot(alpha = 0.8, outliers = F) +
  geom_point(position = position_jitter(width = 0.05, height = 0), 
             shape = 21, size = 2, stroke = 0.2) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  labs(x = "", y = "Index Value") +
  scale_fill_brewer(palette = "Set2") +
  theme_test(base_size = 14) +
  theme(legend.position = "none", 
        plot.title = element_text(size = 16, face = "bold")) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13),
    legend.position = "bottom",
    strip.text = element_text(size = 15)) +
  facet_wrap(~exp)


p22 <- index_long2 %>%
  filter(Index == "simpson") %>% mutate(exp = "Simpson Diversity")  %>%  
  ggplot(aes(x = Group, y = Diversity, fill = Group)) +
  geom_boxplot(alpha = 0.8, outliers = F) +
  geom_point(position = position_jitter(width = 0.05, height = 0), 
             shape = 21, size = 2, stroke = 0.2) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  labs(x = "", y = "Index Value") +
  scale_fill_brewer(palette = "Set2") +
  theme_test(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(size = 16, face = "bold")) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13),
    legend.position = "bottom",
    strip.text = element_text(size = 15)) +
  facet_wrap(~exp)


top_row = ggarrange(panel_A, p02, labels = c("A", "B"), ncol = 2, widths = c(1.4, 1))

bottom_row = ggarrange(p21, p22, labels = c("C", "D"), ncol = 2, widths = c(.6, .6))

final_plot = ggarrange(top_row, bottom_row, nrow = 2, heights = c(1.2, 0.8), widths = c(1, .5), font.label = list(size = 18))


ggsave("raref_ago_diversity.pdf", device = "pdf", width = 11, height = 13, dpi = 300, 
       path = "~/storage/Data/ms_evo_exwago/analysis/figures/", plot = final_plot)

stro = index2[,grepl("A_cey|H_bake|H_poly|N_bra|T_circ", colnames(index2))]
Cele = index2[,!grepl("A_cey|H_bake|H_poly|N_bra|T_circ", colnames(index2))]

wilcox.test(stro[1,], Cele[1,])
wilcox.test(stro[2,], Cele[2,])
