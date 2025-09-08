setwd("~/storage/Data/ms_evo_exwago/analysis/de_tables/ip_fam_tabs/")
pacman::p_load(ggplot2, dplyr)

files = list.files(pattern = "fam_prop.txt")
lst = lapply(files, read.delim)
names(lst) = sub(".fam_prop.txt", "", files)


lst2plot = list()
for (ip in names(lst)) {
  tmp.ip = lst[[ip]]
  tmp.ip$norm = (tmp.ip$raw_IP_prop+1)/(tmp.ip$fam_gprop+1)
  tmp.ip$IP = ip
  lst2plot[[ip]] = tmp.ip
}


nlst = lapply(lst2plot, function(x){
  sim_tbl = x[order(x$norm, decreasing = T),]
  sim_tbl = sim_tbl[!sim_tbl$Ave_CPM == 0 & !sim_tbl$raw_IP_prop == 0,]
  sim_tbl = sim_tbl[1:5,]
  sim_tbl$family = rownames(sim_tbl)
  return(sim_tbl)
  
})

names(nlst) = NULL
tmp.lts = do.call(rbind, nlst)
sim_fams = unique(tmp.lts$family)

lst2plot = lapply(lst2plot, function(x){
  nt = x[sim_fams,]
  rownames(nt) = sim_fams
  nt$raw_IP_prop = ifelse(is.na(nt$raw_IP_prop), 0, nt$raw_IP_prop)
  nt$norm = ifelse(is.na(nt$norm), 0, nt$norm)
  nt$IP = na.omit(nt$IP) %>% unique()
  return(nt)
})

df = do.call(rbind, lst2plot)
df$Feature = sub(".*\\.", "", rownames(df))
df$IP = sub("caenorhabditis_elegans_", "", df$IP)
df$IP = sub("sago1", "SAGO-1", df$IP)
df$IP = sub("sago2", "SAGO-2", df$IP)
df$IP = sub("ppw1", "PPW-1", df$IP)
df$IP = sub("heligmosomoides_bakeri", "H. bakeri", df$IP)
df$IP = sub("heligmosomoides_polygyrus", "H. polygyrus", df$IP)
df$IP = sub("nippostrongylus_brasiliensis", "N. brasiliensis", df$IP)
df$IP = sub("teladorsagia_circumcincta", "T. circumcincta", df$IP)
df$IP = sub("ancylostoma_ceylanicum", "A. ceylanicum", df$IP)
df$IP = factor(df$IP, levels = c("H. bakeri", "H. polygyrus", "N. brasiliensis", "T. circumcincta", "A. ceylanicum", "PPW-1", "SAGO-1", "SAGO-2"))

df

df$str = ifelse(grepl("_As", df$Feature), 1, ifelse(grepl("_S", df$Feature), 2, 3))
df$str = factor(df$str)



df$exp = "Relative enrcihment per superfamily"
norm = ggplot(df, aes(x=IP, y=norm %>% log2(),  group = Feature, color = Feature, shape = str )) +
  geom_line(colour=adjustcolor("black", alpha.f = .6), linetype="dashed") + 
  geom_point(aes(color=Feature, shape = str), size = 5) + theme_test() +
  facet_wrap(~exp)+
  scale_color_manual(values= gplots::rich.colors(20)) +
  theme(axis.title.x = element_text(face = "bold", size = 20), axis.text.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(face = "bold", size = 20),
        axis.text.y = element_text(size = 20, angle = 90, hjust = .5)) +
  xlab("") + ylab("log2(sRNA density in IP)") +
  theme(legend.position = "bottom") +
  theme( strip.text = element_text(size = 20, colour = "black")) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                           "inches")) + ylim(c(-6,4))

ggsave("nematode_ip_adult_fam_dens_paper.pdf", device = "pdf", plot = norm,
       path = "~/storage/Data/ms_evo_exwago/analysis/figures/", width = 14, height = 8, dpi = 300)  
