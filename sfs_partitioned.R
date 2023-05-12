setwd("~/workdir/6910_project")

library(tidyverse)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)

slimid_map <- read_csv("slim_output_recomb/slim_map.txt", 
                       col_types = cols(.default = "character")) %>%
  mutate_at(vars(-names(.)[1]), as.numeric)

slimid_map$se = c("1e-1","1e-1","1e-1","1e-2","1e-2","1e-2","1e-3","1e-3","1e-3")

plot_list <- list()
for (i in 1:9){
  slimid = slimid_map$slim_id[[i]]
  m0_slimtxt = paste0("slim_output_recomb/chr6_run1_gen10000_m0_s",slimid_map$se[[i]],"_",slimid, ".txt")
  m1_slimtxt = paste0("slim_output_recomb/chr6_run1_gen10000_m1_s",slimid_map$se[[i]],"_",slimid, ".txt")
  
  m1 <- read_delim(m1_slimtxt, col_names = F) %>% 
    select(X5,X6,X7,X12) %>% 
    rename(id=X5, type=X6, prevalence=X12, location=X7) %>% 
    mutate(min_range_50kb = location - 50000, max_range_50kb = location + 50000) %>% 
    mutate(min_range_10kb = location - 10000, max_range_10kb = location + 10000) 
  
  m0 <- read_delim(m0_slimtxt, col_names = F) %>% 
    select(X5,X6,X7,X12) %>% 
    rename(id=X5, type=X6, prevalence=X12, location=X7) %>% 
    mutate(var_freq = prevalence/(2*1000))
  
  filter_range_50kb <- function(min_range_50kb, max_range_50kb, locations) {
    locations[locations %in% seq(min_range_50kb, max_range_50kb)] #inside 10kb
  }
  
  filtered_locations_50kb <- pmap(m1[, c("min_range_50kb", "max_range_50kb")], filter_range_50kb, m0$location) %>% 
    reduce(union) #within 50 kb
  
  filter_range_10kb <- function(min_range_10kb, max_range_10kb, locations) {
    locations[locations %in% seq(min_range_10kb, max_range_10kb)] #inside 10kb
  }
  
  filtered_locations_10kb <- pmap(m1[, c("min_range_10kb", "max_range_10kb")], filter_range_10kb, filtered_locations_50kb) %>% 
    reduce(union) #within 10 kb
  
  filtered_locations_10kb_50kb <- setdiff(filtered_locations_50kb, filtered_locations_10kb)
  # between 10kb and 50kb
  
  filtered_m0_10kb <- m0 %>% 
    filter(location %in% filtered_locations_10kb) %>% 
    group_by(var_freq) %>% 
    summarise(count = n()) %>% 
    mutate(freq = count/sum(count)) %>% 
    ungroup() %>% 
    mutate(region = "<10kb")
  
  filtered_m0_outside_50kb <- m0 %>%
    filter(!location %in% filtered_locations_50kb) %>% 
    group_by(var_freq) %>% 
    summarise(count = n()) %>% 
    mutate(freq = count/sum(count)) %>% 
    ungroup() %>% 
    mutate(region = ">50kb")
  
  filtered_m0_10kb_50kb <- m0 %>% 
    filter(location %in% filtered_locations_10kb_50kb) %>% 
    group_by(var_freq) %>% 
    summarise(count = n()) %>% 
    mutate(freq = count/sum(count)) %>% 
    ungroup() %>% 
    mutate(region = "10-50kb")
  
  m0_partitioned <- rbind(filtered_m0_10kb,filtered_m0_10kb_50kb,filtered_m0_outside_50kb)
  
  s=slimid_map$s[[i]]
  Ud=slimid_map$m1[[i]]*slimid_map$m[[i]]
  lambda = Ud*slimid_map$length[[i]]/s
  Ne = slimid_map$ne[[i]]*(2.71**(-lambda))
  
  Un=slimid_map$m0[[i]]*slimid_map$m[[i]]
  Un_per_pop = Un*slimid_map$length[[i]]*slimid_map$ne[[i]]
  total_expected = sum((2*slimid_map$ne[[i]]*Un_per_pop)/m0_partitioned$var_freq)
  
  m0_partitioned$region <- factor(m0_partitioned$region, levels = c("<10kb", "10-50kb", ">50kb"))
  
  plot_list[[i]] <- m0_partitioned %>% 
    ggplot(aes(x=var_freq, y=freq, group = region))+
    stat_smooth(aes(color = region))+
    scale_color_manual(values = c("#F8766D","#00BA38","#619CFF"))+
    geom_line(aes(x=var_freq, y=2*slimid_map$ne[[i]]*Un_per_pop/(var_freq*total_expected)), linetype="solid") +
    #geom_line(aes(x=var_freq, y=2*Ne*Un_per_pop/(var_freq*total_expected)), linetype="dashed") +
    labs(title=paste0("s=", slimid_map$s[[i]],
                      ", Ud=", sprintf("%0.3f", round(Ud*slimid_map$length[[i]], digits = 3)), 
                      ", Un=", sprintf("%0.3f", round(Un*slimid_map$length[[i]], digits = 3))),
         x="",y="")+
    scale_x_continuous(trans="logit",
                       breaks = c(0.00001, 0.001, 0.1, 0.5,1-0.1, 1-0.001,1-0.005)) +
    scale_y_log10() +
    theme_classic()
}

plot_list_reorder = list(plot_list[[3]],plot_list[[6]],plot_list[[9]],
                         plot_list[[2]],plot_list[[5]],plot_list[[8]],
                         plot_list[[1]],plot_list[[4]],plot_list[[7]])

figure <- ggarrange(plotlist = plot_list_reorder, ncol=3, nrow=3, common.legend = TRUE, legend="right")

grid.arrange(figure, 
            top=textGrob("Site Frequency Spectrum, recombination, length = chr6, N=1000, h=0", 
                         gp=gpar(fontsize=15,font=8)),
            left="Relative fraction of SNPs", 
            #left="Count", 
            bottom="Variant frequency")


# h=0.5 -------------------------------------------------------------------

slimid_map <- read_csv("slim_output_recomb/slimid_map_h0.5.txt", 
                       col_types = cols(.default = "character")) %>%
  mutate_at(vars(-names(.)[1]), as.numeric)

slimid_map$se = c("1e-1","1e-1","1e-1","1e-2","1e-2","1e-2","1e-3","1e-3","1e-3")

plot_list_h0.5 <- list()
for (i in 1:9){
  slimid = slimid_map$slim_id[[i]]
  m0_slimtxt = paste0("slim_output_recomb/chr6_run1_gen10000_m0_s",slimid_map$se[[i]],"_",slimid, ".txt")
  m1_slimtxt = paste0("slim_output_recomb/chr6_run1_gen10000_m1_s",slimid_map$se[[i]],"_",slimid, ".txt")
  
  m1 <- read_delim(m1_slimtxt, col_names = F) %>% 
    select(X5,X6,X7,X12) %>% 
    rename(id=X5, type=X6, prevalence=X12, location=X7) %>% 
    mutate(min_range_50kb = location - 50000, max_range_50kb = location + 50000) %>% 
    mutate(min_range_10kb = location - 10000, max_range_10kb = location + 10000) 
  
  m0 <- read_delim(m0_slimtxt, col_names = F) %>% 
    select(X5,X6,X7,X12) %>% 
    rename(id=X5, type=X6, prevalence=X12, location=X7) %>% 
    mutate(var_freq = prevalence/(2*1000))
  
  filter_range_50kb <- function(min_range_50kb, max_range_50kb, locations) {
    locations[locations %in% seq(min_range_50kb, max_range_50kb)] #inside 10kb
  }
  
  filtered_locations_50kb <- pmap(m1[, c("min_range_50kb", "max_range_50kb")], filter_range_50kb, m0$location) %>% 
    reduce(union) #within 50 kb
  
  filter_range_10kb <- function(min_range_10kb, max_range_10kb, locations) {
    locations[locations %in% seq(min_range_10kb, max_range_10kb)] #inside 10kb
  }
  
  filtered_locations_10kb <- pmap(m1[, c("min_range_10kb", "max_range_10kb")], filter_range_10kb, filtered_locations_50kb) %>% 
    reduce(union) #within 10 kb
  
  filtered_locations_10kb_50kb <- setdiff(filtered_locations_50kb, filtered_locations_10kb)
  # between 10kb and 50kb
  
  filtered_m0_10kb <- m0 %>% 
    filter(location %in% filtered_locations_10kb) %>% 
    group_by(var_freq) %>% 
    summarise(count = n()) %>% 
    mutate(freq = count/sum(count)) %>% 
    ungroup() %>% 
    mutate(region = "<10kb")
  
  filtered_m0_outside_50kb <- m0 %>%
    filter(!location %in% filtered_locations_50kb) %>% 
    group_by(var_freq) %>% 
    summarise(count = n()) %>% 
    mutate(freq = count/sum(count)) %>% 
    ungroup() %>% 
    mutate(region = ">50kb")
  
  filtered_m0_10kb_50kb <- m0 %>% 
    filter(location %in% filtered_locations_10kb_50kb) %>% 
    group_by(var_freq) %>% 
    summarise(count = n()) %>% 
    mutate(freq = count/sum(count)) %>% 
    ungroup() %>% 
    mutate(region = "10-50kb")
  
  m0_partitioned <- rbind(filtered_m0_10kb,filtered_m0_10kb_50kb,filtered_m0_outside_50kb)
  
  s=slimid_map$s[[i]]
  Ud=slimid_map$m1[[i]]*slimid_map$m[[i]]
  lambda = Ud*slimid_map$length[[i]]/s
  Ne = slimid_map$ne[[i]]*(2.71**(-lambda))
  
  Un=slimid_map$m0[[i]]*slimid_map$m[[i]]
  Un_per_pop = Un*slimid_map$length[[i]]*slimid_map$ne[[i]]
  total_expected = sum((2*slimid_map$ne[[i]]*Un_per_pop)/m0_partitioned$var_freq)
  
  m0_partitioned$region <- factor(m0_partitioned$region, levels = c("<10kb", "10-50kb", ">50kb"))
  
  plot_list_h0.5[[i]] <- m0_partitioned %>% 
    ggplot(aes(x=var_freq, y=freq, group = region))+
    geom_smooth(aes(color = region))+
    scale_color_manual(values = c("#F8766D","#00BA38","#619CFF"))+
    geom_line(aes(x=var_freq, y=2*slimid_map$ne[[i]]*Un_per_pop/(var_freq*total_expected)), linetype="solid") +
    #geom_line(aes(x=var_freq, y=2*Ne*Un_per_pop/(var_freq*total_expected)), linetype="dashed") +
    labs(title=paste0("s=", slimid_map$s[[i]],
                      ", Ud=", sprintf("%0.3f", round(Ud*slimid_map$length[[i]], digits = 3)), 
                      ", Un=", sprintf("%0.3f", round(Un*slimid_map$length[[i]], digits = 3))),
         x="",y="")+
    scale_x_continuous(trans="logit",
                       breaks = c(0.00001, 0.001, 0.1, 0.5,1-0.1, 1-0.001,1-0.005)) +
    scale_y_log10() +
    theme_classic()
}

plot_list_h0.5_reorder = list(plot_list_h0.5[[3]],plot_list_h0.5[[6]],plot_list_h0.5[[9]],
                         plot_list_h0.5[[2]],plot_list_h0.5[[5]],plot_list_h0.5[[8]],
                         plot_list_h0.5[[1]],plot_list_h0.5[[4]],plot_list_h0.5[[7]])

figure <- ggarrange(plotlist = plot_list_h0.5_reorder, ncol=3, nrow=3, common.legend = TRUE, legend="right")

grid.arrange(figure, 
             top=textGrob("Site Frequency Spectrum, recombination, length = chr6, N=1000, h=0.5", 
                          gp=gpar(fontsize=15,font=8)),
             left="Relative fraction of SNPs", 
             #left="Count", 
             bottom="Variant frequency")
