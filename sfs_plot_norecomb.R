setwd("~/workdir/6910_project")

library(tidyverse)
library(ggplot2)
library(grid)
library(gridExtra)

# h=0 -------------------------------------------------------------------

slimid_map <- read_csv("slim_output/slimid_map.txt", 
                       col_types = cols(.default = "character")) %>%
  mutate_at(vars(-names(.)[1]), as.numeric)

plot_list <- list()
for (i in 1:9){
  slimid = slimid_map$slim_id[[i]]
  slimtxt = paste0("slim_output/chr6_run1_gen10000_m0_", slimid, ".txt")
  
  s=slimid_map$s[[i]]
  Ud=slimid_map$m1[[i]]*slimid_map$m[[i]]
  lambda = Ud*slimid_map$length[[i]]/s
  Ne = slimid_map$ne[[i]]*(2.71**(-lambda))
  
  Un=slimid_map$m0[[i]]*slimid_map$m[[i]]
  Un_per_pop = Un*slimid_map$length[[i]]*slimid_map$ne[[i]]
  
  #print(c(s, Ud, lambda, Ne, slimid_map$ne[[i]], Un, Un_per_pop))
  
  #sigma = sqrt(s*Ud)
  #1/(slimid_map$ne[[i]]*sigma)
  #1-1/(slimid_map$ne[[i]]*sigma)

  sfs <- read_delim(slimtxt, col_names = F) %>%
    select(X13) %>%
    rename(prevalence=X13) %>%
    #filter(prevalence!=1) %>%
    mutate(var_freq = prevalence/(2*slimid_map$ne[[i]])) %>%
    group_by(var_freq) %>%
    summarise(count = n()) %>%
    mutate(freq = count/sum(count)) %>%
    ungroup()

  total_expected = sum((2*slimid_map$ne[[i]]*Un_per_pop)/sfs$var_freq)
  #total_expected_Ne = sum((2*Ne*Un_per_pop)/sfs$var_freq)

  plot_list[[i]] <- sfs %>%
    ggplot(aes(x=var_freq, y=freq))+
    stat_smooth()+
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
  
  #print(2*Ne*Un_per_pop/(sfs$var_freq*total_expected))
}

plot_list_sub = list(plot_list[[3]],plot_list[[6]],plot_list[[9]])

grid.arrange(grobs=plot_list_sub, nrow=1, 
             top=textGrob("Site Frequency Spectrum, no recombination, length = chr6, N=1000, h=0", 
                          gp=gpar(fontsize=15,font=8)),
             left="Relative fraction of SNPs", 
             #left="Count", 
             bottom="Variant frequency (m0)")

plot_list_reorder = list(plot_list[[3]],plot_list[[6]],plot_list[[9]],
                         plot_list[[2]],plot_list[[5]],plot_list[[8]],
                         plot_list[[1]],plot_list[[4]],plot_list[[7]])
grid.arrange(grobs=plot_list_reorder, nrow=3, 
             top=textGrob("Site Frequency Spectrum, no recombination, length = chr6, N=1000, h=0", 
                          gp=gpar(fontsize=15,font=8)),
             left="Relative fraction of SNPs", 
             #left="Count", 
             bottom="Variant frequency")



# h=0.5 -------------------------------------------------------------------

slimid_map <- read_csv("slim_output/slimid_map_h0.5.txt", 
                       col_types = cols(.default = "character")) %>%
  mutate_at(vars(-names(.)[1]), as.numeric)

plot_list <- list()
for (i in 1:9){
  slimid = slimid_map$slim_id[[i]]
  slimtxt = paste0("slim_output/chr6_run1_gen10000_m0_", slimid, ".txt")
  
  s=slimid_map$s[[i]]
  Ud=slimid_map$m1[[i]]*slimid_map$m[[i]]
  lambda = Ud*slimid_map$length[[i]]/s
  Ne = slimid_map$ne[[i]]*(2.71**(-lambda))
  
  Un=slimid_map$m0[[i]]*slimid_map$m[[i]]
  Un_per_pop = Un*slimid_map$length[[i]]*slimid_map$ne[[i]]
  
  #print(c(s, Ud, lambda, Ne, slimid_map$ne[[i]], Un, Un_per_pop))
  
  #sigma = sqrt(s*Ud)
  #1/(slimid_map$ne[[i]]*sigma)
  #1-1/(slimid_map$ne[[i]]*sigma)
  
  sfs <- read_delim(slimtxt, col_names = F) %>%
    select(X13) %>%
    rename(prevalence=X13) %>%
    #filter(prevalence!=1) %>%
    mutate(var_freq = prevalence/(2*slimid_map$ne[[i]])) %>%
    group_by(var_freq) %>%
    summarise(count = n()) %>%
    mutate(freq = count/sum(count)) %>%
    ungroup()
  
  total_expected = sum((2*slimid_map$ne[[i]]*Un_per_pop)/sfs$var_freq)
  #total_expected_Ne = sum((2*Ne*Un_per_pop)/sfs$var_freq)
  
  plot_list[[i]] <- sfs %>%
    ggplot(aes(x=var_freq, y=freq))+
    stat_smooth()+
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
  
  #print(2*Ne*Un_per_pop/(sfs$var_freq*total_expected))
}

plot_list_sub = list(plot_list[[3]],plot_list[[6]],plot_list[[9]])

grid.arrange(grobs=plot_list_sub, nrow=1, 
             top=textGrob("Site Frequency Spectrum, no recombination, length = chr6, N=1000, h=0.5", 
                          gp=gpar(fontsize=15,font=8)),
             left="Relative fraction of SNPs", 
             #left="Count", 
             bottom="Variant frequency (m0)")

plot_list_reorder = list(plot_list[[3]],plot_list[[6]],plot_list[[9]],
                         plot_list[[2]],plot_list[[5]],plot_list[[8]],
                         plot_list[[1]],plot_list[[4]],plot_list[[7]])
grid.arrange(grobs=plot_list_reorder, nrow=3, 
             top=textGrob("Site Frequency Spectrum, no recombination, length = chr6, N=1000, h=0.5", 
                          gp=gpar(fontsize=15,font=8)),
             left="Relative fraction of SNPs", 
             #left="Count", 
             bottom="Variant frequency")
# 12x9 pdf landscape

