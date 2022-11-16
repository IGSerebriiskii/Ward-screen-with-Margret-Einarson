setwd("~/Work/Ward screen with Margret Einarson")
library(tidyverse)
library(readxl)
library(pheatmap)
fadu_data = read_excel("FADU_4_plates.xlsx",)
fadu_data %>% nrow()
fadu_wells = matrix(fadu_data$Well,16,24, byrow = T)

fadu_controls = c()
for (i in c(3,5,7,9,11,13)) {
          fadu_controls = c(fadu_controls, fadu_wells[i,5:9])
          fadu_controls = c(fadu_controls, fadu_wells[i,12:22])
}
fadu_controls

fadu_data %>% filter(Well %in% fadu_controls)
fadu_v1 = matrix(fadu_data$V_pl_1,16,24, byrow = T)
pheatmap(fadu_v1,cluster_rows = F,cluster_cols = F)
fadu_v2 = matrix(fadu_data$V_pl_2,16,24, byrow = T)
pheatmap(fadu_v2,cluster_rows = F,cluster_cols = F)
fadu_cbd1 = matrix(fadu_data$CBD_pl_1,16,24, byrow = T)
pheatmap(fadu_cbd1,cluster_rows = F,cluster_cols = F)
fadu_cbd2 = matrix(fadu_data$CBD_pl_2,16,24, byrow = T)
pheatmap(fadu_cbd2,cluster_rows = F,cluster_cols = F)

ggplot(fadu_data, aes(x = V_pl_1)) +
          geom_density()
fadu_v = fadu_data %>% select(Well,V_pl_1, V_pl_2)
fadu_v %>% pivot_longer(!Well, names_to = "plate", values_to = "CTB") %>% 
          ggplot(aes(x = CTB,fill = plate, col = plate)) +
          geom_density(alpha = 0.5)
fadu_data %>% pivot_longer(!Well, names_to = "plate", values_to = "CTB") %>% 
          ggplot(aes(x = CTB,fill = plate, col = plate)) +
          geom_density(alpha = 0.2)
ggsave("FADU CTB density by plates.png")

fadu_data_calc = fadu_data %>% mutate(veh = rowMeans(select(fadu_data,V_pl_1,V_pl_2)) , cbd = rowMeans(select(fadu_data,CBD_pl_1, CBD_pl_2)), 
                                      fr_change = cbd/veh, log_fr_change = round(log2(fr_change),3))
fadu_data_calc %>% write_csv("fadu_data_calc.csv")
fr_change = matrix(fadu_data_calc$fr_change,16,24, byrow = T)
pheatmap(fr_change,cluster_rows = F,cluster_cols = F)
log_fr_change = matrix(fadu_data_calc$log_fr_change,16,24, byrow = T)
pheatmap(log_fr_change,cluster_rows = F,cluster_cols = F)

fadu_data_calc %>% filter(Well %in% fadu_controls) %>% select(fr_change) %>% min()
fadu_data_calc %>% filter(Well %in% fadu_controls) %>% select(fr_change) %>% max()
fadu_data_calc %>% filter(Well %in% fadu_controls, fr_change> 1.4) %>% nrow()
fadu_data_calc %>% filter(Well %in% fadu_controls, fr_change < 0.7) %>% nrow()
fadu_data_calc %>% filter(Well %in% fadu_controls) %>% select(log_fr_change) %>% min()
fadu_data_calc %>% filter(Well %in% fadu_controls) %>% select(log_fr_change) %>% max()
fadu_data_calc %>% filter(Well %in% fadu_controls, log_fr_change> 0.5) %>% nrow()
fadu_data_calc %>% filter(Well %in% fadu_controls, log_fr_change < -0.5) %>% nrow()
fadu_controls %>% length() #96
fadu_inner_compounds = c()
for (i in c(2,4,6,8,10,12,14)) {
          fadu_inner_compounds = c(fadu_inner_compounds, fadu_wells[i,3:22])
          
}
fadu_inner_compounds %>% length() #140
fadu_data_calc %>% filter(Well %in% fadu_inner_compounds) %>% select(log_fr_change) %>% min()
fadu_data_calc %>% filter(Well %in% fadu_inner_compounds) %>% select(log_fr_change) %>% max()
# https://www.programmingr.com/statistics/z-score-in-r/