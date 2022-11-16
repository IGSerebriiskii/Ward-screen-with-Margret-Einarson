setwd("~/Work/Ward screen with Margret Einarson")
library(tidyverse)
library(readxl)
library(pheatmap)
fadu_data = read_excel("FADU_4_plates_run_1.xlsx")
fadu_data %>% nrow()
fadu_wells = matrix(fadu_data$Well,16,24, byrow = T)


# Defining control wells
fadu_controls = c()
for (i in c(3,5,7,9,11,13)) {
          fadu_controls = c(fadu_controls, fadu_wells[i,5:9])
          fadu_controls = c(fadu_controls, fadu_wells[i,12:22])
}
fadu_controls %>% length() #96
fadu_controls %>% as_tibble_col(column_name = "control_wells") %>%  write_csv("control_wells.csv")

# Defining compound wells - 
fadu_inner_compounds = c()
for (i in c(2,4,6,8,10,12,14)) {
          fadu_inner_compounds = c(fadu_inner_compounds, fadu_wells[i,3:22])
          
}
fadu_inner_compounds %>% length() #140
fadu_inner_compounds %>% as_tibble_col(column_name = "compound_wells") %>%  write_csv("inner_compound_wells.csv")

fadu_gem = c()
for (i in c(3,5,7,9,11,13)) {
          fadu_gem = c(fadu_gem, fadu_wells[i,c(4,11)])
}
fadu_gem %>% length() #96
fadu_gem %>% as_tibble_col(column_name = "gem_wells") %>%  write_csv("gem_wells.csv")

fadu_doc = c()
for (i in c(3,5,7,9,11,13)) {
          fadu_doc = c(fadu_doc, fadu_wells[i,c(3,10)])
}
fadu_doc %>% length() #96
fadu_doc %>% as_tibble_col(column_name = "doc_wells") %>%  write_csv("doc_wells.csv")

# - - - - - - - - - - - - - - - - - - -
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
# same as above, a bit shorter
fadu_data_calc = fadu_data %>% mutate(veh = rowMeans(.[2:3]) , cbd = rowMeans(.[4:5]), 
                                      fr_change = cbd/veh, log_fr_change = round(log2(fr_change),3))

fadu_data_calc %>% write_csv("fadu_data_calc.csv")
fr_change = matrix(fadu_data_calc$fr_change,16,24, byrow = T)
pheatmap(fr_change,cluster_rows = F,cluster_cols = F)
log_fr_change = matrix(fadu_data_calc$log_fr_change,16,24, byrow = T)
pheatmap(log_fr_change,cluster_rows = F,cluster_cols = F)

fadu_data_calc %>% filter(Well %in% fadu_controls) %>% select(fr_change) %>% min()
fadu_data_calc %>% filter(Well %in% fadu_controls) %>% select(fr_change) %>% max()
fadu_data_calc %>% filter(Well %in% fadu_controls, fr_change> 1.28) %>% nrow() # 2
fadu_data_calc %>% filter(Well %in% fadu_controls, fr_change < 0.76) %>% nrow() # 1
fadu_data_calc %>% filter(Well %in% fadu_controls) %>% select(log_fr_change) %>% min()
fadu_data_calc %>% filter(Well %in% fadu_controls) %>% select(log_fr_change) %>% max()
fadu_data_calc %>% filter(Well %in% fadu_controls, log_fr_change> 0.5) %>% nrow()
fadu_data_calc %>% filter(Well %in% fadu_controls, log_fr_change < -0.5) %>% nrow()

fadu_data_calc %>% filter(Well %in% fadu_inner_compounds) %>% select(log_fr_change) %>% min()
fadu_data_calc %>% filter(Well %in% fadu_inner_compounds) %>% select(log_fr_change) %>% max()
# https://www.programmingr.com/statistics/z-score-in-r/

fadu_data_calc %>% filter(Well %in% fadu_controls) %>% pull(fr_change) %>% shapiro.test()
# and it is NOT a normal distribution although looks not too distorted
fadu_data_calc %>% filter(Well %in% fadu_controls) %>% 
  ggplot(aes(x = fr_change)) + geom_density()
fadu_data_calc %>% filter(Well %in% fadu_controls) %>% pull(fr_change) %>% mean()
fadu_data_calc %>% filter(Well %in% fadu_controls) %>% pull(fr_change) %>% sd()
fadu_data_calc %>% filter(Well %in% fadu_controls) %>% pull(veh) %>% mean()


fadu_data_calc %>% filter(Well %in% fadu_inner_compounds) %>% select(fr_change) %>% min()
fadu_data_calc %>% filter(Well %in% fadu_inner_compounds) %>% select(fr_change) %>% max()

fadu_hits = fadu_data_calc %>% filter(Well %in% fadu_inner_compounds,(fr_change > 1.27844 | fr_change < 0.7594916)) 
fadu_hits %>% write_csv("fadu_hits.csv")


fadu_data %>% filter(Well %in% fadu_controls) %>% select(V_pl_1,V_pl_2) %>% cor() #0.3645661
fadu_data %>% filter(Well %in% fadu_controls) %>% select(CBD_pl_1,CBD_pl_2) %>% cor() # 0.3401433
fadu_data %>% filter(Well %in% fadu_inner_compounds) %>% select(V_pl_1,V_pl_2) %>% cor() #0.7407735
fadu_data %>% filter(Well %in% fadu_inner_compounds) %>% select(CBD_pl_1,CBD_pl_2) %>% cor() # 0.6920442
fadu_data %>% filter(Well %in% fadu_controls) %>% 
  ggplot(aes(x = V_pl_1, y = V_pl_2))+
  geom_point()
fadu_data %>% filter(Well %in% fadu_inner_compounds) %>% 
  ggplot(aes(x = V_pl_1, y = V_pl_2))+
  geom_point()

fadu_data %>% filter(Well %in% fadu_controls) %>% 
  ggplot(aes(x = V_pl_1, y = V_pl_2))+
  geom_point()
fadu_data %>% filter(Well %in% fadu_inner_compounds) %>% 
  ggplot(aes(x = V_pl_1, y = CBD_pl_1))+
  geom_point()
fadu_data_calc %>% filter(Well %in% fadu_inner_compounds) %>% 
  ggplot(aes(x = veh, y = cbd))+
  geom_point()
fadu_data_calc %>% filter(Well %in% fadu_inner_compounds) %>% 
  ggplot(aes(x = veh, y = cbd))+
  geom_point() + geom_smooth()
fadu_data_calc %>% filter(Well %in% fadu_inner_compounds) %>% 
          ggplot(aes(x = veh, y = cbd))+
          geom_point() + geom_smooth(method = 'lm')
fadu_data_calc_for_plot_comp = fadu_data_calc %>% filter(Well %in% fadu_inner_compounds) %>% 
          select(veh,cbd) %>% mutate(well = "compounds")
fadu_data_calc_for_plot_control = fadu_data_calc %>% filter(Well %in% fadu_controls) %>% 
          select(veh,cbd) %>% mutate(well = "controls")
fadu_data_calc_for_plot_gem = fadu_data_calc %>% filter(Well %in% fadu_gem) %>% 
          select(veh,cbd) %>% mutate(well = "gem")
fadu_data_calc_for_plot_doc = fadu_data_calc %>% filter(Well %in% fadu_doc) %>% 
          select(veh,cbd) %>% mutate(well = "doc")

fadu_data_calc_for_plot = bind_rows(fadu_data_calc_for_plot_control,
                                    fadu_data_calc_for_plot_comp,
                                    fadu_data_calc_for_plot_gem,
                                    fadu_data_calc_for_plot_doc)
fadu_data_calc_for_plot  %>% 
          ggplot(aes(x = veh, y = cbd, color = well))+
          geom_point() + geom_smooth(method = 'lm')

#  * * *
fadu_data2 = read_excel("FADU_6_plates_run_2.xlsx")
fadu_data2 %>% nrow() # 384 which is correct



# fadu_v1 = matrix(fadu_data$V_pl_1,16,24, byrow = T)
# pheatmap(fadu_v1,cluster_rows = F,cluster_cols = F)
# fadu_v2 = matrix(fadu_data$V_pl_2,16,24, byrow = T)
# pheatmap(fadu_v2,cluster_rows = F,cluster_cols = F)
# fadu_cbd1 = matrix(fadu_data$CBD_pl_1,16,24, byrow = T)
# pheatmap(fadu_cbd1,cluster_rows = F,cluster_cols = F)
# fadu_cbd2 = matrix(fadu_data$CBD_pl_2,16,24, byrow = T)
# pheatmap(fadu_cbd2,cluster_rows = F,cluster_cols = F)
# 
# ggplot(fadu_data, aes(x = V_pl_1)) +
#   geom_density()
# fadu_v = fadu_data %>% select(Well,V_pl_1, V_pl_2)
# fadu_v %>% pivot_longer(!Well, names_to = "plate", values_to = "CTB") %>% 
#   ggplot(aes(x = CTB,fill = plate, col = plate)) +
#   geom_density(alpha = 0.5)
# fadu_data %>% pivot_longer(!Well, names_to = "plate", values_to = "CTB") %>% 
#   ggplot(aes(x = CTB,fill = plate, col = plate)) +
#   geom_density(alpha = 0.2)
# ggsave("FADU CTB density by plates.png")

# fadu_data2_calc = fadu_data2 %>% mutate(veh = rowMeans(select(fadu_data,V_pl_1,V_pl_2)) , cbd = rowMeans(select(fadu_data,CBD_pl_1, CBD_pl_2)), 
#                                       fr_change = cbd/veh, log_fr_change = round(log2(fr_change),3))
# same as above, a bit shorter
fadu_data2_calc = fadu_data2 %>% mutate(veh = rowMeans(.[2:4]) , cbd = rowMeans(.[5:7]), 
                                      fr_change = cbd/veh, log_fr_change = round(log2(fr_change),3))

fadu_data2_calc %>% write_csv("fadu_data2_calc.csv")
fr_change2 = matrix(fadu_data2_calc$fr_change,16,24, byrow = T)
pheatmap(fr_change2,cluster_rows = F,cluster_cols = F)
# log_fr_change = matrix(fadu_data_calc$log_fr_change,16,24, byrow = T)
# pheatmap(log_fr_change,cluster_rows = F,cluster_cols = F)

fadu_data2_calc %>% filter(Well %in% fadu_controls) %>% select(fr_change) %>% min() # 0.6357736
fadu_data2_calc %>% filter(Well %in% fadu_controls) %>% select(fr_change) %>% max() # 1.136606
fadu_data2_calc %>% filter(Well %in% fadu_controls, fr_change> 1.28) %>% nrow() # 2
fadu_data2_calc %>% filter(Well %in% fadu_controls, fr_change < 0.76) %>% nrow() # 1
fadu_data2_calc %>% filter(Well %in% fadu_controls) %>% select(log_fr_change) %>% min() # -0.653
fadu_data2_calc %>% filter(Well %in% fadu_controls) %>% select(log_fr_change) %>% max() #0.185
fadu_data2_calc %>% filter(Well %in% fadu_controls, log_fr_change> 0.5) %>% nrow()
fadu_data2_calc %>% filter(Well %in% fadu_controls, log_fr_change < -0.5) %>% nrow()
# fadu_controls %>% length() #96
# fadu_inner_compounds = c()
# for (i in c(2,4,6,8,10,12,14)) {
#   fadu_inner_compounds = c(fadu_inner_compounds, fadu_wells[i,3:22])
#   
# }
# fadu_inner_compounds %>% length() #140
fadu_data2_calc %>% filter(Well %in% fadu_inner_compounds) %>% select(log_fr_change) %>% min() # -1.132
fadu_data2_calc %>% filter(Well %in% fadu_inner_compounds) %>% select(log_fr_change) %>% max()
# https://www.programmingr.com/statistics/z-score-in-r/

fadu_data2_calc %>% filter(Well %in% fadu_controls) %>% pull(fr_change) %>% shapiro.test()
# normally distributed
fadu_data2_calc %>% filter(Well %in% fadu_controls) %>% 
  ggplot(aes(x = fr_change)) + geom_density()
fadu2_controls_mean = fadu_data2_calc %>% filter(Well %in% fadu_controls) %>% pull(fr_change) %>% mean()
fadu2_controls_sd = fadu_data2_calc %>% filter(Well %in% fadu_controls) %>% pull(fr_change) %>% sd()
# fadu_data_calc %>% filter(Well %in% fadu_controls) %>% pull(veh) %>% mean()
fadu2_controls_mean
fadu2_controls_sd
fadu_data2_calc %>% filter(Well %in% fadu_inner_compounds) %>% select(fr_change) %>% min() 0.456178
fadu_data2_calc %>% filter(Well %in% fadu_inner_compounds) %>% select(fr_change) %>% max() 1.544866

fadu_hits2 = fadu_data2_calc %>% filter(Well %in% fadu_inner_compounds,(fr_change > (fadu2_controls_mean + 2* fadu2_controls_sd)| fr_change < (fadu2_controls_mean - 2* fadu2_controls_sd))) 
fadu_hits2 %>% write_csv("fadu_hits2.csv")
# fadu_data %>% filter(Well %in% fadu_controls) %>% select(V_pl_1,V_pl_2) %>% cor() #0.3645661
# fadu_data %>% filter(Well %in% fadu_controls) %>% select(CBD_pl_1,CBD_pl_2) %>% cor() # 0.3401433
# fadu_data %>% filter(Well %in% fadu_inner_compounds) %>% select(V_pl_1,V_pl_2) %>% cor() #0.7407735
# fadu_data %>% filter(Well %in% fadu_inner_compounds) %>% select(CBD_pl_1,CBD_pl_2) %>% cor() # 0.6920442
# fadu_data %>% filter(Well %in% fadu_controls) %>% 
#   ggplot(aes(x = V_pl_1, y = V_pl_2))+
#   geom_point()
# fadu_data %>% filter(Well %in% fadu_inner_compounds) %>% 
#   ggplot(aes(x = V_pl_1, y = V_pl_2))+
#   geom_point()
# 
# fadu_data %>% filter(Well %in% fadu_controls) %>% 
#   ggplot(aes(x = V_pl_1, y = V_pl_2))+
#   geom_point()
# fadu_data %>% filter(Well %in% fadu_inner_compounds) %>% 
#   ggplot(aes(x = V_pl_1, y = CBD_pl_1))+
#   geom_point()
fadu_data2_calc %>% filter(Well %in% fadu_inner_compounds) %>% 
  ggplot(aes(x = veh, y = cbd))+
  geom_point()
fadu_data2_calc %>% filter(Well %in% fadu_inner_compounds) %>% 
  ggplot(aes(x = veh, y = cbd))+
  geom_point() + geom_smooth(method = "lm")
ggsave("veh vs cbd run2.png")

fadu_data2_calc_for_plot_comp = fadu_data2_calc %>% filter(Well %in% fadu_inner_compounds) %>% 
          select(veh,cbd) %>% mutate(well = "compounds")
fadu_data2_calc_for_plot_control = fadu_data2_calc %>% filter(Well %in% fadu_controls) %>% 
          select(veh,cbd) %>% mutate(well = "controls")
fadu_data2_calc_for_plot_gem = fadu_data2_calc %>% filter(Well %in% fadu_gem) %>% 
          select(veh,cbd) %>% mutate(well = "gem")
fadu_data2_calc_for_plot_doc = fadu_data2_calc %>% filter(Well %in% fadu_doc) %>% 
          select(veh,cbd) %>% mutate(well = "doc")

fadu_data2_calc_for_plot = bind_rows(fadu_data2_calc_for_plot_control,
                                    fadu_data2_calc_for_plot_comp,
                                    fadu_data2_calc_for_plot_gem,
                                    fadu_data2_calc_for_plot_doc)
fadu_data2_calc_for_plot  %>% 
          ggplot(aes(x = veh, y = cbd, color = well))+
          geom_point() + geom_smooth(method = 'lm')
ggsave("veh vs cbd run2 detailed.png")


# - - - - - - - - - - - - - -
#  * * *
fadu_data2n = read_excel("FADU run2 nuclei count.xlsx")
fadu_data2n %>% nrow() # 2304  which is 384 * 6 which is correct
fadu_data2n = fadu_data2n %>% pivot_wider(names_from = Plate, values_from = Nuclei)


# fadu_v1 = matrix(fadu_data$V_pl_1,16,24, byrow = T)
# pheatmap(fadu_v1,cluster_rows = F,cluster_cols = F)
# fadu_v2 = matrix(fadu_data$V_pl_2,16,24, byrow = T)
# pheatmap(fadu_v2,cluster_rows = F,cluster_cols = F)
# fadu_cbd1 = matrix(fadu_data$CBD_pl_1,16,24, byrow = T)
# pheatmap(fadu_cbd1,cluster_rows = F,cluster_cols = F)
# fadu_cbd2 = matrix(fadu_data$CBD_pl_2,16,24, byrow = T)
# pheatmap(fadu_cbd2,cluster_rows = F,cluster_cols = F)
# 
# ggplot(fadu_data, aes(x = V_pl_1)) +
#   geom_density()
# fadu_v = fadu_data %>% select(Well,V_pl_1, V_pl_2)
# fadu_v %>% pivot_longer(!Well, names_to = "plate", values_to = "CTB") %>% 
#   ggplot(aes(x = CTB,fill = plate, col = plate)) +
#   geom_density(alpha = 0.5)
# fadu_data %>% pivot_longer(!Well, names_to = "plate", values_to = "CTB") %>% 
#   ggplot(aes(x = CTB,fill = plate, col = plate)) +
#   geom_density(alpha = 0.2)
# ggsave("FADU CTB density by plates.png")

# fadu_data2_calc = fadu_data2 %>% mutate(veh = rowMeans(select(fadu_data,V_pl_1,V_pl_2)) , cbd = rowMeans(select(fadu_data,CBD_pl_1, CBD_pl_2)), 
#                                       fr_change = cbd/veh, log_fr_change = round(log2(fr_change),3))
# same as above, a bit shorter
fadu_data2n_calc = fadu_data2n %>% mutate(veh = rowMeans(.[5:7]) , cbd = rowMeans(.[2:4]), 
fr_change = cbd/veh, log_fr_change = round(log2(fr_change),3))
# note the column order is different
# 
# fadu_data2n_calc = fadu_data2n %>% rowwise() %>% mutate(veh = mean(c_across(5:7)) , cbd = mean(c_across(2:4)),
#                                                         veh_sd = sd(c_across(5:7)) , cbd_sd = sd(c_across(2:4)),
#                        fr_change = cbd/veh, log_fr_change = round(log2(fr_change),3))
# yet a third way

# maybe the easiest way - on the top
fadu_data2n_calc$veh_sd = apply(fadu_data2n_calc[5:7],1,sd)
fadu_data2n_calc$cbd_sd = apply(fadu_data2n_calc[2:4],1,sd)
fadu_data2n_calc = fadu_data2n_calc %>% mutate(cbd_min = cbd - 2*cbd_sd, cbd_max = cbd + 2*cbd_sd)
fadu_data2n_calc %>% head()

fadu_data2n_calc %>% write_csv("fadu_data2_nuclei_calc.csv")
fr_change2n = matrix(fadu_data2n_calc$fr_change,16,24, byrow = T)
pheatmap(fr_change2n,cluster_rows = F,cluster_cols = F)
# saved as "FADU run 2 nuclei count fr_change.png"

fadu_data2n_calc %>% filter(Well %in% fadu_controls) %>% select(fr_change) %>% min() # 0.3801222
fadu_data2n_calc %>% filter(Well %in% fadu_controls) %>% select(fr_change) %>% max() # 0.634187
fadu_data2n_calc %>% filter(Well %in% fadu_inner_compounds) %>% select(log_fr_change) %>% min() # -1.617
fadu_data2n_calc %>% filter(Well %in% fadu_inner_compounds) %>% select(log_fr_change) %>% max() # 0.035

# https://www.programmingr.com/statistics/z-score-in-r/

fadu_data2n_calc %>% filter(Well %in% fadu_controls) %>% pull(fr_change) %>% shapiro.test()
# NOT normally distributed but looks just fine

fadu_data2n_calc %>% filter(Well %in% fadu_controls) %>% 
          ggplot(aes(x = fr_change)) + geom_density()
ggsave("FADU run2 nuclei count controls fr_change distribution.png")


fadu2n_controls_mean = fadu_data2n_calc %>% filter(Well %in% fadu_controls) %>% pull(fr_change) %>% mean()
fadu2n_controls_sd = fadu_data2n_calc %>% filter(Well %in% fadu_controls) %>% pull(fr_change) %>% sd()
# fadu_data_calc %>% filter(Well %in% fadu_controls) %>% pull(veh) %>% mean()
fadu2n_controls_mean # 0.5182667
fadu2n_controls_sd # .04414386

fadu_hits2n = fadu_data2n_calc %>% filter(Well %in% fadu_inner_compounds,(fr_change > (fadu2n_controls_mean + 2* fadu2n_controls_sd)| fr_change < (fadu2n_controls_mean - 2* fadu2n_controls_sd))) 
fadu_hits2n %>% write_csv("fadu_hits2_nuclei.csv")


# fadu_data %>% filter(Well %in% fadu_controls) %>% select(V_pl_1,V_pl_2) %>% cor() #0.3645661
# fadu_data %>% filter(Well %in% fadu_controls) %>% select(CBD_pl_1,CBD_pl_2) %>% cor() # 0.3401433
# fadu_data %>% filter(Well %in% fadu_inner_compounds) %>% select(V_pl_1,V_pl_2) %>% cor() #0.7407735
# fadu_data %>% filter(Well %in% fadu_inner_compounds) %>% select(CBD_pl_1,CBD_pl_2) %>% cor() # 0.6920442
# fadu_data %>% filter(Well %in% fadu_controls) %>% 
#   ggplot(aes(x = V_pl_1, y = V_pl_2))+
#   geom_point()
# fadu_data %>% filter(Well %in% fadu_inner_compounds) %>% 
#   ggplot(aes(x = V_pl_1, y = V_pl_2))+
#   geom_point()
# 
# fadu_data %>% filter(Well %in% fadu_controls) %>% 
#   ggplot(aes(x = V_pl_1, y = V_pl_2))+
#   geom_point()
# fadu_data %>% filter(Well %in% fadu_inner_compounds) %>% 
#   ggplot(aes(x = V_pl_1, y = CBD_pl_1))+
#   geom_point()

fadu_data2n_calc_for_plot_comp = fadu_data2n_calc %>% filter(Well %in% fadu_inner_compounds) %>% 
          select(veh,cbd,cbd_min,cbd_max,cbd_sd, veh_sd ) %>% mutate(well = "compounds")
fadu_data2n_calc_for_plot_control = fadu_data2n_calc %>% filter(Well %in% fadu_controls) %>% 
          select(veh,cbd,cbd_min,cbd_max,cbd_sd, veh_sd) %>% mutate(well = "controls")
fadu_data2n_calc_for_plot_gem = fadu_data2n_calc %>% filter(Well %in% fadu_gem) %>% 
          select(veh,cbd,cbd_min,cbd_max,cbd_sd, veh_sd) %>% mutate(well = "gem")
fadu_data2n_calc_for_plot_doc = fadu_data2n_calc %>% filter(Well %in% fadu_doc) %>% 
          select(veh,cbd,cbd_min,cbd_max,cbd_sd, veh_sd) %>% mutate(well = "doc")

fadu_data2n_calc_for_plot = bind_rows(fadu_data2n_calc_for_plot_control,
                                     fadu_data2n_calc_for_plot_comp,
                                     fadu_data2n_calc_for_plot_gem,
                                     fadu_data2n_calc_for_plot_doc)
fadu_data2n_calc_for_plot  %>% 
          ggplot(aes(x = veh, y = cbd, color = well))+
          geom_point() + geom_smooth(method = 'lm')

fadu_data2n_calc_for_plot  %>%  
          ggplot(aes(x = veh, y = cbd, color = well))+
          geom_point() + geom_smooth(method = 'lm') +
          geom_errorbar(aes(ymin = cbd_min, ymax = cbd_max), width = 0.2)

fadu_data2n_calc_for_plot  %>% 
          ggplot(aes(x = veh, y = cbd, color = well))+
          geom_point() + geom_smooth(method = 'lm') +
          geom_errorbar(aes(ymin = cbd - 2*cbd_sd, ymax = cbd + 2*cbd_sd), width = 0.2)
# same result
ggsave("veh vs cbd run2 nuclei detailed_2.png")

fadu_data2n_calc_for_plot  %>% 
          ggplot(aes(x = veh, y = cbd, color = well))+
          geom_point() + geom_smooth(method = 'lm') +
          geom_errorbar(aes(ymin = cbd - 2*cbd_sd, ymax = cbd + 2*cbd_sd), width = 0.2) +
          geom_errorbarh(aes(xmin = veh - 2*veh_sd, xmax = veh + 2*veh_sd), height = 0.2) +
          coord_cartesian(xlim = c(2000,16000))
ggsave("veh vs cbd run2 nuclei detailed_3.png")

         
# -------------------


ggsave("veh vs cbd run2 nuclei detailed_2.png")

fadu_hits2_comb_wells = append(fadu_hits2$Well,fadu_hits2n$Well) %>% unique()
fadu_hits2_extr = fadu_data2_calc %>% filter(Well %in% fadu_hits2_comb_wells) %>% 
          select(Well, fr_change) %>% rename(fr_change_CTB = fr_change)
fadu_hits2n_extr = fadu_data2n_calc %>% filter(Well %in% fadu_hits2_comb_wells) %>% 
          select(Well, fr_change) %>% rename(fr_change_Nuc = fr_change)
fadu_hits2_comparison = left_join(fadu_hits2_extr,fadu_hits2n_extr)
pheatmap(as.matrix(select(fadu_hits2_comparison, -Well)),cluster_rows = F,cluster_cols = F)
fadu_hits2_comparison %>% write_csv("fadu_hits2_comparison.csv")
