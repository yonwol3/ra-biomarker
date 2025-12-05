library(corrplot)
library(tidyverse)

# Calculate the pearson correlation among sample A biomarkers
# source("clean-data-A.R")
biomarkers_A<-c("RF IgA","RF IgM","RF IgG","ACPA IgA","ACPA IgM","ACPA IgG")
data_A<- data_A %>% 
  rename_with(~ biomarkers_A, .cols = c("igarfconc_","igmrfconc_" , "iggrfconc_"
                                      ,"igaccpavgconc", "igmccpavgconc", "iggccpavgconc"))
# Pearson correlation matrix
cor_matrix_A<- round (cor(data_A[ ,biomarkers_A]), 2)
png(file = "figures/Sample A correlation.png", width = 10,height = 10,units = "in",     
    res = 300)      
corrplot( cor_matrix_A, method="color")
dev.off()
# Calculate the pearson correlation among sample B biomarkers
# source ("clean-data-B.R")
old_biomarker_labels<-c("aptivaccp3igg_≥5#00flu", "aptiva_acpafsiggvimentin2_≥5#00au",
                       "aptiva_acpafsiggfibrinogen_≥5#00au","aptiva_acpafsigghistone1_≥5#00au",
                       "aptivaccp3iga_≥5#00flu", "aptiva_acpafsigavimentin2_≥5#00au",
                       "aptiva_acpafsigafibrinogen_≥5#00au","aptiva_acpafsigahistone1_≥5#00au")
biomarkers_B<- c("anti-CCP3 (IgG)", "anti-citVim2 (IgG)", "anti-citFib (IgG)",
                 "anti-citHis1 (IgG)", "anti-CCP3 (IgA)", "anti-citVim2 (IgA)",
                 "anti-citFib (IgA)", "anti-citHis1 (IgA)")

data_B<- data_B %>% 
  rename_with(~ biomarkers_B, .cols = all_of(old_biomarker_labels))

cor_matrix_B<- round (cor(data_B[ ,biomarkers_B]), 2)
png(file = "figures/Sample B correlation.png", width = 10,height = 10,units = "in",     
    res = 300)      
corrplot( cor_matrix_B, method="color")
dev.off()

