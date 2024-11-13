################################
## code to create table-1 
# By: Yonatan Wolde
###############################


library(tidyverse)
library(labelled)
library(table1)
library(naniar)

source("clean-data.R")


dat<- clean

dat<-dat%>%
  dplyr::rename('RF IgA'= igarfconc_ , 
                'RF IgM'= igmrfconc_ , 
                'RF IgG'= iggrfconc_ , 
                'ACPA IgA'= igaccpavgconc, 
                'ACPA IgM'= igmccpavgconc, 
                'ACPA IgG'= iggccpavgconc)


dat$fem<-factor(as.character(dat$fem), levels=c(0,1), labels=c("M", "F")) # sex 
dat$nw<-factor(as.character(dat$nw), levels = c(0,1), labels=c("White", "Non-White")) # nowhite status
dat$famhx<- factor(as.character(dat$famhx), levels = c(0,1), labels=c("No", "Yes")) #family RA history
dat$diagnosis<- factor(dat$diagnosis, levels = c(0,1), labels=c("No-RA", "RA"))
dat <- set_variable_labels(dat,
                            fem = "Sex",
                            nw= "Race",
                            famhx = "Family RA History")


# Then create your Table 1
covariates <- c("fem", "nw")
biomarkers <- c("RF IgA", "RF IgM", "RF IgG", "ACPA IgA","ACPA IgM","ACPA IgG" )
t1_variables<-c(covariates,biomarkers)
categorical_vars <-covariates


# Enclose variable names with backticks to handle spaces
t1_variables_backticked <- paste0("`", t1_variables, "`")

# Create the table with the corrected formula
t1 <- table1(
  as.formula(paste0("~", paste0(t1_variables_backticked, collapse = "+"), "|", "diagnosis")),
  data = dat,
  caption = "Table1: Descriptive summary of data by RA status"
)

t1  # Display the table







