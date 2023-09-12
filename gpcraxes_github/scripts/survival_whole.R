# Custom addition function to handle character concatenation
`+` <- function(e1, e2) {
  if (is.character(e1) | is.character(e2)) {
    paste0(e1, e2)
  } else {
    base::`+`(e1, e2)
  }
}

# Custom function to calculate maximum with handling NA values
my.max <- function(x) ifelse(!all(is.na(x)), max(x, na.rm = TRUE), NA)

# Load required libraries
library(survival)
library(dplyr)

# Load receptor-ligand OR receptor-enzyme pairs CSV file

#receptor_enzyme file
#rl_list <- read.csv('/data/survival/RE_pairs.csv')

#receptor_ligand file
rl_list <- read.csv('/data/survival/RL_pairs.csv')

# Extract receptor and ligand columns
recep <- rl_list$Recep
ligand <- rl_list$Ligand

# List of cancer types
can=c('Bile duct',
     'Adrenal gland',
     'White blood cell',
     'Head and Neck region',
      'Cervix',
     'Lining of body cavities',
     'Paraganglia',
     'Eye',
     'Endometrium',
     'Rectum',
     'Lymphatic tissue',
     'Thymus',
     'Skin',
     'Prostate',
     'Testis',
     'Lung',
     'Breast',
     'Kidney',
     'Colon',
     'Esophagus',
     'Stomach',
     'Bladder',
     'Uterus',
     'Brain',
     'thyroid',
     'Pancreas',
     'Liver',
     'Ovary')

# Loop through each cancer type
for (cancer in can) {
  
  # Read data CSV file for the current cancer type
  data1 <- read.csv('/data/survival/data/' + cancer + '.csv')
  
  # Remove rows with missing vital status
  data1 <- subset(data1, vital_status != '[NA]')
  
  # Define columns for the output data frame
  x <- c("Receptor-Ligand", "Receptor_ref", "HR", "logrank-p", "%95I", "Ligand_ref", "HR", "logrank-p", "%95I", "Receptor-Ligand_ref", "HR", "logrank-p", "%95I", "N1", "N2", "C")
  d <- setNames(data.frame(matrix(ncol = 16, nrow = 0)), x)
  
  # Loop through each receptor-ligand pair
  for(i in seq(from=1, to=length(recep), by=1))
    #
  {  
    tryCatch({
    data2=data1
    r=data2 %>% select(recep[i])
    l=data2 %>% select(ligand[i])
    

# ... Data preprocessing ...

    data2 <- mutate(data2, r= ifelse(select(data2, c(recep[i]))<=median(as.numeric(unlist(select(data2, c(recep[i]))))), "Low", "High"))
    data2 <- mutate(data2, l = ifelse(select(data2, c(ligand[i]))<=median(as.numeric(unlist(select(data2, c(ligand[i]))))), "Low", "High"))
    
    
    
    data2$r=factor(data2$r)
    data2$l=factor(data2$l)
    
 
# Perform Cox proportional hazards survival analysis
   
    #RECEPTOR
    surv_object <- Surv(time = data2$OS.time/(30*12), event = as.numeric(data2$vital_status))
    fit1 <- survfit(surv_object~r,data=data2)
    fit1.coxph <- coxph(surv_object~r,data=data2)
    a<-summary(fit1.coxph)
    b<-coef(a)
    
    r_ref=paste(rownames(b))
    r_hr=b[2]
    r_p=b[5]
    r_ci=paste(a$conf.int[3],a$conf.int[4])
    
    
    
    #LIGAND
    surv_object <- Surv(time = data2$OS.time/(30*12), event = as.numeric(data2$vital_status))
    fit1 <- survfit(surv_object~l,data=data2)
    fit1.coxph <- coxph(surv_object~l,data=data2)
    a<-summary(fit1.coxph)
    b<-coef(a)
    
    l_ref=paste(rownames(b),collapse ='|')
    l_hr=b[2]
    l_p=b[5]
    l_ci=paste(a$conf.int[3],a$conf.int[4])
    
    #table(data2$l)
    
    
    
    #RECEPTOR-LIGAND HH vs LL
    data2<-subset(data2,((r=='High')&(l=='High'))|((r=='Low')&(l=='Low')))
    
    surv_object <- Surv(time = data2$OS.time/(30*12), event = as.numeric(data2$vital_status))
    fit1 <- survfit(surv_object~r+l,data=data2)
    fit1.coxph <- coxph(surv_object~r+l,data=data2)
    a<-summary(fit1.coxph)
    b<-coef(a)
    
    rl_ref=paste(rownames(b),collapse ='|')
    rl_hr=b[1,2]
    rl_p=b[1,5]
    rl_ci=paste(a$conf.int[1,3],a$conf.int[1,4])
    n1=fit1$n[1]
    n2=fit1$n[2]
    #rl_w_p=a$waldtest[3]
    rl_c=a$concordance


    
      # Append the row to the output data frame
      d[nrow(d) + 1, ] = c(recep[i] + '-' + ligand[i], r_ref, r_hr, r_p, r_ci, l_ref, l_hr, l_p, l_ci, rl_ref, rl_hr, rl_p, rl_ci, n1, n2, rl_c)
      
      # Remove temporary data
      rm(data2)
    }, error = function(e) {})
  }
  
  # Write the output data frame to a CSV file
  write.csv(d, '/data/survival/output/' + cancer + '.csv', row.names = FALSE)
  
  # Remove temporary data
  rm(data1)

}

