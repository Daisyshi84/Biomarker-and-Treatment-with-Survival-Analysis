rm(list = ls())
library(tidyverse); library(survival); library(survminer); library(plotly)
library(ggplot2); library(psych); library(biostat3); library(lubridate);library(data.table)
library(kableExtra);library(gtsummary);
memory.limit(90000000)

###########################################################
# 1. Explore the data and report on any potential data issues
##########################################################
setwd("C:/Users/SSE6/Desktop/CoxPH Model")
data<- readRDS("VeracyteToyData.rds")
summary(data)
str(data)
# list of variable names to convert to numeric
vars_to_convert <- c("CLIN_COV2", "STATUS_1", "STATUS_2","STATUS_3","STATUS_4","STATUS_5")
# for loop to convert variables to numeric
for (var in vars_to_convert) {
  data <- data %>% 
    mutate(!!sym(var) := as.numeric(!!sym(var)))
}

# list of date variable to convert to date format
df_dates <-  data %>% select(starts_with("VISIT"),BASELINE_VISIT)
# loop through each column in the original data frame
for (col in names(df_dates)) {
  # convert the column to a date format
  data[[col]] <-  mdy(data[[col]])
}


str(data)

summary(data$CLIN_COV2)
summary(data$MARKER)
 
# Check for duplicated observations
sum(duplicated(data))

# Check for unique IDs
length(unique(data$ID))

data %>% group_by(ID) %>% 
  count(ID) %>% filter(n>1)

 table(data$STATUS_1)
 table(data$STATUS_2)
 table(data$STATUS_3)
 table(data$STATUS_4)
 table(data$STATUS_5)

dup<- data %>% filter(ID %in% c("P125","P15"))  

cat("We noticed that ID number:",unique(dup$ID),
    "has more than 1 record, the fact that participant P15 has 2 records with different baseline dates and follow-up visit data not being the same raises a data quality issue.
This could be due to errors in data entry or data management. participant P125 has 3 records, it could be a real occurrence of the participant having the event of interest twice with a long time gap between the two events.
To address this issue, the investigator should review the data carefully and investigate whether the two records for P15 are accurate. They could also contact the study site to obtain additional information about the participant's medical history to confirm whether they had experienced the event of interest twice for participant P125.
If the two records are accurate, the investigator should decide whether to exclude one of the records.")

#missingness checking
missing_check<-data %>%
  summarize(across(all_of(names(data)),list(".Missing_n"=function(x) sum(is.na(x)),
                                            ".Reported_n"=function(x) sum(is.na(x)==F),
                                            ".Percent_missing"=function(x) round(mean(is.na(x)*100),1))))

# reformat table from wide to long
missing_check.tidy <- missing_check %>% 
  gather(stat, val) %>%  
  separate(stat, into = c("Variable", "stat"), sep = "[.]") %>%  
  mutate(Variable=substr(Variable,1,nchar(Variable)-1)) %>%    
  spread(stat,val) 

### Check 1 - missingness

# output nicely with kable package
kable(missing_check.tidy,format="html" ,booktabs=T,align="lrrr",toprule = "\\toprule") %>%
  row_spec(as.numeric(row.names(missing_check.tidy )[is.na(missing_check.tidy$Percent_missing)==F & missing_check.tidy$Missing_n==0]),background="#a1d99b") %>% # color row dark green if none missing
  row_spec(as.numeric(row.names(missing_check.tidy )[is.na(missing_check.tidy$Percent_missing)==F & missing_check.tidy$Missing_n>0 & missing_check.tidy$Percent_missing<=1]),background="#c7e9c0") %>% # color row light green if =<1% missing (can change this to whatever threshold)
  row_spec(as.numeric(row.names(missing_check.tidy )[is.na(missing_check.tidy$Percent_missing)==F & missing_check.tidy$Percent_missing>1 & missing_check.tidy$Reported_n!=0]),background="#ffeda0") %>% # color row yellow if >1% and <100% missing (can change this to whatever threshold)
  row_spec(as.numeric(row.names(missing_check.tidy )[is.na(missing_check.tidy$Percent_missing)==F & missing_check.tidy$Reported_n==0]),background="#fc9272") %>% # color row red if all missing?
  column_spec(1:4,width="1.5in") %>%
  row_spec(0, align = "c")



### Check 2 -   we notice STATUS_2,STATUS_4 reported number not coupled with the corresponding VISIT_2,VISIT_4
data_ex <- data %>%
  mutate(not_coupled_2 = ifelse(
    is.na(STATUS_2)==F & is.na(VISIT_2)==F | 
      ifelse(is.na(STATUS_2)==T & is.na(VISIT_2)==T, TRUE, FALSE),
    "OK",
    "Issue"
  ),
  not_coupled_4 = ifelse(
    is.na(STATUS_4)==F & is.na(VISIT_4)==F | 
      ifelse(is.na(STATUS_4)==T & is.na(VISIT_4)==T, TRUE, FALSE),
    "OK",
    "Issue"
  ))

issue<- data_ex %>% filter(not_coupled_2=="Issue"|not_coupled_4=="Issue")
cat("we also noticed 2 patients ID number:", issue$ID,  ",their STATUS not coupled with the corresponding VISIT.")


###Check 3 - CLIN_COV2 has value 2500 outlier, so we replace using mean value;

data <- data %>%
  mutate(CLIN_COV2 = if_else(ID == "P43", mean(data$CLIN_COV2), CLIN_COV2)) %>% 
  rowwise() %>%
  mutate(LAST_DATE = max(ymd(c_across(starts_with('VISIT'))), na.rm = TRUE),
         DAYS =  as.duration(BASELINE_VISIT %--% LAST_DATE) / ddays(1),
         STATUS = ifelse(STATUS_1 == 1|STATUS_2 == 1|STATUS_3 == 1|STATUS_4 == 1|STATUS_5 == 1,1,0))

# Create some plots
ggplot(data, aes(x=MARKER,y=STATUS, fill=STATUS)) +
  geom_point()+
  ggtitle("Biomarker by Status")
ggplot(data, aes(x=MARKER,y=TREATMENT, fill=STATUS)) +
  geom_point()+
  ggtitle("Biomarker by Treatment")
ggplot(data, aes(x=CLIN_COV2, fill=STATUS)) +
  geom_histogram(bins =30 ) +
  ggtitle("Distribution for CLIN_COV2")
ggplot(data, aes(y=MARKER,x=DAYS,fill=STATUS)) +
  geom_point()+
  ggtitle("Biomarker")

plot_ly(data, x = ~DAYS, y = ~CLIN_COV2, color = ~TREATMENT, 
        colors = c('#636EFA', '#EF553B'), type = "scatter", mode = "markers")
plot_ly(data, x = ~DAYS, y = ~MARKER, color = ~TREATMENT, 
        colors = c('#636EFA', '#EF553B'), type = "scatter", mode = "markers")
plot_ly(data, x = ~DAYS, y = ~STATUS, color = ~TREATMENT, 
        colors = c('#636EFA', '#EF553B'), type = "scatter", mode = "markers")
### descriptive stats
describeBy(data$DAYS, data$TREATMENT)
describeBy(data$DAYS, data$CLIN_COV1)

############################################################################################
# 2. Provide a basic "Table 1" describing the distribution of the biomarker and clinical 
# variables (CLIN_COV1, CLIN_COV2) by treatment group
############################################################################################        

#Table 1: Distribution of Biomarker and Clinical Variables by Treatment Group
table1<- data %>%
  select(TREATMENT,MARKER, CLIN_COV1, CLIN_COV2) %>%
  tbl_summary(by = TREATMENT,
              statistic = list(all_continuous() ~ "{mean} ({sd})")) %>%
  add_p() %>%
  add_overall() %>% 
  bold_labels()

 

############################################################################################
# 3. Fit and briefly interpret Cox PH models to address the investigators hypotheses. Provide 
# a rough visualization to aid in the interpretation of all models.
############################################################################################ 


# Fit Cox PH model  
model1 <- coxph(Surv(DAYS, STATUS) ~ MARKER, data = data)
summary(model1)

model2 <- coxph(Surv(DAYS, STATUS) ~ MARKER+TREATMENT, data = data)
summary(model2)

model3 <- coxph(Surv(DAYS, STATUS) ~ MARKER*TREATMENT, data = data)
summary(model3)


model4 <- coxph(Surv(DAYS, STATUS) ~ MARKER*TREATMENT + CLIN_COV2, data = data)
summary(model4)

model5 <- coxph(Surv(DAYS, STATUS) ~ MARKER*TREATMENT + CLIN_COV1, data = data)
summary(model5)

model6 <- coxph(Surv(DAYS, STATUS) ~ MARKER*TREATMENT + CLIN_COV1 + CLIN_COV2, data = data)
summary(model6)


data$rediduals<- residuals(model3)
plot(x=data$CLIN_COV2,y=data$rediduals)
plot(x=data$MARKER,y=data$rediduals)


ggforest(model3, data = data)
coxph(Surv(DAYS, STATUS) ~ MARKER*TREATMENT, data = data) %>% 
  tbl_regression(exp = TRUE) 

ggsurvplot(survfit(model3), data = data, conf.int = T)


ggadjustedcurves(model3, data = data, method= "conditional")
test_ph_m<- cox.zph(model2)
plot(test_ph_m) 



############################################################################################
# 4. Adjusting for multiple comparisons:
# My suggestion is When conducting multiple tests, it is important to adjust for the increased risk of type 1 errors (false positives) due to chance alone. One common method for controlling the family-wise error rate (FWER) is the Bonferroni correction.

#5. Yes. if CLIN_COV2 is a potential confounder or effect modifier, it should be included in the analysis to improve the accuracy of the estimates, and it should be added as a covariate or interaction term in the Cox proportional hazards model accordingly.
############################################################################################ 
                