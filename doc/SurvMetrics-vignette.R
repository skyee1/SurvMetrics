## ----fig.align='center', message=FALSE, warning=FALSE, results='hide'---------
library(SurvMetrics)
library(caret)
library(randomForestSRC)
library(survival)  
library(pec)
library(ggplot2)

#Initialization
metrics_cox = 0
metrics_rsf = 0
for (i in 1:20) {
  
  mydata = SDGM4(N = 200, p = 20, c_step = 0.5)
  index_data = createFolds(1:nrow(mydata), 2)
  train_data = mydata[index_data[[1]],]
  test_data = mydata[index_data[[2]],]
  
  #fit the models
  fitrsf = rfsrc(Surv(time, status)~., data = train_data, nsplit = 3, ntree = 500)
  mat_rsf = predict(fitrsf, test_data)$survival
  
  dis_time = fitrsf$time.interest
  fitcox = coxph(Surv(time, status)~., data = train_data, x = TRUE)
  mat_cox = predictSurvProb(fitcox, test_data, dis_time)
  
  #calculate the C index
  med_index = median(1:length(dis_time))
  surv_obj = Surv(test_data$time, test_data$status)
  
  #C index for Cox
  metrics_cox[i] = Cindex(surv_obj, predicted = mat_cox[, med_index])
  #C index for RSF
  metrics_rsf[i] = Cindex(surv_obj, predicted = mat_rsf[, med_index])
  
}

data_CI = data.frame('Cindex' = c(metrics_cox, metrics_rsf),
                     'model' = c(rep('Cox', 20), rep('RSF', 20)))

ggplot(data_CI, aes(x = model, y = Cindex, fill = model)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FFBBCC", "#88CCFF"))


## ----fig.align='center', warning=FALSE, results='hide'------------------------
#Initialization
metrics_cox = 0
metrics_rsf = 0
for (i in 1:20) {
  
  mydata = SDGM4(N = 200, p = 20, c_step = 0.5)
  index_data = createFolds(1:nrow(mydata), 2)
  train_data = mydata[index_data[[1]],]
  test_data = mydata[index_data[[2]],]
  
  #fit the models
  fitrsf = rfsrc(Surv(time, status)~., data = train_data, nsplit = 3, ntree = 500)
  mat_rsf = predict(fitrsf, test_data)$survival
  
  dis_time = fitrsf$time.interest
  fitcox = coxph(Surv(time, status)~., data = train_data, x = TRUE)
  mat_cox = predictSurvProb(fitcox, test_data, dis_time)
  
  #calculate the C index
  med_index = median(1:length(dis_time))
  surv_obj = Surv(test_data$time, test_data$status)
  t_star = median(fitrsf$time.interest)
  
  #Brier Score for Cox
  metrics_cox[i] = Brier(surv_obj, pre_sp = mat_cox[, med_index], t_star)
  #Brier Score for RSF
  metrics_rsf[i] = Brier(surv_obj, pre_sp = mat_rsf[, med_index], t_star)
  
}

data_BS = data.frame('BS' = c(metrics_cox, metrics_rsf),
                     'model' = c(rep('Cox', 20), rep('RSF', 20)))

ggplot(data_BS, aes(x = model, y = BS, fill = model)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FFBBCC", "#88CCFF"))

## ----fig.align='center', message=FALSE, warning=FALSE, results='hide'---------
#Initialization
metrics_cox = 0
metrics_rsf = 0
for (i in 1:20) {
  
  mydata = SDGM4(N = 100, p = 20, c_step = 0.5)
  index_data = createFolds(1:nrow(mydata), 2)
  train_data = mydata[index_data[[1]],]
  test_data = mydata[index_data[[2]],]
  
  #fit the models
  fitrsf = rfsrc(Surv(time, status)~., data = train_data, nsplit = 3, ntree = 500)
  mat_rsf = predict(fitrsf, test_data)$survival
  
  dis_time = fitrsf$time.interest
  fitcox = coxph(Surv(time, status)~., data = train_data, x = TRUE)
  mat_cox = predictSurvProb(fitcox, test_data, dis_time)
  
  #calculate the C index
  med_index = median(1:length(dis_time))
  surv_obj = Surv(test_data$time, test_data$status)

  
  #IBS for Cox
  metrics_cox[i] = IBS(surv_obj, sp_matrix = mat_cox, dis_time)
  #IBS for RSF
  metrics_rsf[i] = IBS(surv_obj, sp_matrix = mat_rsf, dis_time)
  
}

data_IBS = data.frame('IBS' = c(metrics_cox, metrics_rsf),
                     'model' = c(rep('Cox', 20), rep('RSF', 20)))

ggplot(data_IBS, aes(x = model, y = IBS, fill = model)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FFBBCC", "#88CCFF"))

## ----fig.align='center', message=FALSE, warning=FALSE, results='hide'---------
#Initialization
metrics_cox_IAE = 0
metrics_rsf_IAE = 0
metrics_cox_ISE = 0
metrics_rsf_ISE = 0
for (i in 1:20) {
  
  mydata = SDGM4(N = 100, p = 20, c_step = 0.5)
  index_data = createFolds(1:nrow(mydata), 2)
  train_data = mydata[index_data[[1]],]
  test_data = mydata[index_data[[2]],]
  
  #fit the models
  fitrsf = rfsrc(Surv(time, status)~., data = train_data, nsplit = 3, ntree = 500)
  mat_rsf = predict(fitrsf, test_data)$survival
  
  dis_time = fitrsf$time.interest
  fitcox = coxph(Surv(time, status)~., data = train_data, x = TRUE)
  mat_cox = predictSurvProb(fitcox, test_data, dis_time)
  
  #calculate the C index
  med_index = median(1:length(dis_time))
  surv_obj = Surv(test_data$time, test_data$status)
  
  
  #C index for Cox
  temp1 = IAEISE(surv_obj, sp_matrix = mat_cox, dis_time)
  metrics_cox_IAE[i] = temp1[1]
  metrics_cox_ISE[i] = temp1[2]
  #C index for RSF
  temp2 = IAEISE(surv_obj, sp_matrix = mat_rsf, dis_time)
  metrics_rsf_IAE[i] = temp2[1]
  metrics_rsf_ISE[i] = temp2[2]
  
}

data_IAE = data.frame('IAE' = c(metrics_cox_IAE, metrics_rsf_IAE),
                     'model' = c(rep('Cox', 20), rep('RSF', 20)))

data_ISE = data.frame('ISE' = c(metrics_cox_ISE, metrics_rsf_ISE),
                     'model' = c(rep('Cox', 20), rep('RSF', 20)))

P1 = ggplot(data_IAE, aes(x = model, y = IAE, fill = model)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FFBBCC", "#88CCFF")) +
  theme(legend.position = 'none')

P2 = ggplot(data_ISE, aes(x = model, y = ISE, fill = model)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FFBBCC", "#88CCFF")) +
  theme(legend.position = 'none')

library(ggpubr)
ggarrange(P1, P2, ncol = 2)



