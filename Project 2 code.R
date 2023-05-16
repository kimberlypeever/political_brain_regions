#Loading in dataset
brain_responses=read.csv("n90_pol.csv")

#Exploratory data analysis
#Univariate analysis
hist(brain_responses$amygdala, main="Histogram of Amygdala", xlab="Amygdala")
summary(brain_responses$amygdala)
hist(brain_responses$acc, main="Histogram of Acc", xlab="Acc")
summary(brain_responses$acc)
hist(brain_responses$orientation, breaks = 5:0, labels=TRUE, 
     main="Histogram of Orientation",xlab="Orientation")
summary(brain_responses$orientation)

#The correlation between volumes of amygdala and acc
plot(brain_responses$amygdala, brain_responses$acc, main="Amygdala vs. Acc"
     ,xlab="Amygdala", ylab="Acc")
abline(0,0)
cor(brain_responses$amygdala, brain_responses$acc)

#The correlation between orientation and the volume of the amygdala
palette=colorRampPalette(c("#FF3030", "#0000FF"))(5)
plot(brain_responses$amygdala, brain_responses$orientation,
     main="Amygdala vs. Orientation",xlab="Amygdala", ylab="Orientation",
     col=palette[brain_responses$orientation])
cor(brain_responses$amygdala, brain_responses$orientation,method="s")

#The correlation between orientation and the volume of the acc. 
plot(brain_responses$acc, brain_responses$orientation,
     main="Acc vs. Orientation",xlab="Acc", ylab="Orientation", 
     col=palette[brain_responses$orientation])
cor(brain_responses$acc, brain_responses$orientation, method="s")

#95% bootstrap confidence intervals for above correlations. 
#For amygdala
nboot=500
boot_vec=rep(NA,nboot)
set.seed(444)
for (b in 1:nboot) {
  boot_data=brain_responses[sample(1:90,replace = T),] 
  boot_cor=cor(boot_data$amygdala, boot_data$orientation,method = "s")
  boot_vec[b]=boot_cor
}
quantile(boot_vec,probs=c(0.025,0.975))

#For acc
boot_vec=rep(NA,nboot)
for (b in 1:nboot) {
  boot_data=brain_responses[sample(1:90,replace = T),] 
  boot_cor=cor(boot_data$acc, boot_data$orientation,method = "s")
  boot_vec[b]=boot_cor
}
quantile(boot_vec,probs=c(0.025,0.975))

#Linear regression model for orientation on the volumes of amygdala and acc
brain_lm=lm(orientation~amygdala+acc, data=brain_responses)
summary(brain_lm)

#Main analysis
#creating a binary response variable 'conservative', which is 1 when the 
#student has orientation less than or equal to 2, and 0 otherwise.
brain_responses$conservative=as.numeric(brain_responses$orientation<=2)
mean(brain_responses$conservative)*90

#Fitting a logistic regression of conservative on the linear effects of the 
#volumes of amygdala and acc.
brain_responses_glm=glm(conservative~amygdala+acc, data = brain_responses, 
                        family = binomial)
summary(brain_responses_glm)
#Performing model calibration for logistic model
brain_responses$pred=predict(brain_responses_glm,type="response")
brain_responses$pred_bin=cut(brain_responses$pred,breaks = seq(0,1,0.05))
cali_table=aggregate(brain_responses$conservative,
                     by=list(bin=brain_responses$pred_bin)
                     ,FUN=function(x)c(mean=mean(x),count=length(x)))
cali_table$pred=c(0.025, 0.075, 0.125, 0.175, 0.225, 0.325, 0.375, 0.425,
                  0.525, 0.575, 0.625)
cali_table$se= sqrt(cali_table$pred * 
                      (1 - cali_table$pred)/cali_table$x[,"count"])
#Plotting calibration plot for logistic model
plot(cali_table$x[,"mean"]~cali_table$pred, 
     main="Predicted Probability vs. Actual Frequency for Logistic Model",
     xlab="Predicted Probability", ylab="Actual Frequency", 
     xlim = c(0, 0.7), ylim = c(0, 1))
abline(0, 1, col = "grey")
rug(fitted(brain_responses_glm), col = "grey")
segments(x0 = cali_table$pred, y0 = cali_table$pred - 1.96 * cali_table$se, 
         y1 = cali_table$pred +  1.96 * cali_table$se)

#Finding in-sample mis-classification rate using 0.5 cutoff for logistic
#model
cons_pred_prob=data.frame(predict(brain_responses_glm,type="response"))
brain_responses_merge=merge(x=brain_responses,y=cons_pred_prob,all=TRUE,
                            by="row.names")
brain_responses_merge$cons_pred=
  as.numeric(brain_responses_merge$
               predict.brain_responses_glm..type....response..>0.5)
table(brain_responses_merge$conservative,brain_responses_merge$cons_pred)
(3+12)/90

#Recalculating the classification error rate for logistic model
#using cross-validation 
set.seed(444)
folds=matrix(sample(nrow(brain_responses)),ncol=5)
cv_error=function(f_fit,form,...){
  cv_cer=rep(NA,5)
  for (i in 1:5) {
    mod=f_fit(form, data = brain_responses[c(folds[,-i]),], 
              family = binomial,...)
    pred=as.numeric(predict(mod,brain_responses,type="response")
                    [folds[,i]]>0.5)
    cons=brain_responses$conservative[folds[,i]]
    cv_cer[i]=mean(pred!=cons)
  }
  return(mean(cv_cer))
}
cv_error(glm,conservative~amygdala+acc)                        


#Fitting a generalized additive model for conservative on the volumes of 
#amygdala and acc.
library(mgcv)
brain_responses_gam=gam(conservative~s(amygdala)+s(acc), 
                        data=brain_responses, family=binomial)
summary(brain_responses_gam)
plot(brain_responses_gam,select=1, main="Estimated Effect of Amygdala",
     ylab="Estimate", xlab="Amygdala")
plot(brain_responses_gam,select=2, main="Estimated Effect of Acc",
     ylab="Estimate", xlab="Acc")

#Performing model calibration for GAM model
brain_responses$pred=predict(brain_responses_gam,type="response")
brain_responses$pred_bin=cut(brain_responses$pred,breaks = seq(0,1,0.05))
cali_table=aggregate(brain_responses$conservative, 
                     by=list(bin=brain_responses$pred_bin)
                     ,FUN=function(x)c(mean=mean(x),count=length(x)))
cali_table$pred=c(0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 
                  0.425, 0.475, 0.525, 0.625, 0.675, 0.775)
cali_table$se= sqrt(cali_table$pred *
                      (1 - cali_table$pred)/cali_table$x[,"count"])

#Plotting calibration plot for Gam model.
plot(cali_table$x[,"mean"]~cali_table$pred, 
     main="Predicted Probability vs. Actual Frequency for GAM Model",
     xlab="Predicted Probability", ylab="Actual Frequency",
     xlim = c(0, 0.8), ylim = c(0, 1))
abline(0, 1, col = "grey")
rug(fitted(brain_responses_gam), col = "grey")
segments(x0 = cali_table$pred, y0 = cali_table$pred - 1.96 
         * cali_table$se, y1 = cali_table$pred +  1.96 * cali_table$se)


#Finding in-sample mis-classification rate using 0.5 cutoff for GAM model
cons_pred_prob2=data.frame(predict(brain_responses_gam,type="response"))
brain_responses_merge2=merge(x=brain_responses,y=cons_pred_prob2,all=TRUE,
                             by="row.names")
brain_responses_merge2$cons_pred2=
  as.numeric(brain_responses_merge2$
               predict.brain_responses_gam..type....response..>0.5)
table(brain_responses_merge2$conservative,brain_responses_merge2$cons_pred2)
(10+3)/90

#Recalculating the classification error rate for GAM model using 
#cross-validation 
cv_error(gam, conservative~s(amygdala)+s(acc) )

#Calculating deviance
set.seed(444)
simulate.from.logr <- function(df, mdl) {
  probs <- predict(mdl, newdata = df, type = "response")
  df$conservative <- rbinom(n = nrow(df), size = 1, prob = probs)
  return(df)
}
delta.deviance.sim <- function(df, mdl) {
  sim.df <- simulate.from.logr(df, mdl)
  GLM.dev <- glm(conservative~amygdala+acc, data = sim.df, 
                  family = "binomial")$deviance
  GAM.dev <- gam(conservative~s(amygdala)+s(acc), data=sim.df,
                 family='binomial')$deviance
  return(GLM.dev - GAM.dev)
}
(delta.dev.observed <- brain_responses_glm$deviance - 
    brain_responses_gam$deviance)
delta.dev <- replicate(100, delta.deviance.sim(brain_responses, 
                                               brain_responses_glm))
mean(delta.dev.observed <= delta.dev)