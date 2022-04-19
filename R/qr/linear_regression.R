require(data.table)


## Diabetes data
diabetes<-fread('/home/anne/Documents/featurecloud/pca/qr/data/diabetes/diabetes.tsv')
diabetes<-diabetes[0:(nrow(diabetes)-20),]

lables <-fread('/home/anne/Documents/featurecloud/pca/qr/data/diabetes/diabetes_lab.tsv')
lables<-lables[0:(nrow(diabetes))]$V1

f<-formula(lables~(-1)+.)
model<-lm(f, data = diabetes)
su<-summary(model)
coeff.diabetes<-model$coefficients
rsquared.diabetes<-su$r.squared
pvals.diabetes<-su$coefficients[,4]

pyres<-fread('/home/anne/Documents/featurecloud/pca/qr/results/diabetes.txt', fill=T)

sum(pyres[1, 2:ncol(pyres)]-coeff.diabetes)
sum(pyres[2,2]-rsquared.diabetes)
sum(pyres[3,2:ncol(pyres)]-pvals.diabetes)




### Fish data
fish<-fread('/home/anne/Documents/featurecloud/pca/qr/data/fish/Fish.csv')
weight<-fish[,2]
fish<-fish[,3:ncol(fish)]
ffish<-formula(weight$Weight~(-1)+.)
model.fish<-lm(ffish, data = fish)
su.fish<-summary(model.fish)
coeff.fish <-model.fish$coefficients
rsquared.fish<-su.fish$r.squared
pvals.fish<-su.fish$coefficients[,4]

pyres<-fread('/home/anne/Documents/featurecloud/pca/qr/results/fish.txt', fill = T)

sum(pyres[1, 2:ncol(pyres)]-coeff.fish)
pyres[2,2]-rsquared.fish
sum(pyres[3,2:ncol(pyres)]-pvals.fish)



  ### 
  who<-fread('/home/anne/Documents/featurecloud/pca/qr/data/who/who2.csv')
  who<-who[Year==2012 & Status == 'Developing']
  colnames(who)<-sapply(colnames(who), function(x) gsub(' ', '_',x))
  le<-who$Life_expectancy
  who<-who[, c(5:ncol(who)), with=F]
  who.formula<-formula(le~(-1)+.)
  model.who<-lm(who.formula, data = who)
  su.who<-summary(model.who)
  coeff.who <-model.who$coefficients
  rsquared.who<-su.who$r.squared
  pvals.who<-su.who$coefficients[,4]
  
  pyres<-fread('/home/anne/Documents/featurecloud/pca/qr/results/who.txt', fill=T)
  
  sum(pyres[1, 2:ncol(pyres)]-coeff.who)
  pyres[2,2]-rsquared.who
  abs(sum(pyres[3,2:ncol(pyres)]-pvals.who))

