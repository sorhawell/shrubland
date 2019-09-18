rm(list=ls())
setwd("~/cpp/shrubland/R")
library(randomForest)
set.seed(1)
makeXY = function(N=500) {
  out = list()
  out$X = data.frame(replicate(5,rnorm(N)))
  out$y = with(out$X,X1)
  return(out)
}
train = makeXY();
test  = makeXY();

Time = system.time({
  rf = randomForest(train$X[,1:3],train$y,ntree=1,sampsize=250,nodesize = 7,mtry =3,replace = FALSE)
})
plot(rf$predicted,train$X$X1)
#diagnosttics
rf
Time
preds = predict(rf,test$X)
plot(test$y,preds)
cor(preds,test$y)^2
rf

yt= train$y
sqrt(sum(
  (test$y-mean(train$y))^2 / (length(yt)-1)
))

#write datasets
# writeLines("#ifndef _FOO_H \n#define _FOO_H \nstatic float DData[] = {",con = "../data/train1.h")
# write.table(x=cbind(train$X,train$y),file="../data/train1.h",col.names=FALSE,row.names = FALSE,sep=",",append = TRUE)
# writeLines("\n} \n",con = file("../data/train1.h",open="at"))

file.remove("../shrubland_arduino/lib/testdata/src/testdata.h")
conn = file("../shrubland_arduino/lib/testdata/src/testdata.h",open="at")
writeLines("#ifndef _FOO_H \n#define _FOO_H \nstatic float DData[] = {",con = conn)
write.table(x=t(cbind(train$X,train$y)),file=conn,col.names=FALSE,row.names = FALSE,sep=",", eol = ",\n")
writeLines("}; \n\n static float DData_test[] = {",con = conn)
write.table(x=t(cbind( test$X, test$y)),file= conn,col.names=FALSE,row.names = FALSE,sep=",",eol = ",\n")
writeLines("}; \n#endif \n",con = conn)
close.connection(conn)

