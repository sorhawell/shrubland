rm(list=ls())
library(randomForest)
set.seed(1)
makeXY = function(N=1000) {
  out = list()
  out$X = data.frame(replicate(5,runif(N)))
  out$y = with(out$X,X1 + sin(X2) + 5 * X3*X4)
  return(out)
}
train = makeXY();
test  = makeXY();

Time = system.time({
  rf = randomForest(train$X,train$y,ntree=50,maxnodes = 500)
})

#diagnosttics
rf
Time
preds = predict(rf,test$X)
cor(preds,test$y)^2

#write datasets
write.table(x=cbind(train$X,train$y),file="../data/train1.csv",col.names=FALSE,row.names = FALSE,sep=",")
write.table(x=cbind( test$X, test$y),file= "../data/test1.csv",col.names=FALSE,row.names = FALSE,sep=",")
