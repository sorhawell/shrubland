library(xgboost)
library(data.table)

dtrain = fread("../data/train1.csv") 
dtrain[ , V6 := V6 + rnorm(1000)*.02]
mtrain = as.matrix(dtrain[,-6])
dtest  = fread("../data/test1.csv")
dtest[ , V6 := V6 + rnorm(1000)*.02]
mtest  = as.matrix(dtest[,-6])

log(10)
etas = exp(seq(log(.05),log(1),length.out = 25))
out = lapply(etas, function(x) {
  xg = xgboost::xgboost(mtrain,dtrain$V6,nrounds=50, params = list(
    eta=x,max_depth=2
  ))
  preds = predict(xg,mtest)
  out = cor(preds,dtest$V6)
  print(out)
  return(out)
})
8*50*4

md=4
ntree=50
xg = xgboost::xgboost(mtrain,dtrain$V6,nrounds=ntree, params = list(
  eta=.85,max_depth=md
))
preds = predict(xg,mtest)
plot(etas,unlist(out))


xl = xgboost::xgb.dump(xg)
object.size(xl)/1400

idx = cumsum(grepl("booster",xl))
td = lapply(1:max(idx),function(i) {
  xl[i==idx][-1]
})


msize = 2^(md+1)
leaf_start = 2^md

trees = lapply(1:ntree,function(i) {

  td1 = td[[i]]
  c1  =1
  c2 = regexpr(":",td1)-1
  td1 = td1[order(as.numeric(substr(td1,1,c2)))]
  
  
  
  
  
  splt= rep(NA,msize)
  vars = rep(NA,msize)
  i_leaf=1
  leafs = rep(NA,msize)
  lrvec = rep(NA,msize)
  dvec = rep(NA,msize)
  lrval = rep(NA,msize)
  
  parser = function(j,d) {
      
    
    x =td1[j]
    #print(x)
    v0 = regexpr("\\[f",x)[1]
    #print(j)
    dvec[j]<<-d
    if(v0>=1) {
      v1 = v0+2
      v2 = regexpr("<",x)-1
      vars[j] <<- as.numeric(substr(x,v1,v2))
      s1 = v2+2
      s2 = regexpr("\\] ",x)-1
      splt[j] <<- as.numeric(substr(x,s1,s2))
      i<<-i+1
      
      temp_val= 1
      if(d>0) {
        for(ii in 1:d) temp_val = temp_val + (1+lrvec[ii])*2^(d-ii)
      }
      # cat("\n",temp_val,d,"\n")
      # print(lrvec)
      # 
       lrval[j] <<- temp_val  
      
      
      ##go left
      l1 = regexpr("yes=",x)+4
      l2 = regexpr(",no=",x)-1
      newj=as.numeric(substr(x,l1,l2))+1
      lrvec[d+1] <<- 0
      newd =d+1
      parser(newj,newd)
      
      
      #go right
      r1 =l2+5
      r2 = regexpr(",missing=",x)-1
      newj=as.numeric(substr(x,r1,r2))+1
      lrvec[d+1] <<- 1
      parser(newj,newd)
      
    } else {
      lf1  = regexpr("leaf=",x)+5
      lf2  = nchar(x)
      leaf = as.numeric(substr(x,lf1,lf2))+1
      leafs[j]<<-leaf
      
      
      temp_val= leaf_start
      if(d>0) {
        for(ii in 1:d) temp_val = temp_val + (lrvec[ii])*2^(d-ii)
      }
      # cat("\n",temp_val,d,"\n")
      # print(lrvec)
      lrval[j] <<- temp_val  
      
    }
  }
  
  parser(j,0)
  
  # dres = data.frame(
  # splt,
  # vars,
  # leafs,
  # lrval,
  # dvec
  # )[order(lrval),]
  
  list(
  splitvals = splt[match(1:(2^md-1),dres$lrval)],
  splitvar  = as.integer(vars[match(1:(2^md-1),dres$lrval)]),
  leafs     = leafs[match(leaf_start:(msize-1),dres$lrval)]
  )
})


trees[[32]]
