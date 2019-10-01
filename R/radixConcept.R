rm(list=ls())
b=10
ind=1:99
x1 = sample(ind)

i1 = x1%%b
ind1 = lapply(0:(b-1),function(j) ind[i1==j])
x2 = (x1-i1)/b


ind2 = lapply(ind1,function(ind)
  lapply(0:(b-1),function(j) {
    ind[i1==j]
  }
))
   
i2 = x2%%b
x3 = (x2-i2)/b




