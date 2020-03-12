Args<-commandArgs(T)
species=Args[1]
out_dir=Args[2]
datasetfile=Args[3]
setwd(out_dir)
dataset=read.table(file = datasetfile,header = F) 
#coverage=dataset[,4]
fraction=dataset[,6]
#fold_enrichment=dataset[,8]
ipdratio=dataset[,7]
#f_mean=mean(fraction)
dataset=subset(dataset,fraction>=0.7)
data1 = sample(x=dataset[,7],size=30)
t1=t.test(data1)$conf.int[1]
data2 = sample(x=dataset[,7],size=30)
t2=t.test(data1)$conf.int[1]
data3 = sample(x=dataset[,7],size=30)
t3=t.test(data1)$conf.int[1]
data4 = sample(x=dataset[,7],size=30)
t4=t.test(data1)$conf.int[1]
data5 = sample(x=dataset[,7],size=30)
t5=t.test(data1)$conf.int[1]
t=(t1+t2+t3+t4+t5)/5
print(t[1])



