# test bsub for running a job array
args <- commandArgs(TRUE)
task.id <- as.integer(args[1])
id <- paste("Juan",task.id,'.png',sep='')
print(id)
png(filename=id)
plot(1:100)
dev.off()
