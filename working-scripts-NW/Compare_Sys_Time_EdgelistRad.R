start <- Sys.time()
testRad_old <- ComputeEdgelist(sub,'x','y',rad.set = 25, k = NULL, nn.method = NULL)
end <- Sys.time()
old.time <- end - start

start <- Sys.time()
testRad_new <- CreateEdgelist(sub,radius = 25)
end <- Sys.time()
new.time <- end - start

old.time
new.time
