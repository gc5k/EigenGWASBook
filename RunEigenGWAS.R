# make sure the gear.jar is in the path of your own data
RunEigenGWAS <- function(dataName, pc, breed){
  gear='java -jar gear.jar'
  EigenGWAS=paste(gear, "eigengwas", paste0("--", breed), "--bfile", dataName, "--ev", pc,  "--out", dataName)
  system(EigenGWAS) 
}