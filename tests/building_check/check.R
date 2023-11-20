install.packages('distro', repos='https://cran.stat.unipd.it/')

# get compiler -----------------------------------------------------------------
compiler <- "gcc"
if(!distro::distro()$id == "fedora"){
  system("dpkg --list | grep -i -c --color clang > compiler.txt")
  if(as.integer(read.table("compiler.txt", header = F)) > 0)
    compiler <- "clang"
}else{ # fedora
  system("yum list installed | grep -i -c --color clang > compiler.txt")
  if(as.integer(read.table("compiler.txt", header = F)) > 0)
    compiler <- "clang"
}
unlink("compiler.txt")
 
# get r status -----------------------------------------------------------------
r_status <- "release"
if(version$status == "Patched" || version$status == "RC"){
  r_status <- "patched"
}else if( version$status == "Under development (unstable)"){
  r_status <- "devel"
}

# run R CMD check --------------------------------------------------------------
rcmdcheck::rcmdcheck('fdaPDE/', 
                     check_dir= paste0('fdaPDE/tests/building_check/', 
                                       distro::distro()$id, '_', 
                                       compiler, '_',
                                       r_status))