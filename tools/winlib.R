# Build against mingw-w64 build of libnlopt2.4.2
if(!file.exists("../windows/nlopt/include/nlopt.h")){
  download.file("https://github.com/rwinlib/nlopt/archive/v2.4.2.zip", "lib.zip", quiet = TRUE)
  dir.create("../windows", showWarnings = FALSE)
  unzip("lib.zip", exdir = "../windows")
  unlink("lib.zip")
}

