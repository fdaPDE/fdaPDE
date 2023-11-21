msg_file <- file(test_path("../data/msg.txt"), open = "wt")
sink(msg_file, type = "message")

# Building all reference solutions 
cat("-------------- \t \t \t Building Data \t \t \t --------------\n")
pb = txtProgressBar(min = 0, max = length(list.files(test_path("../data/setup/"))), 
                    initial = 0, style = 3)
i = 0
for( file in list.files(test_path("../data/setup/"))){
  path <- test_path("../data/setup", file)
  source(path)
  i <- i + 1
  setTxtProgressBar(pb, i)
}
close(pb)
cat("-------------- \t \t \t ------------- \t \t \t --------------\n")

sink(type="message")
close(msg_file)
unlink(test_path("../data/msg.txt"))