
# Building all reference solutions 
for( file in list.files("tests/data/setup/")){
  path <- test_path("../data/setup", file)
  source(path)
}
