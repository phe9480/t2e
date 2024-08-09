
#USE THE FOLLOWING FOR A NEW FOLDER

rm(list = ls())
remove.packages("t2e")

devtools::create("C:\\Users\\phe\\OneDrive - Daiichi Sankyo\\Training\\t2e")

devtools::document(pkg="C:\\Users\\phe\\OneDrive - Daiichi Sankyo\\Training\\t2e")
devtools::install(pkg="C:\\Users\\phe\\OneDrive - Daiichi Sankyo\\Training\\t2e") 

library(t2e)
help(package="t2e")


#######################################

#USE THE FOLLOWING FOR A NEW FOLDER

rm(list = ls())
remove.packages("corrTests")

devtools::create("C:/phe/corrTests")

devtools::document(pkg="C:/phe/corrTests")
devtools::install(pkg="C:/phe/corrTests") 

library(corrTests)
help(package="corrTests")

