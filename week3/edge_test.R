library(edge)
data(kidney)
age <- kidney$age
sex <- kidney$sex
kidexpr <- kidney$kidexpr
de_obj <- build_study(data = kidexpr, adj.var = sex,
                      tme = age, sampling = "timecourse")
full_model <- fullModel(de_obj)
null_model <- nullModel(de_obj)