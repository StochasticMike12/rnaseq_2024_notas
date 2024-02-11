## Pasos para crear este RStudio project con este archivo de notas.
usethis::create_project("~/rnaseq_2024_notas")
usethis::use_r("01-notas.R")


## Ejercicio 2
usethis::use_r("02-visualizar-mtcars.R")


## Ligar este proyecto a github
usethis::use_git()
usethis::use_github()
