app:
	Rscript -e "shiny::runApp()"

render: 
	Rscript -e "rmarkdown::render('project2.Rmd')"