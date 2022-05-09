# Some vignettes are too slow to run on CRAN, etc.
# These vignettes are given a .orig suffix
# Running this makefile will knit .orig to .Rmd
# Result is an .Rmd with output but no running code so it will compile instantly 
# You must include the generated figure files in the R package, git repo, etc.

.PHONY: all clean

all: colext.Rmd

colext.Rmd: colext.Rmd.orig
	Rscript -e "knitr::knit('colext.Rmd.orig', output='colext.Rmd')"
	# skip vignette name check
	echo '' >> colext.Rmd
	echo '```{r, echo=FALSE}' >> colext.Rmd
	echo 'options(rmarkdown.html_vignette.check_title = FALSE)' >> colext.Rmd
	echo '```' >> colext.Rmd

clean:
	rm -rf colext.Rmd

