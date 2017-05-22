Randomize
CreateObject("Wscript.Shell").Run "R CMD BATCH --vanilla --slave phylomad.Rscript" & " " & RND & " ", 0, False