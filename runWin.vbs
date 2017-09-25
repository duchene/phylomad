Randomize
CreateObject("Wscript.Shell").Run "R-Portable\App\R-Portable\bin\R.exe CMD BATCH --vanilla --slave phylomad.Rscript" & " " & RND & " ", 0, False