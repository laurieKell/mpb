t0=Sys.time()
od="C:/Users/lauri/AppData/Local/Temp/spictcheck"
dir.create(od, showWarnings=FALSE)
res=tryCatch(rmarkdown::render("spict.Rmd", output_dir=od,
                               output_format=rmarkdown::md_document(),
                               intermediates_dir=od, clean=FALSE),
             error=function(e) {cat("RENDER ERROR:\n"); print(e); NULL})
print(res)
print(Sys.time()-t0)
cat("DONE\n")
