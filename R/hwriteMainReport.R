hwriteMainReport <- function(dirname="result", title="RBDmap", links=list()) {
  dir.create(dirname, showWarnings=FALSE)
  if (!file.exists(file.path(dirname,"hwriter.css"))) {
    file.copy(system.file("images","hwriter.css",package="hwriter"),
              file.path(dirname,"hwriter.css"))
  }
  file.copy(system.file(file.path("MainReport","index.html"), package="RBDmap"),
            file.path(dirname,"index.html"), overwrite=TRUE)
  L = readLines(file.path(dirname,"index.html"))
  L = gsub("TITLE", title, L, fixed=TRUE)
  writeLines(L, file.path(dirname,"index.html"))
  file.copy(system.file(file.path("MainReport","frameLeft.html"), package="RBDmap"),
            file.path(dirname,"frameLeft.html"), overwrite=TRUE)
  Links = sprintf("<a href=%s target=main>%s</a><br>",links, names(links))
  Links[is.na(links)] = "<br>"
  Links = paste(Links, collapse="\n")
  L = readLines(file.path(dirname,"frameLeft.html"))
  L = gsub("TITLE", title, L, fixed=TRUE)
  L = gsub("LINKS", Links, L, fixed=TRUE)
  writeLines(L, file.path(dirname,"frameLeft.html"))
  invisible(NULL)
}
