aaa <- seq(1, 90*24)
out <- base::cut(x = aaa, breaks = seq(1, 90*24, by = 24))

as.numeric(out)
