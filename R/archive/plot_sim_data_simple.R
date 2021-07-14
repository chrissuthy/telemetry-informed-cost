plot(landscape[[sim]][[1]])

plot(scr_ss[[sim]], add=T)
points(traps)

lines(spatdata[inds][[1]])
lines(teldata[[sim]][inds][[1]], col = 4)

# started at 12:50