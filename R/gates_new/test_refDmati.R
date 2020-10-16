test_sbar_cell <- extract(local_ss_r, sbar, cellnumbers=T)[1]

Dac <- costDistance(
  tr1Corr,
  c(sbar_x, sbar_y),
  local_ss %>% select(x,y) %>% as.matrix())

table(Dmat_i[test_sbar_cell,] - Dac)




ref.Dmat_i <- function(from, to, local_ss_r, Dmat_i){
  from_cell <- as.numeric(raster::extract(x = local_ss_r, y = from, cellnumbers=T)[,1])
  to_cell <- as.numeric(raster::extract(x = local_ss_r, y = to, cellnumbers=T)[,1])
  result <- matrix(Dmat_i[from_cell, to_cell], nrow = nrow(from))
  return(result)
}

ref.Dmat_i(
  from = rbind(sbar, 
  to = local_ss %>% select(x,y) %>% as.matrix(), 
  local_ss_r, Dmat_i)

  
ref.Dmat_i(
    from = local_ss %>% filter(cell == s.grid[i-1]) %>% select(x,y) %>% as.matrix(),
    to = local_ss %>% select(x,y) %>% as.matrix(),
    local_ss_r,
    Dmat_i
  )
  