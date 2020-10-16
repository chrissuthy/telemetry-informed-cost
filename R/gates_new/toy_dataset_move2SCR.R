library(dplyr)
library(tidyr)
library(magrittr)

# Simulate some data (long-form)
input.data <- data.frame(
  id =    c(1,1,1, 2,2,2, 3,3,3),
  trap =  c(1,2,4, 3,3,3, 4,4,4),
  event = c(1,2,2, 2,2,2, 1,1,1)
)

capt.hist <- input.data %>%
  group_by(id, trap, event) %>%
  summarise(caps = n()) %>%
  distinct() # %>% mutate(caps = 1) If you want caps off

acast(capt.hist, id~trap~event, value.var = "caps", fill=0)




