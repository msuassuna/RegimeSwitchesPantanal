library(tidyverse)
library(dplyr)

evtp <- read.table("data/evtp/1d_evapo_ssebop_67100000_PORTO_MURTINHO.txt",
           sep = "\t") %>%
  rename(DataHora = V1,
         EVTP = V2,
         DP_EVTP = V3) %>%
  dplyr::select("DataHora","EVTP", "DP_EVTP") %>%
  mutate(DataHora = as.Date(DataHora, "%Y-%m-%dT%H:%M:%S"))
  
# Criar o gráfico
p  <- ggplot(evtp, aes(x = as.Date(DataHora), y = EVTP)) +
  geom_line(color = "darkblue", linewidth = 1.5) +
  geom_smooth(method = "loess", color = "red", fill = "pink", alpha = 0.3) + # Linha de tendência LOESS
  labs(x = "Time",
       y = "Evapotranspiration (mm)") +
  theme_minimal()


png("plots/evtp.png", height = 500, width = 800)
print(p)
dev.off()
