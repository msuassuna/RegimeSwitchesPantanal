# Read HidroWeb data
Est <- read.csv2(unz(paste0("data/HidroWeb/",67100000,"/vazoes_T_",67100000,".zip"),
                     paste0("vazoes_T_",67100000,".txt")),
                 skip = 12, sep = ";", row.names=NULL)

# It requires adjustments on the header
RepNomes <- names(Est)[2:48]
Est <- Est[1:47]
names(Est) <- RepNomes

# From wide to long, and other adjustments
Est <- FlowsHidroWeb(Est)
plot(Est$Data, Est$Vazao, type = "l")

Est <- Est %>%
  filter(complete.cases(Vazao)) %>%
  mutate(hydroYear = cut(Data, breaks, labels=1900:2023)) %>%
  group_by(hydroYear) %>%
  dplyr::summarize(Minima = min(Vazao, na.rm = TRUE),
                   Media = mean(Vazao, na.rm = TRUE),
                   Dias = length(Vazao),
                   Data = Data[which.min(Vazao)]) %>%
  ungroup() %>%
  filter(Dias >= 210) %>%
  filter(as.numeric(as.character(hydroYear)) < 2022) %>%
  mutate(Site = 67100000)
abline(h = mean(Est$Media), col = 2)

Est <- merge(YearRainfall, Est %>% dplyr::select(hydroYear, Media), by = "hydroYear")
Est$Media_mm <- Est$Media * 86400 / (streamGauges$DA[streamGauges$Site == 67100000] * 10^3) * 365

plot(Est$Media_mm/Est$Rainfall, type = "l")
mean(Est$Media_mm/Est$Rainfall)

2350 * 86400 / (streamGauges$DA[streamGauges$Site == 67100000] * 10^3) * 365
