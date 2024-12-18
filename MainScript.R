##########################################################
##########################################################
##
## Marcus Suassuna Santos, November 2024
##  
## Perform analysis of hidden states for the Minimum flows time series
## in the Paraguay river basin
##
##########################################################
##########################################################

########################################################################
########################################################################
##
## 1: Prepare the environment
##
########################################################################
########################################################################

# Clean Global Environment
rm(list = ls(all = TRUE))

# Load the required libraries
library(tidyverse)
library(gridExtra)
library(lubridate)
library(ForestFit)
library(caret)
library(nnet)
library(randomForest)
library(kernlab)
library(MASS)
library(cowplot)
library(markovchain)
library(e1071)
library(rpart)
library(rpart.plot)
library(dismo)
library(mlbench)
library(gbm)
library(xgboost)
library(varSelRF)
library(gbm)
library(Gmisc)
library(pdp)
library(depmixS4)
library(plyr)
library(cowplot)
library(ggpubr)
library(diagram)
library(nnet) # For multinomial logistic regression
library(zoo) # Windowed standard deviation
library(ggdist) # For the rain plots
library(tidyquant)
cat('\014')

# Load used functions
source("scripts/select_data.R")
source("scripts/select_data_flow.R")
source("scripts/add_legend.R")

########################################################################
########################################################################
##
## 1: Read and correct water levels data
##
########################################################################
########################################################################

# Hidroweb data was downloaded in June, 14th, 2022
# The data was saved in the folder "data/hidroWeb"
sites <- substr(dir("data/HidroWeb"), 1, 8)

# Add the metadata
# Read the ANA inventory of the Brazilian stream gauges
streamGauges <- read_csv2("data/streamGauges.txt",
                          locale = readr::locale(encoding = "latin1"))

names(streamGauges) <- c("Site", "Name", "lat", "lon", "elevation", "DA")
streamGauges <- streamGauges %>%
  arrange(factor(Site, levels = sites))

# Read hydrological data and build daily water levels plots
# Annual miminum water levels is included
# In the stream gauge 66870000, stream flow data was used
# In all of the others, water level time series were used
df <- NULL

breaks_inicio <- rep(as.Date("1900-06-01"), length(sites))
for(i in 1:length(sites)){
  
  if(sites[i] == 66260002) next
  
  print(i)
  # Beginning of the dry-hydro-year is the month of March
  if(sites[i] %in% c(66270000, 66280000)){
    breaks <- seq(as.Date("1900-03-01"), length=135, by="year")
  } else {
    breaks <- seq(as.Date("1900-06-01"), length=135, by="year")
  }
  
  if(sites[i] == 66870000){
    
    # Read HidroWeb data
    Est <- read.csv2(paste0("data/HidroWeb/",sites[i],"_txt/",sites[i],"_Vazoes.txt"),
              skip = 15, sep = ";", row.names=NULL)
    
    # From wide to long, and other adjustments
    Est <- FlowsHidroWeb(Est)
    
    Est_Sup <- read.table(paste0("data/talemetry/",sites[i],".txt"),
                          header = TRUE) %>%
      mutate(DataHora = as.POSIXct(DataHora, tz = Sys.timezone())) %>%
      group_by(Data = as.Date(DataHora)) %>%
      dplyr::summarize(Vazao = mean(Vazao, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(NivelConsistencia = 1,
             Dia = day(Data),
             Mes = month(Data),
             Ano = year(Data),
             Dia_Ano = as.numeric((as.character(format(Data, "%j"))))) %>%
      dplyr::select("Data","NivelConsistencia","Dia","Mes","Ano","Vazao","Dia_Ano") %>%
      as.data.frame() %>%
      filter(Data > tail(Est$Data,1))
    
    Est <- rbind.data.frame(Est, Est_Sup) %>% arrange(Data)
    
    # Create plot
    png(paste0("plots/daily_recent/",sites[i],".png"),
        height = 500 * 5, width = 800 * 5, res = 75 * 10)
    par(mar = c(2.5,4.5,2.5,1))
    with(Est %>% filter(Data >= as.Date("2010-01-01")),
         plot(Data, Vazao,
              type = "l", col = "grey", bty = "n",
              main = "",
              xlab = "", ylab = "Discharge (m3/s)",
              xlim = c(as.Date("2010-01-01"),as.Date("2025-01-01"))))
    
    title(paste(streamGauges$Name[i]), adj = 0, line = 0)
    
    # Summarizes annual minima data
    Est_2 <- Est %>%
      filter(complete.cases(Vazao)) %>%
      mutate(hydroYear = cut(Data, breaks, labels=1900:2033)) %>%
      group_by(hydroYear) %>%
      dplyr::summarize(Minima = min(Vazao, na.rm = TRUE),
                       Dias = length(Vazao),
                       Data = Data[which.min(Vazao)]) %>%
      ungroup() %>%
      filter(Dias >= 20) %>%
      #filter(as.numeric(as.character(hydroYear)) < 2022) %>%
      mutate(Site = sites[i])
    
    # Include annual minima data in the daily plots
    with(Est_2, points(Data, Minima, pch = 20))
    
    dev.off()
    
    # Create plot
    png(paste0("plots/daily/",sites[i],".png"),
        height = 500 * 5, width = 800 * 5, res = 75 * 10)
    par(mar = c(2.5,4.5,2.5,1))
    with(Est,
         plot(Data, Vazao,
              type = "l", col = "grey", bty = "n",
              main = "",
              xlab = "", ylab = "Discharge (m3/s)"))
    
    title(paste(streamGauges$Name[i]), adj = 0, line = 0)
    
    # Summarizes annual minima data
    Est_2 <- Est %>%
      filter(complete.cases(Vazao)) %>%
      mutate(hydroYear = cut(Data, breaks, labels=1900:2033)) %>%
      group_by(hydroYear) %>%
      dplyr::summarize(Minima = min(Vazao, na.rm = TRUE),
                       Dias = length(Vazao),
                       Data = Data[which.min(Vazao)]) %>%
      ungroup() %>%
      filter(Dias >= 210) %>%
      filter(as.numeric(as.character(hydroYear)) < 2022) %>%
      mutate(Site = sites[i])
    
    # Include annual minima data in the daily plots
    with(Est_2, points(Data, Minima, pch = 20))
    abline(v = breaks,
           lty = 2, col = "lightgrey")
    
    dev.off()
    
    
  } else {
    
    # Read HidroWeb data
    Est <- read.csv2(paste0("data/HidroWeb/",sites[i],"_txt/",sites[i],"_Cotas.txt"),
                     skip = 15, sep = ";", row.names=NULL)
    
    # From wide to long, and other adjustments
    Est <- WaterLevelsHidroWeb(Est)
    
    telemetry_data <- substr(dir("data/talemetry"), 1, 8)
    telemetry_data <- telemetry_data[-which(telemetry_data == "66870000")]
    
    if(sites[i] %in% telemetry_data){
      
      Est_Sup <- read.table(paste0("data/talemetry/",sites[i],".txt"),
                 header = TRUE) %>%
        mutate(DataHora = as.POSIXct(DataHora, tz = Sys.timezone())) %>%
        group_by(Data = as.Date(DataHora)) %>%
        dplyr::summarize(Cota = mean(Nivel, na.rm = TRUE)) %>%
        ungroup() %>%
        mutate(NivelConsistencia = 1,
               Dia = day(Data),
               Mes = month(Data),
               Ano = year(Data),
               Dia_Ano = as.numeric((as.character(format(Data, "%j"))))) %>%
        dplyr::select("Data","NivelConsistencia","Dia","Mes","Ano","Cota","Dia_Ano") %>%
        as.data.frame() %>%
        filter(Data > tail(Est$Data,1))
      
      Est <- rbind.data.frame(Est, Est_Sup) %>% arrange(Data)
    }
    
    
    if(sites[i] %in% c(66825000, 67100000, 66970000)){
      
      Est_Sup <- read.csv2(paste0("data/talemetry/correction.csv")) %>%
        filter(Estacao == sites[i]) %>%
        mutate(Data = as.Date(Data, "%d/%m/%Y")) %>%
        group_by(Data) %>%
        dplyr::summarize(Cota = mean(Cota, na.rm = TRUE)) %>%
        ungroup() %>%
        mutate(NivelConsistencia = 1,
               Dia = day(Data),
               Mes = month(Data),
               Ano = year(Data),
               Dia_Ano = as.numeric((as.character(format(Data, "%j"))))) %>%
        dplyr::select("Data","NivelConsistencia","Dia","Mes","Ano","Cota","Dia_Ano") %>%
        as.data.frame() %>%
        filter(Data > tail(Est$Data,1))
      
      Est <- rbind.data.frame(Est, Est_Sup) %>% arrange(Data)
    }
    
    
    #with(Est_Sup[Est_Sup$Data > as.Date("1022-01-01"),],
    #     plot(Data, Cota, type = "l"))
    
    # Clear inconsistent data manually
    if(sites[i] == 66090000){
      Est[Est$Data > as.Date("2014-07-01") & Est$Data < as.Date("2014-09-01"),]$Cota <- NA
      Est[Est$Data > as.Date("2015-10-31") & Est$Data < as.Date("2015-11-15"),]$Cota <- NA
    } else if(sites[i] == 66125000) {
      Est[Est$Data > as.Date("1977-01-01") & Est$Data < as.Date("1977-02-01"),]$Cota <- NA
      Est[Est$Data == as.Date("2023-10-02"),]$Cota <- NA
      Est[Est$Data == as.Date("2023-10-24"),]$Cota <- NA
    } else if(sites[i] == 66750000) {
      Est[Est$Data > as.Date("1991-05-01") & Est$Data < as.Date("1991-07-01"),]$Cota <- NA
    } else if(sites[i] == 66910000) {
      Est[Est$Data > as.Date("1992-07-01") & Est$Data < as.Date("1992-09-01"),]$Cota <- NA
    }
    
    # Create plot
    png(paste0("plots/daily_recent/",sites[i],".png"),
        height = 500 * 5, width = 800 * 5, res = 75 * 10)
    par(mar = c(2.5,4.5,2,1))
    with(Est %>% filter(Data >= as.Date("2010-01-01")),
         plot(Data, Cota,
              type = "l", col = "grey", bty = "n",
              main = "",
              xlab = "", ylab = "Water level (cm)",
              xlim = c(as.Date("2010-01-01"),as.Date("2025-01-01"))))
    
    title(paste(streamGauges$Name[i]), adj = 0, line = 0)
    
    # Summarizes annual minima data
    Est_2 <- Est %>%
      filter(complete.cases(Cota)) %>%
      mutate(hydroYear = cut(Data, breaks, labels=1900:2033)) %>%
      group_by(hydroYear) %>%
      dplyr::summarize(Minima = min(Cota, na.rm = TRUE),
                       Dias = length(Cota),
                       Data = Data[which.min(Cota)]) %>%
      ungroup() %>%
      filter(complete.cases(hydroYear)) %>%
      filter(Dias >= 100) %>%
      #filter(as.numeric(as.character(hydroYear)) < 2024) %>%
      mutate(Site = sites[i],
             hydroYear = as.numeric(as.character(hydroYear)))
    
    # Include 2021 data manually in the sites were HidroWeb data is not updated
    # to 2021
    if(sites[i] == 67100000){
      
      Est_2 <- rbind.data.frame(Est_2,
                              data.frame("hydroYear" = c(2024),
                                         "Minima" = c(53),
                                         "Dias" = c(365),
                                         "Data" = c(as.Date("2024-10-24")),
                                         "Site" = c(67100000)))
      
    } else if(sites[i] == 66825000){
      
      Est_2 <- rbind.data.frame(Est_2,
                              data.frame("hydroYear" = c(2024),
                                         "Minima" = c(-69),
                                         "Dias" = c(365),
                                         "Data" = c(as.Date("2024-10-17")),
                                         "Site" = c(66825000)))
      
    } else if(sites[i] == 66070004){
      
      Est_2 <- rbind.data.frame(Est_2,
                              data.frame("hydroYear" = c(2024),
                                         "Minima" = c(33),
                                         "Dias" = c(365),
                                         "Data" = c(as.Date("2024-10-057")),
                                         "Site" = c(66070004)))
      
    } else if(sites[i] == 66260001){
      
      Est_2 <- rbind.data.frame(Est_2,
                              data.frame("hydroYear" = c(2024),
                                         "Minima" = c(91),
                                         "Dias" = c(365),
                                         "Data" = as.Date("2024-10-03"),
                                         "Site" = 66260001))
      
    } else if(sites[i] == 66470000){
      
      Est_2 <- rbind.data.frame(Est_2,
                              data.frame("hydroYear" = c(2021),
                                         "Minima" = c(36),
                                         "Dias" = c(365),
                                         "Data" = c(as.Date("2024-09-29")),
                                         "Site" = c(66470000)))
      
    } else if(sites[i] == 66910000){
      Est_2 <- Est_2
      Est_2 <- rbind.data.frame(Est_2,
                              data.frame("hydroYear" = c(2024),
                                         "Minima" = c(88),
                                         "Dias" = c(365),
                                         "Data" = c(as.Date("2024-09-14")),
                                         "Site" = c(66910000)))
      
    } else if(sites[i] == 66125000){
      Est_2 <- Est_2
      Est_2 <- rbind.data.frame(Est_2,
                                data.frame("hydroYear" = c(2024),
                                           "Minima" = c(226),
                                           "Dias" = c(365),
                                           "Data" = c(as.Date("2024-10-13")),
                                           "Site" = c(66125000)))
      
    } else if(sites[i] == 66945000){
      
      Est_2 <- rbind.data.frame(Est_2,
                              data.frame("hydroYear" = c(2024),
                                         "Minima" = c(150),
                                         "Dias" = c(365),
                                         "Data" = c(as.Date("2024-09-11")),
                                         "Site" = c(66945000)))
      
    } 
    
    # Include annual minima data in the daily plots
    with(Est_2, points(Data, Minima, pch = 20))
    #abline(v = breaks,
    #       lty = 2, col = "lightgrey")
    # with(Est, lines(Data, Minima, lty = 2))
    
    dev.off()
    
    
    # Create plot
    png(paste0("plots/daily/",sites[i],".png"),
        height = 500 * 5, width = 800 * 5, res = 75 * 10)
    par(mar = c(2.5,4.5,2,1))
    with(Est,
         plot(Data, Cota,
              type = "l", col = "grey", bty = "n",
              main = "",
              xlab = "", ylab = "Water level (cm)"))
    
    title(paste(streamGauges$Name[i]), adj = 0, line = 0)
    
    # Summarizes annual minima data
    Est_2 <- Est %>%
      filter(complete.cases(Cota)) %>%
      mutate(hydroYear = cut(Data, breaks, labels=1900:2033)) %>%
      group_by(hydroYear) %>%
      dplyr::summarize(Minima = min(Cota, na.rm = TRUE),
                       Dias = length(Cota),
                       Data = Data[which.min(Cota)]) %>%
      ungroup() %>%
      filter(complete.cases(hydroYear)) %>%
      filter(Dias >= 100) %>%
      filter(as.numeric(as.character(hydroYear)) < 2024) %>%
      mutate(Site = sites[i],
             hydroYear = as.numeric(as.character(hydroYear)))
    
    # Include 2021 data manually in the sites were HidroWeb data is not updated
    # to 2021
    if(sites[i] == 67100000){
      
      Est_2 <- rbind.data.frame(Est_2,
                                data.frame("hydroYear" = c(2024),
                                           "Minima" = c(53),
                                           "Dias" = c(365),
                                           "Data" = c(as.Date("2024-10-24")),
                                           "Site" = c(67100000)))
      
    } else if(sites[i] == 66825000){
      
      Est_2 <- rbind.data.frame(Est_2,
                                data.frame("hydroYear" = c(2024),
                                           "Minima" = c(-69),
                                           "Dias" = c(365),
                                           "Data" = c(as.Date("2024-10-17")),
                                           "Site" = c(66825000)))
      
    } else if(sites[i] == 66070004){
      
      Est_2 <- rbind.data.frame(Est_2,
                                data.frame("hydroYear" = c(2024),
                                           "Minima" = c(33),
                                           "Dias" = c(365),
                                           "Data" = c(as.Date("2024-10-07")),
                                           "Site" = c(66070004)))
      
    } else if(sites[i] == 66260001){
      
      Est_2 <- rbind.data.frame(Est_2,
                                data.frame("hydroYear" = c(2024),
                                           "Minima" = c(91),
                                           "Dias" = c(365),
                                           "Data" = as.Date("2024-10-03"),
                                           "Site" = 66260001))
      
    } else if(sites[i] == 66470000){
      
      Est_2 <- rbind.data.frame(Est_2,
                                data.frame("hydroYear" = c(2021),
                                           "Minima" = c(36),
                                           "Dias" = c(365),
                                           "Data" = c(as.Date("2024-09-29")),
                                           "Site" = c(66470000)))
      
    } else if(sites[i] == 66910000){
      Est_2 <- Est_2
      Est_2 <- rbind.data.frame(Est_2,
                                data.frame("hydroYear" = c(2024),
                                           "Minima" = c(88),
                                           "Dias" = c(365),
                                           "Data" = c(as.Date("2024-09-14")),
                                           "Site" = c(66910000)))
      
    } else if(sites[i] == 66945000){
      
      Est_2 <- rbind.data.frame(Est_2,
                                data.frame("hydroYear" = c(2024),
                                           "Minima" = c(150),
                                           "Dias" = c(365),
                                           "Data" = c(as.Date("2024-09-11")),
                                           "Site" = c(66945000)))
      
    } 
    
    # Include annual minima data in the daily plots
    with(Est_2, points(Data, Minima, pch = 20))
    # with(Est, lines(Data, Minima, lty = 2))
    
    abline(v = breaks,
           lty = 2, col = "lightgrey")
    
    dev.off()
    
    
    
  }
  
  df <- rbind(df, Est_2)

}

df <- df %>%
  mutate(hydroYear = as.numeric(as.character(hydroYear))) %>%
  arrange(hydroYear)

highlight <- data.frame(
  Site = "66825000",
  xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
)

# Useful plot for data visualization
p <- ggplot(df, aes(x = hydroYear, y = Minima)) +
  geom_line() +
  facet_wrap(~Site, scales = "free_y") +
  geom_rect(data = highlight, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = "red", size = 1.5, inherit.aes = FALSE) +
  labs(x = "Water year", y = "Annual Minima (cm)") +
  theme_bw()

png(paste0("plots/AnnualMinima.png"),
    height = 500 * 4, width = 800 * 4, res = 75 * 5)
print(p)
dev.off()

# Add the number length of the times series
siteN <- ddply(df,.(Site), summarise, N=length(Site))
streamGauges <- merge(streamGauges, siteN, by = "Site")

sum(streamGauges$N)
# 1373 minimum flows in the study area

df_unique <- df %>%
  distinct(hydroYear, Site, .keep_all = TRUE)

# Convert the data from long to wide
df_wide <- df_unique[,-c(3,4)] %>% 
  pivot_wider(names_from = Site, values_from = Minima, values_fill = list(Minima = NA)) %>%
  arrange(hydroYear)

# Scales wide data set
df_wide_scaled <- df_wide
df_wide_scaled[,-1] <- df_wide_scaled[,-1] %>% apply(2, scale)

########################################################################
########################################################################
##
## 2: Fit HMM - number of latent states
##
########################################################################
########################################################################

AIC_Obs <- array(NA, 5)
BIC_Obs <- array(NA, 5)
for(i in 1:5){
  # Removed inconsistent data 66710000 and 66260002
  set.seed(123)
  mod <- depmix(list(`66825000` ~ 1,`67100000` ~ 1,`66970000` ~ 1,
                     `66960008` ~ 1,`66070004` ~ 1,`66870000` ~ 1,
                     `66910000` ~ 1,`66941000` ~ 1,`67030000` ~ 1,
                     `66270000` ~ 1,`66600000` ~ 1,`66900000` ~ 1,
                     `66750000` ~ 1,`66810000` ~ 1,`66945000` ~ 1,
                     `66470000` ~ 1,`66650000` ~ 1,`66125000` ~ 1,
                     `66120000` ~ 1,`66090000` ~ 1,`66460000` ~ 1,
                     `66260001` ~ 1,`66010000` ~ 1,`66340000` ~ 1),
                data = df_wide_scaled,
                nstates = i,
                na.allow = TRUE,
                family = list(gaussian(),gaussian(),gaussian(),
                              gaussian(),gaussian(),gaussian(),
                              gaussian(),gaussian(),gaussian(),
                              gaussian(),gaussian(),gaussian(),
                              gaussian(),gaussian(),gaussian(),
                              gaussian(),gaussian(),gaussian(),
                              gaussian(),gaussian(),gaussian(),
                              gaussian(),gaussian(),gaussian()),
                instart = runif(i))
  
  fm_total <- fit(mod, verbose = FALSE, emc=em.control(rand=FALSE))
  AIC_Obs[i] <- AIC(fm_total)
  BIC_Obs[i] <- BIC(fm_total)
}

CompareFits <- data.frame("Classes" = c(1:5),
                          "BIC" = BIC_Obs)

CompareFits <- CompareFits %>%
  pivot_longer(!Classes,
               names_to = "InfoCriteria", values_to = "Value")

p <- ggplot(CompareFits, aes(x = Classes, y = Value, colour = InfoCriteria)) +
  geom_line(size = 3) +
  theme_bw() +
  ggtitle("Bayesian Information Criterion with different HMM fits",
          subtitle = "Information criteria evolution as the number of classes increases") +
  xlab("Number of low-water level regimes") +
  scale_color_manual(values=c("#023047")) +
  theme(legend.position="none")

p <- ggplot(CompareFits, aes(x = Classes, y = Value, colour = InfoCriteria)) +
  geom_line(size = 3) +
  theme_bw() +
  xlab("Number of low-water level regimes") +
  ylab("BIC estimate") +
  scale_color_manual(values=c("#023047")) +
  theme(legend.position="none")

png(paste0("plots/InfoCriteria.png"),
    height = 500 * 3, width = 800 * 3, res = 75 * 6)
print(p)
dev.off()

########################################################################
########################################################################
##
## 2: Fit HMM - actual fit
##
########################################################################
########################################################################

# Removed inconsistent data 66710000 and 66260002
set.seed(1)
mod <- depmix(list(`66825000` ~ 1,`67100000` ~ 1,`66970000` ~ 1,
                   `66960008` ~ 1,`66070004` ~ 1,`66870000` ~ 1,
                   `66910000` ~ 1,`66941000` ~ 1,`67030000` ~ 1,
                   `66270000` ~ 1,`66600000` ~ 1,`66900000` ~ 1,
                   `66750000` ~ 1,`66810000` ~ 1,`66945000` ~ 1,
                   `66470000` ~ 1,`66650000` ~ 1,`66125000` ~ 1,
                   `66120000` ~ 1,`66090000` ~ 1,`66460000` ~ 1,
                   `66260001` ~ 1,`66010000` ~ 1,`66340000` ~ 1),
              data = df_wide_scaled,
              nstates = 3,
              na.allow = TRUE,
              family = list(gaussian(),gaussian(),gaussian(),
                            gaussian(),gaussian(),gaussian(),
                            gaussian(),gaussian(),gaussian(),
                            gaussian(),gaussian(),gaussian(),
                            gaussian(),gaussian(),gaussian(),
                            gaussian(),gaussian(),gaussian(),
                            gaussian(),gaussian(),gaussian(),
                            gaussian(),gaussian(),gaussian()),
              instart = runif(3))

fm_total <- fit(mod, verbose = FALSE, emc=em.control(rand=FALSE))

# Visualize the fitted model
summary(fm_total)

# Removed inconsistent data 66710000 and 66260002
set.seed(1)
modLad <- depmix(list(`66825000` ~ 1),
              data = df_wide_scaled,
              nstates = 3,
              na.allow = TRUE,
              family = list(gaussian()),
              instart = runif(3))

fm_total_Lad <- fit(modLad, verbose = FALSE, emc=em.control(rand=FALSE))

# Check the number of stream gauges considered in each year to the multisite HMM
rowSums(!is.na(df_wide_scaled[,-1]))

# Check which class of drought was observed each year in the past
prstates_total <- apply(posterior(fm_total, type = "viterbi")[,c("S1", "S2","S3")], 1, which.max)
prstates_total_Lad <- apply(posterior(fm_total_Lad, type = "viterbi")[,c("S1", "S2","S3")], 1, which.max)

confusionMatrix(as.factor(prstates_total),
                as.factor(prstates_total_Lad))

# Check the runs of drought classes
Runs <- rle(prstates_total)
Runs$lengths[Runs$values==1]
Runs$lengths[Runs$values==2]
Runs$lengths[Runs$values==3]

# Estimate drought statistics at Ladário stram gauges
# The same can be done in any site, changing the code of the stream gauge
df_wide$class <- prstates_total
Medias <- df_wide %>%
  group_by(class) %>%
  dplyr::summarize(Medias = mean(`66825000`, na.rm = TRUE),
                   Sd = sd(`66825000`, na.rm = TRUE),
                   N = length(`66825000`)) %>%
  ungroup()

Medias$PClass <- Medias$N / sum(Medias$N)
Probability <- pnorm(-61, Medias$Medias[1], Medias$Sd[1]) * Medias$PClass[1]
ReturnPeriod <- 1/(pnorm(-61, Medias$Medias[1], Medias$Sd[1]) * Medias$PClass[1])

# Test if the residuals area Gaussian
TestNormality <- df_wide %>%
  group_by(class) %>%
  dplyr::summarize(Medias = mean(`66825000`, na.rm = TRUE),
                   Sd = sd(`66825000`, na.rm = TRUE),
                   Residual = `66825000` - mean(`66825000`, na.rm = TRUE)) %>%
  ungroup()

ggqqplot(TestNormality$Residual)
shapiro.test(TestNormality$Residual)
# Residuals look normally distributed

# Plot Ladário stream gauge
Paleta <- c(rgb(0.2,0.2,0.8,0.8),
            rgb(1,0.5,0.2,0.8),
            rgb(0.2,0.8,0.2,0.8),
            rgb(0.8,0.2,0.8,0.8))

Paleta <- Paleta[c(2,3,1,4)]

par(mfrow = c(1,1))

prstates_valores <- as.numeric(prstates_total)
prstates_valores[prstates_valores == 1] <- as.numeric(Medias[Medias$class == 1,2])
prstates_valores[prstates_valores == 2] <- as.numeric(Medias[Medias$class == 2,2])
prstates_valores[prstates_valores == 3] <- as.numeric(Medias[Medias$class == 3,2])
prstates_valores[prstates_valores == 4] <- as.numeric(Medias[Medias$class == 4,2])

png(paste0("plots/LadarioRegimes.png"),
    height = 500 * 5, width = 800 * 5, res = 75 * 5)
par(mar = c(4.5,4.5,1,7), oma = c(0,0,0,0), mfrow = c(1,1))
plot(df_wide$hydroYear, df_wide$`66825000`, type = "l", bty = "n", col = "grey",
     xlab = "Year", ylab = "Water level (cm)",
     main = "")

abline(h = seq(-100,400,100), col = "grey", lty = 2)
abline(v = seq(1900,2020,20), col = "grey", lty = 2)

points(df_wide$hydroYear, df_wide$`66825000`, pch = 20, cex = 2.5)
points(df_wide$hydroYear, df_wide$`66825000`, pch = 20, cex = 2,
       col = Paleta[prstates_total])

lines(df_wide$hydroYear, prstates_valores)
points(df_wide$hydroYear, prstates_valores,
       col = Paleta[prstates_total], pch = 20)

add_legend("right",
           legend = paste0(c("Dry Years: \n", "Normal Years: \n", "Wet Years: \n"),
                           "Mean = ", round(Medias$Medias,1), " cm \n",
                           "StdDev = ", round(Medias$Sd,1), " cm \n"),
           pch=20, 
           pt.cex = 2,
           text.font = 2,
           col=Paleta[1:3],
           horiz=FALSE, bty='n', cex=0.8)
dev.off()


png(paste0("plots/OnlyLadarioRegimes.png"),
    height = 500 * 5, width = 800 * 5, res = 75 * 5)
par(mar = c(4.5,4.5,1,7), oma = c(0,0,0,0), mfrow = c(1,1))
plot(df_wide$hydroYear, df_wide$`66825000`, type = "l", bty = "n", col = "grey",
     xlab = "Year", ylab = "Water level (cm)",
     main = "")

abline(h = seq(-100,400,100), col = "grey", lty = 2)
abline(v = seq(1900,2020,20), col = "grey", lty = 2)

points(df_wide$hydroYear, df_wide$`66825000`, pch = 20, cex = 2.5)
points(df_wide$hydroYear, df_wide$`66825000`, pch = 20, cex = 2,
       col = Paleta[prstates_total_Lad])

lines(df_wide$hydroYear, prstates_valores)
points(df_wide$hydroYear, prstates_valores,
       col = Paleta[prstates_total_Lad], pch = 20)

add_legend("right",
           legend = paste0(c("Dry Years: \n", "Normal Years: \n", "Wet Years: \n"),
                           "Mean = ", round(Medias$Medias,1), " cm \n",
                           "StdDev = ", round(Medias$Sd,1), " cm \n"),
           pch=20, 
           pt.cex = 2,
           text.font = 2,
           col=Paleta[1:3],
           horiz=FALSE, bty='n', cex=0.8)
dev.off()



# Rainplots

plotData <- data.frame("Year" = df_wide$hydroYear,
                       "Amf" = df_wide$`66825000`,
                       "Class" = prstates_total)

p <- plotData %>%
  ggplot(aes(x = factor(Class),
             y = Amf,
             fill = factor(Class))) +
  ggdist::stat_halfeye(
    adjust = 1,
    justification = -.1,
    .width = 0,
    point_colour = NA) +
  geom_boxplot(
    width = 0.12,
    alpha = 0.5) +
  ggdist::stat_dots(
    side = "left",
    justification = 1.1,
    binwidth = 5) +
  scale_fill_manual(values=c("#FF8033CC", "#33CC33CC", "#3333CCCC")) +
  theme_tq() +
  labs(
    title = "Annual minimum flows at the Ladário streamgauge",
    subtitle = "Classification on different drought classes was performed using a multisite HMM",
    x = "Drought Classes",
    y = "Annual mimimum flows (m³/s)",
    fill = "Drought classes"
  )
p

png(paste0("plots/RainPlots_Ladario.png"),
    height = 500 * 5, width = 800 * 5, res = 75 * 6)
print(p)
dev.off()


# Plot mixed distribution
R <- 5
png(paste("plots/mixedDistribution.png"), height = 500 * R, width = 800 * R, res = 75 * R)


dev.off()


# Plot mixed distribution
R <- 5
png(paste("plots/TimeSeries_mixedDistribution.png"), height = 600 * R, width = 800 * R, res = 75 * R)

# Dividir a área gráfica para incluir espaço para a legenda
layout(matrix(c(1, 2), ncol = 1), heights = c(3, 2)) # Primeiro gráfico maior, segundo menor
par(mar = c(4.5, 4.5, 1, 10)) # Aumenta margem direita para incluir a legenda

# Primeiro gráfico
plot(df_wide$hydroYear, df_wide$`66825000`, type = "l", bty = "n", col = "grey",
     xlab = "Year", ylab = "Water level (cm)", main = "")

abline(h = seq(-100, 400, 100), col = "grey", lty = 2)
abline(v = seq(1900, 2020, 20), col = "grey", lty = 2)

points(df_wide$hydroYear, df_wide$`66825000`, pch = 20, cex = 2.5)
points(df_wide$hydroYear, df_wide$`66825000`, pch = 20, cex = 2, col = Paleta[prstates_total])

lines(df_wide$hydroYear, prstates_valores)
points(df_wide$hydroYear, prstates_valores, col = Paleta[prstates_total], pch = 20)

mtext("(a)", side = 3, line = -0.5
      , adj = 0, font = 2, cex = 1.2)

# Adicionar a legenda na margem direita
par(xpd = NA) # Permitir desenhar fora da área do gráfico
legend(x = par("usr")[2] + 3, y = par("usr")[4], # Ajustar coordenadas manualmente
       legend = paste0(c("Dry Years: \n", "Normal Years: \n", "Wet Years: \n"),
                       "Mean = ", round(Medias$Medias,1), " cm \n",
                       "StdDev = ", round(Medias$Sd,1), " cm \n"),
       pt.cex = 2,
       col=Paleta[1:3],
       pch = 20, bty = "n", cex=0.8)

# Segundo gráfico
par(mar = c(4.5, 4.5, 1, 10)) # Mesmas margens para consistência
A <- c(-100:300)
plot(A, dnorm(A, Medias$Medias[1], Medias$Sd[1]), type = "l", bty = "n", ylim = c(0, 0.013),
     ylab = "Density", lwd = 3, xlab = "Water level (cm)", col = "#FF8033CC")

lines(A, dnorm(A, Medias$Medias[2], Medias$Sd[2]), col = "#33CC33CC", lwd = 3)
lines(A, dnorm(A, Medias$Medias[3], Medias$Sd[3]), col = "#3333CCCC", lwd = 3)

Dlad <- density(df_wide$`66825000`)
polygon(c(Dlad$x, rev(Dlad$x)), c(Dlad$y, rep(0, length(Dlad$x))), col = rgb(0.5, 0.5, 0.5, 0.3), border = "white")
lines(A, fixDistLad(A), col = 1, lwd = 3, lty = 2)

add_legend("right",
           legend = c("Dry years' PDF", "Normal years' PDF", "Wet years' PDF",
                      "Mixed distribution", "Empirical density"),
           pch=NA,
           col = c("#FF8033CC", "#33CC33CC", "#3333CCCC", rgb(0.5, 0.5, 0.5, 0.3), 1),
           lty = c(1,1,1,1,2), lwd = c(3, 3, 3, 8, 3),
           text.font = 2,
           horiz=FALSE, bty='n', cex=0.8)

mtext("(b)", side = 3, line = 0, adj = 0.1, font = 2, cex = 1.2)

dev.off()



# Violin plots of annual minima separated by drought class
df_prstates_total <- data.frame(hydroYear = df_wide$hydroYear,class = prstates_total)
df2 <- merge(df,df_prstates_total, by = "hydroYear", all.x = TRUE)
df2$class <- as.factor(df2$class)

for(i in 1:length(streamGauges$Site)){
  
  if(streamGauges$Site[i] == "66870000"){
    
    p <- ggplot(df2 %>% filter(Site == streamGauges$Site[i]),
                aes(x = as.factor(class), y = Minima, fill = as.factor(class))) +
      geom_violin() + theme_bw() +
      ylab("Discharge (m3/s)") + xlab("Drought classes") +
      ggtitle("Violin plots of annual minimum flows",
              subtitle = paste(streamGauges$Name[i], "stream gauge")) +
      scale_fill_discrete(name = "Drought class")
    
    png(paste0("plots/minimaViolin/",streamGauges$Site[i],".png"),
        height = 500 * 5, width = 800 * 5, res = 75 * 10)
    print(p)
    dev.off()
    
  } else {
    
    p <- ggplot(df2 %>% filter(Site == streamGauges$Site[i]),
                aes(x = as.factor(class), y = Minima, fill = as.factor(class))) +
      geom_violin() + theme_bw() +
      ylab("Water levels (cm)") + xlab("Drought classes") +
      ggtitle("Violin plots of annual minimum water levels",
              subtitle = paste(streamGauges$Name[i], "stream gauge")) +
      scale_fill_discrete(name = "Drought class")
    
    png(paste0("plots/minimaViolin/",streamGauges$Site[i],".png"),
        height = 500 * 5, width = 800 * 5, res = 75 * 10)
    print(p)
    dev.off()
  }
}


# Evaluate if differences between classes are meaningful
Estimates <- list()
p_values <- matrix(NA, nrow = length(streamGauges$Site), ncol = 3)
same_order <- array(NA, length(streamGauges$Site))

for(i in 1:length(streamGauges$Site)){
  
  if(streamGauges$Site[i] == 66260002) next
  
  fit_comparison <- lm(Minima ~ as.factor(class),
                       data = df2 %>% filter(Site == streamGauges$Site[i]))
  
  Estimates[[i]] <- fit_comparison %>% summary()
  p_values[i,] <- round(Estimates[[i]]$coefficients[,4], 4)
  
  values <- df2 %>%
    filter(Site == streamGauges$Site[i]) %>%
    group_by(class) %>%
    dplyr::summarize(Media = mean(Minima, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(class = as.numeric(as.character(class))) %>%
    arrange(class)
  
  values <- values$Media
  
  same_order[i] <- ifelse(values[1] < values[2] & values[3] > values[2], 1, 0)
  
}

names(Estimates) <- streamGauges$Site

streamGauges$Regimes <- rowSums(ifelse(p_values > 0.05,0,1))
names(streamGauges)[7] <- "N"

highlight <- data.frame(
  Site = unique(c(streamGauges$Site[same_order == 0],
                  streamGauges$Site[p_values[,2] > 0.05])),
  xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
)

# Useful plot for data visualization
ggplot(df, aes(x = hydroYear, y = Minima)) +
  geom_line() +
  facet_wrap(~Site, scales = "free_y") +
  geom_rect(data = highlight, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = "red", size = 1.5, inherit.aes = FALSE)

p <- ggplot(df2 %>% filter(Site != 66260002),
            aes(x = class, y = Minima, fill = class)) +
  geom_violin(colour = "white") + facet_wrap(~Site, scale = "free") +
  theme(panel.background = element_rect(fill = "white")) +
  scale_fill_manual(values = Paleta[1:3],
                    name = "Drought classes") +
  ylab("Annual minima (cm)") +
  xlab("Drought classes") +
  ggtitle("Violin plots of annual minimum water levels in all stream gauges for each drought class") +
  geom_rect(data = highlight, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = "red", size = 1.5, inherit.aes = FALSE)

png(paste0("plots/Regimes.png"),
    height = 500 * 5, width = 800 * 5, res = 75 * 5)
print(p)
dev.off()

# Classify drought regimes according to regional patterns
streamGauges$Regimes2 <- ifelse(streamGauges$Site %in% highlight$Site,
       "Different", "Regional")

write.csv(streamGauges, "outputTables/streamGauges.txt", row.names=FALSE)
table(streamGauges$Regimes2)

########################################################################
########################################################################
##
## 3: Evaluate the transition matrices of HMM
##
########################################################################
########################################################################

df_prstates_total$hydroYear <- as.numeric(as.character(df_prstates_total$hydroYear))

set.seed(123)
Drought_Regional_Total <- markovchainFit(data = df_prstates_total$class,
                                         method = "bootstrap", nboot = 50)

set.seed(123)
Drought_Regional_1900_1960 <- markovchainFit(data = df_prstates_total$class[df_prstates_total$hydroYear <= 1960],
                                             method = "bootstrap", nboot = 50)

set.seed(123)
Drought_Regional_1960_2020 <- markovchainFit(data = df_prstates_total$class[df_prstates_total$hydroYear > 1960],
                                             method = "bootstrap", nboot = 50)

Total <- rbind(Drought_Regional_Total$estimate[,],
               Drought_Regional_Total$confidenceInterval[2]$lowerEndpointMatrix,
               Drought_Regional_Total$confidenceInterval[3]$upperEndpointMatrix)

Regional_1900_1960 <- rbind(Drought_Regional_1900_1960$estimate[,],
                            Drought_Regional_1900_1960$confidenceInterval[2]$lowerEndpointMatrix,
                            Drought_Regional_1900_1960$confidenceInterval[3]$upperEndpointMatrix)

Regional_1960_2020 <- rbind(Drought_Regional_1960_2020$estimate[,],
                            Drought_Regional_1960_2020$confidenceInterval[2]$lowerEndpointMatrix,
                            Drought_Regional_1960_2020$confidenceInterval[3]$upperEndpointMatrix)

write.csv2(cbind(Total, Regional_1900_1960, Regional_1960_2020),
           "outputTables/transition.csv", row.names = FALSE)

Drought_Regional_1900_1960$estimate
Drought_Regional_1900_1960$confidenceInterval

Drought_Regional_1960_2020$estimate
Drought_Regional_1960_2020$confidenceInterval

# Explore transition probabilities
transitionProbability(Drought_Regional_Total$estimate, "3", "2")
transitionProbability(Drought_Regional_1900_1960$estimate, "1", "2")
transitionProbability(Drought_Regional_1960_2020$estimate, "2", "1")

# Build the transition diagrams
M1 <- round(Drought_Regional_Total$estimate[], 2)
M2 <- round(Drought_Regional_1900_1960$estimate[], 2)
M3 <- round(Drought_Regional_1960_2020$estimate[], 2)

colSums(t(M1))
colSums(t(M2))
colSums(t(M3))

Col <- M1
Col[] <- "darkred"
diag(Col) <- "darkred"

png(paste0("plots/Transition.png"),
    height = 500 * 5, width = 800 * 10, res = 75 * 10)
par(mar = c(3, 1, 1, 1), mfrow = c(1, 3), oma = c(0,0,1,0))
names <- c("Dry","Normal","Wet")
plotmat(t(M1), pos = c(2, 1), name = names, lwd = 1,
        box.lwd = 2, cex.txt = 0.8, box.size = 0.1,
        box.type = "circle", box.prop = 0.5,
        box.col = Paleta[c(1:3)],
        self.lwd = 2,
        self.shiftx = c(-0.1, 0.1, -0.1), self.cex = 0.5,
        arr.col = Col)

mtext("Transition matrix 1900 - 2021")
title(main = list("(a)", cex = 2, font = 2), adj = 0, line = 0, outer = TRUE)

plotmat(t(M2), pos = c(2, 1), name = names, lwd = 1,
        box.lwd = 2, cex.txt = 0.8, box.size = 0.1,
        box.type = "circle", box.prop = 0.5,
        box.col = Paleta[c(1:3)],
        self.lwd = 2,
        self.shiftx = c(-0.1, 0.1, -0.1), self.cex = 0.5,
        arr.col = Col)
mtext("Transition matrix 1900 - 1960")
title(main = list("(b)", cex = 2, font = 2), adj = 0.33, line = 0, outer = TRUE)

plotmat(t(M3), pos = c(2, 1), name = names, lwd = 1,
        box.lwd = 2, cex.txt = 0.8, box.size = 0.1,
        box.type = "circle", box.prop = 0.5,
        box.col = Paleta[c(1:3)],
        self.lwd = 2,
        self.shiftx = c(-0.1, 0.1, -0.1), self.cex = 0.5,
        arr.col = Col)
mtext("Transition matrix 1961 - 2021")
title(main = list("(c)", cex = 2, font = 2), adj = 0.66, line = 0, outer = TRUE)

par(mfrow = c(1,1))
dev.off()

# Explore other aspects of hidden states
initialState <- c(1,0,0)
after7Years <- initialState * (Drought_Regional_Total$estimate ^ 7)
after7Years
initialState * (Drought_Regional_Total$estimate ^ c(1:7))

Drought_Regional_1900_1960$estimate[1,1]^12*100
Drought_Regional_1960_2020$estimate[1,1]^12*100

fvals <- function(mchain,initialstate,n) {
  out<-data.frame()
  names(initialstate)<-names(mchain)
  for (i in 0:n){
    iteration<-initialstate*mchain^(i)
    out<-rbind(out,iteration)
  }
  out<-cbind(out, i=seq(0,n))
  out<-out[,c(4,1:3)]
  return(out)
}

CompareTransition <- data.frame("Years" = c(1:21),
                                "Period_Total" = fvals(mchain = Drought_Regional_Total$estimate, initialstate = c(100,0,0), n=20)[,2],
                                "Period_1900_1960" = fvals(mchain = Drought_Regional_1900_1960$estimate, initialstate = c(100,0,0), n=20)[,2],
                                "Period_1961_2021" = fvals(mchain = Drought_Regional_1960_2020$estimate, initialstate = c(100,0,0), n=20)[,2])

fvals(mchain = Drought_Regional_1960_2020$estimate, initialstate = c(100,0,0), n=20)
plot(fvals(mchain = Drought_Regional_1960_2020$estimate, initialstate = c(100,0,0), n=20)[,c(1,2)])

CompareTransition <- CompareTransition %>%
  pivot_longer(!Years,
               names_to = "Period", values_to = "Probability")

p1 <- ggplot(CompareTransition, aes(x = Years, y = Probability, colour = Period)) +
  geom_line(size = 2) +
  theme_bw() +
  ylim(0,100) +
  ggtitle("Probability that the region will leave extreme drought conditions",
          subtitle = "Different transition matrices were used for different calibration periods") +
  xlab("Years") +
  scale_color_manual(values=c("#2a9d8f", "#f4a261", "#e76f51")) +
  geom_hline(yintercept = 50, linetype="dashed")

png(paste0("plots/TransitionFromDroughts.png"),
    height = 500 * 3, width = 800 * 3, res = 75 * 5)
print(p1)
dev.off()

steadyStates(Drought_Regional_Total$estimate)
steadyStates(Drought_Regional_1900_1960$estimate)
steadyStates(Drought_Regional_1960_2020$estimate)

firstpassageKernel <- function(P, i, n){
  G <- P
  H <- P[i,]
  E <- 1 - diag(size(P)[2])
  for (m in 2:n) {
    G <- P %*% (G * E)
    H <- rbind(H, G[i,])
  }
  return(H)
}

firstPassage(object = Drought_Regional_Total$estimate, state = "1", n = 10)
# Interpretação: A probabilidade de o primeiro ano normal ser o segundo,
# dado que o ano atual é Seco, é de 24%

meanFirstPassageTime(Drought_Regional_1960_2020$estimate)
meanRecurrenceTime(Drought_Regional_Total$estimate)
meanRecurrenceTime(Drought_Regional_1960_2020$estimate)

firstPassagePdF.long <- firstPassage(object = Drought_Regional_Total$estimate, state = "2", n = 50)
P1 <- round(firstPassagePdF.long[47,], 5) * 100
firstPassagePdF.long <- firstPassage(object = Drought_Regional_1900_1960$estimate, state = "2", n = 50)
P2 <- round(firstPassagePdF.long[47,], 5) * 100
firstPassagePdF.long <- firstPassage(object = Drought_Regional_1960_2020$estimate, state = "2", n = 50)
P3 <- round(firstPassagePdF.long[47,], 5) * 100

c(P1[1], P2[1], P3[1])

# Probability with markovchain objects
conditionalDistribution(Drought_Regional_Total$estimate, "1")

########################################################################
########################################################################
##
## 4: Simulate with the Markov Chains
##
########################################################################
########################################################################

par(mar = c(4.5,4.5,2.5,1))
with(df_prstates_total,
     plot(hydroYear, class, bty = "n", pch = 20, cex = 2,
          xlim = c(1960, 2030)))

with(df_prstates_total,
     points(hydroYear, class, pch = 20, col = Paleta[class], cex = 1.5))

N <- 2000
Simulated_Total <- matrix(NA, ncol = 48, nrow = N)
Simulated_1900_1960 <- matrix(NA, ncol = 48, nrow = N)
Simulated_1961_2021 <- matrix(NA, ncol = 48, nrow = N)

for(i in 1:N){
  set.seed(i)
  weathersOfDays_1 <- rmarkovchain(n = 48,
                                   object = Drought_Regional_Total$estimate,
                                   t0 = "2")
  weathersOfDays_2 <- rmarkovchain(n = 48,
                                   object = Drought_Regional_1900_1960$estimate,
                                   t0 = "2")
  set.seed(i)
  weathersOfDays_3 <- rmarkovchain(n = 48,
                                   object = Drought_Regional_1960_2020$estimate,
                                   t0 = "2")
  
  Sim1 <- weathersOfDays_1[1:48] %>% as.character() %>% as.numeric()
  Sim2 <- weathersOfDays_2[1:48] %>% as.character() %>% as.numeric()
  Sim3 <- weathersOfDays_3[1:48] %>% as.character() %>% as.numeric()
  # lines(1974 + c(0:47), Sim, lwd = 4, col = rgb(0.5,0.5,0.5,0.05)
  
  Simulated_Total[i,] <- Sim1
  Simulated_1900_1960[i,] <- Sim2
  Simulated_1961_2021[i,] <- Sim3
  
}


Paleta <- c("#e76f51", "#e9c46a", "#264653", "#2a9d8f")

N_dry_3 <- apply(Simulated_1961_2021[,1:46], 1, function(x) sum(x == 1))
N_dry_3_2 <- ifelse(N_dry_3 == 0, -15, N_dry_3)
H <- hist(N_dry_3,
          nclass = 25)

H$counts[1]/sum(H$counts)


my_breaks <- H$breaks
my_colors <- rep(Paleta[3], length(my_breaks))
my_colors[1] <- rgb(0.8,0.2,0.2)

png(paste0("plots/DroughtYearsCount.png"),
    height = 500 * 3, width = 800 * 6, res = 75 * 6)
par(mfrow = c(1,3),
    mar = c(4.5,2.5,1.5,0),
    oma = c(0,2.5,0.5,1))
N_dry_1 <- apply(Simulated_Total[,1:46], 1, function(x) sum(x == 1))

hist(N_dry_1,
     border = "white",
     nclass = 25,
     freq = FALSE,
     main = "",
     xlab = "Number of simulated drought years",
     ylab = "",
     xlim = c(0,46),
     ylim = c(0,0.13),
     col = Paleta[1])
abline(v = mean(N_dry_1), col = 2, lwd = 3)
lines(c(1:46), dpois(c(1:46), mean(N_dry_1)), col = 2)
mtext("Transition matrix: 1900-2021", line = 0.5, cex = 0.8)

title(main = list("(a)", cex = 2, font = 2), adj = 0, line = 0)

N_dry_2 <- apply(Simulated_1900_1960[,1:46], 1, function(x) sum(x == 1))
hist(N_dry_2,
     border = "white", 
     nclass = 25,
     freq = FALSE,
     main = "",
     xlab = "Number of simulated drought years",
     ylab = "",
     xlim = c(0,46),
     ylim = c(0,0.13),
     col = Paleta[2])
abline(v = mean(N_dry_2), col = 2, lwd = 3)
lines(c(1:46), dpois(c(1:46), mean(N_dry_2)), col = 2)
mtext("Transition matrix: 1900-1960", line = 0.5, cex = 0.8)

title(main = list("(b)", cex = 2, font = 2), adj = 0, line = 0)

plot(H,
     border = "white", 
     freq = FALSE,
     main = "",
     xlab = "Number of simulated drought years",
     ylab = "",
     xlim = c(0,46),
     ylim = c(0,0.13),
     col = my_colors)

abline(v = mean(N_dry_3[N_dry_3!=0]), col = 2, lwd = 3)
#lines(c(1:46), dpois(c(1:46), mean(N_dry_3[N_dry_3!=0])), col = 2)
lines(c(1:46), dnorm(c(1:46),
                     mean(N_dry_3),
                     sd(N_dry_3)), col = 2)

mtext("Transition matrix: 1961-2021", line = 0.5, cex = 0.8)

title(main = list("(c)", cex = 2, font = 2), adj = 0, line = 0)

mtext("Density",
      outer = TRUE, line = 1, side = 2, cex = 0.8)
dev.off()

mean(N_dry_1)
mean(N_dry_2)

table(N_dry_1)[which.max(table(N_dry_1))]/N
table(N_dry_2)[which.max(table(N_dry_2))]/N
table(N_dry_3)[which.max(table(N_dry_3))]/N

ppois(0, mean(N_dry_1))
ppois(0, mean(N_dry_2)) * 100

# Teste
mean(N_dry_3[N_dry_3>0])
sd(N_dry_3[N_dry_3>0])

########################################################################
########################################################################
##
## 5: Precipitation
##
########################################################################
########################################################################

breaks <- seq(as.Date("1900-10-01"), length=130, by="year")

Catchment <- "67100000"
CHIRPSFolder <- dir("data/rainfall")

RainfallCHIRPS <- read.table(paste0("data/rainfall/",
                                    grep(Catchment, CHIRPSFolder, value = TRUE)),
                             header = FALSE,
                             col.names=c("Data","Rainfall","SD_pixels")) %>%
  mutate(Data = as.Date(Data),
         DataDia = as.Date(Data)) 

day(RainfallCHIRPS$Data) <- 1

RainfallCHIRPS <- RainfallCHIRPS %>%
  group_by(Data) %>%
  dplyr::summarize(Chuva = sum(Rainfall)) %>%
  ungroup() %>%
  mutate(Station = Catchment) %>%
  mutate(hydroYear = cut(Data, breaks, labels=1901:2029))

names(RainfallCHIRPS)[2] <- "Rain_CHIRPS"
RainfallCHIRPS$Source <- "CHIRPS"

RainfallGPCC <- read.csv2("data/gpcc/Chuvas_Paraguai.csv") %>%
  mutate(Data = as.Date(Data)) %>%
  filter(Station == Catchment)
names(RainfallGPCC)[2] <- "Rain_GPCC"
RainfallGPCC$Source <- "GPCC"

tail(RainfallCHIRPS)
head(RainfallGPCC)

Compara <- merge(RainfallCHIRPS, RainfallGPCC, by = "Data")
Compara$Season <- quarter(Compara$Data)
head(Compara)
Paleta <- c("#e76f51", "#e9c46a", "#264653", "#2a9d8f")

with(Compara,
     plot(Rain_GPCC, Rain_CHIRPS, pch = 20, bty = "n", col = Paleta[Season]))
fit <- lm(Rain_CHIRPS ~ Rain_GPCC - 1, data = Compara)
abline(fit, col = "red", lwd = 3)
summary(fit)
sqrt(mean(fit$residuals^2))

par(mfrow = c(2,2))
plot(fit)
par(mfrow = c(1,1))

head(RainfallCHIRPS)
head(RainfallGPCC)

names(RainfallCHIRPS)[2] <- "Rainfall"
names(RainfallGPCC)[2] <- "Rainfall"

Rainfall <- rbind.data.frame(RainfallGPCC %>% filter(Data < min(RainfallCHIRPS$Data)),
                             RainfallCHIRPS) %>%
  arrange(Data)

with(Rainfall, plot(Data, Rainfall, type = "l"))

STL_Rain <- stl(ts(Rainfall$Rainfall, frequency = 12, start = c(1891,1)),
                s.window = "periodic")

plot(STL_Rain,
     main = "Seasonal Decomposition of Paraguay river catchment precipitation by Loess")

plot(STL_Rain$time.series[,3])

WindowedSD <- rollapply(data = STL_Rain$time.series[,3], width=24, FUN=sd)


png(paste0("plots/RainfallComparison_1.png"),
    height = 800 * 6 , width = 800 * 6 * 2, res = 75 * 5 * 2)

par(mfrow = c(1,2),
    mar = c(4.5,4.5,2.5,1),
    oma = c(0,0,2,0))
with(Compara,
     plot(Rain_GPCC, Rain_CHIRPS, pch = 20, bty = "n", col = Paleta[Season],
          xlab = "Precipitation with GPCC (mm)",
          ylab = "Precipitation with CHIRPS (mm)",
          main = "CHIRPS against GPCC precipitation"))
fit <- lm(Rain_CHIRPS ~ Rain_GPCC - 1, data = Compara)
abline(0,1, col = 2, lwd = 3)
legend("bottomright", bty = "n",
       legend = c("Identity line", "1st quarter", "2nd quarter", "3rd quarter", "4th quarter"),
       pch = c(NA,20,20,20,20), lty = c(1,NA,NA,NA,NA),
       lwd = c(3,NA,NA,NA,NA),
       col = c(2,Paleta))
title(main = list("(a)", cex = 1.5, font = 2), adj = 0, line = 0, outer = TRUE)


plot(WindowedSD, bty = "n",
     ylab = "Standard Deviation (mm)",
     main = "Windowed SD of GPCC data")
fit <- loess(WindowedSD ~ seq_along(WindowedSD), span = 0.9)
LoessTS <- ts(fit$fitted, frequency = 12, start = c(1891,1))
lines(LoessTS, col = 2, lwd = 3)
legend("bottomright", bty = "n",
       legend = c("GPCC SD", "Loess fit"),
       lty = 1,
       lwd = c(1,3),
       col = c(1,2))
title(main = list("(b)", cex = 1.5, font = 2), adj = 0.5, line = 0, outer = TRUE)

dev.off()

breaks <- seq(as.Date("1899-10-01"), length=125, by="year")

YearRainfall <- Rainfall %>%
  group_by(hydroYear) %>%
  dplyr::summarize(Rainfall = sum(Rainfall),
                   Length = length(hydroYear)) %>%
  ungroup() %>%
  filter(Length == 12)
YearRainfall$hydroYear <- as.numeric(as.character(YearRainfall$hydroYear))

YearRainfall <- merge(YearRainfall, df_prstates_total) %>%
  filter(hydroYear >= 1940)
names(YearRainfall)[4] <- "Regime"
YearRainfall$Regime <- as.factor(YearRainfall$Regime)

png(paste0("plots/RainfallComparison_2.png"),
    height = 800 * 6 * 2, width = 800 * 6 * 2, res = 75 * 5 * 3)

par(mfrow = c(1,1),
    mar = c(4.5,4.5,2.5,1),
    oma = c(2,0,0,0))

layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))

with(YearRainfall,
     plot(hydroYear, Rainfall,
          type = "h", lwd = 6, bty = "n",
          col = "#264653", ylim = c(800,1500),
          ylab = "Precipitation (mm)",
          xlab = "",
          main = "Water year total precipitation"))

fit <- loess(Rainfall ~ hydroYear, data = YearRainfall, span = 0.6)
lines(YearRainfall$hydroYear, predict(fit, newdata = YearRainfall),
      lwd = 4, col = "#e76f51")
mtext("Water year",side=1,outer=F,cex=0.8, line = 2)
mtext("(a)", side=1, outer=F, line = 3.5, cex = 1.2)


with(Compara,
     plot(Rain_GPCC, Rain_CHIRPS, pch = 20, bty = "n", col = Paleta[Season],
          xlab = "Precipitation with GPCC (mm)",
          ylab = "Precipitation with CHIRPS (mm)",
          main = "CHIRPS and GPCC precipitation"))
fit <- lm(Rain_CHIRPS ~ Rain_GPCC - 1, data = Compara)
abline(0,1, col = 2, lwd = 3)
legend("bottomright", bty = "n",
       legend = c("Identity line", "1st quarter", "2nd quarter", "3rd quarter", "4th quarter"),
       pch = c(NA,20,20,20,20), lty = c(1,NA,NA,NA,NA),
       lwd = c(3,NA,NA,NA,NA),
       col = c(2,Paleta))


plot(WindowedSD, bty = "n",
     ylab = "Standard Deviation (mm)",
     main = "Windowed Standard Deviation of GPCC data")
fit <- loess(WindowedSD ~ seq_along(WindowedSD), span = 0.9)
LoessTS <- ts(fit$fitted, frequency = 12, start = c(1891,1))
lines(LoessTS, col = 2, lwd = 3)
legend("bottomright", bty = "n",
       legend = c("GPCC SD", "Loess fit"),
       lty = 1,
       lwd = c(1,3),
       col = c(1,2))

mtext('(b)',at=0.3, side=1,outer=T,cex=1.2)
mtext('(c)',at=0.8, side=1,outer=T,cex=1.2) 

dev.off()

png(paste0("plots/RainfallComparison_3.png"),
    height = 800 * 6, width = 800 * 6 * 2, res = 75 * 5 * 3)

par(mfrow = c(1,1),
    mar = c(4.5,4.5,2.5,1), oma = c(0,0,0,0))

with(YearRainfall,
     plot(hydroYear, Rainfall,
          type = "h", lwd = 6, bty = "n",
          col = "#264653", ylim = c(700,1500),
          ylab = "Precipitation (mm)",
          xlab = "",
          main = "Water year total precipitation"))
fit <- loess(Rainfall ~ hydroYear, data = YearRainfall, span = 0.6)
lines(YearRainfall$hydroYear, predict(fit, newdata = YearRainfall),
      lwd = 4, col = "#e76f51")
mtext("Water year",side=1,outer=F,cex=1, line = 2.5)

dev.off()

png(paste0("plots/RainfallComparison_4.png"),
    height = 800 * 12, width = 800 * 6 * 6, res = 75 * 5 * 6)

layout(matrix(c(1, 1, 2), nrow = 1, byrow = TRUE))
par(mar = c(4.5,4.5,1.5,1))
with(YearRainfall,
     plot(hydroYear, Rainfall,
          type = "h", lwd = 6, bty = "n",
          col = "#4178bc", ylim = c(700,1500),
          ylab = "Precipitation (mm)",
          xlab = "",
          main = ""))
fit <- loess(Rainfall ~ hydroYear, data = YearRainfall, span = 0.6)
lines(YearRainfall$hydroYear, predict(fit, newdata = YearRainfall),
      lwd = 4, col = "#eecc16")

fit <- lm(Rainfall ~ hydroYear, data = YearRainfall)
lines(YearRainfall$hydroYear, predict(fit, newdata = YearRainfall),
      lwd = 4, col = "#e8384f")
AFit <- anova(fit)

legend("topleft", bty = "n",
       legend = c("Observed rainfall", "Loess fit", "Linear trend"),
       lty = c(1,1,1),
       lwd = c(4,4,4),
       col = c("#4178bc","#eecc16", "#e8384f"))

title(main = list("(a)", cex = 1.5, font = 2), adj = 0, line = .5)

mtext("Water year",side=1,outer=F,cex=1, line = 2.5)

with(Compara,
     plot(Rain_GPCC, Rain_CHIRPS, pch = 20, bty = "n", col = Paleta[Season],
          xlab = "Precipitation with GPCC (mm)",
          ylab = "Precipitation with CHIRPS (mm)",
          main = ""))
fit <- lm(Rain_CHIRPS ~ Rain_GPCC - 1, data = Compara)
abline(0,1, col = 2, lwd = 3)
legend("bottomright", bty = "n",
       legend = c("Identity line", "1st quarter", "2nd quarter", "3rd quarter", "4th quarter"),
       pch = c(NA,20,20,20,20), lty = c(1,NA,NA,NA,NA),
       lwd = c(3,NA,NA,NA,NA),
       col = c(2,Paleta))

title(main = list("(b)", cex = 1.5, font = 2), adj = 0, line = .5)

dev.off()

"Water year total precipitation"
"CHIRPS and GPCC precipitation"

# Faz o gráfico
Paleta <- c(rgb(0.2,0.2,0.8,0.8),
            rgb(1,0.5,0.2,0.8),
            rgb(0.2,0.8,0.2,0.8),
            rgb(0.8,0.2,0.8,0.8))

Paleta <- Paleta[c(2,3,1,4)]

p1 <- ggplot(YearRainfall, aes(x = Regime, y = Rainfall, group = Regime)) +
  geom_violin(aes(fill = Regime)) +
  theme_bw() +
  ggtitle("Violin plots of annual cummulative precipitation in each low-water level regime") +
  xlab("Drought classes") +
  scale_fill_manual(values = Paleta[c(1:3)]) +
  theme(plot.title = element_text(size=8))

p1 <- ggplot(YearRainfall, aes(x = Regime, y = Rainfall, group = Regime)) +
  geom_violin(aes(fill = Regime)) +
  theme_bw() +
  xlab("Drought classes") +
  scale_fill_manual(values = Paleta[c(1:3)]) +
  theme(plot.title = element_text(size=8))

p1

p2 <- ggplot(YearRainfall, aes(x = Rainfall, group = Regime)) +
  geom_density(aes(fill = Regime), position = "fill") +
  ggtitle("Evolution of Density Kernels of different drought classes with annual precipitation") +
  ylab("Empirical probabilities") +
  xlab("Annual precipitation (mm)") +
  theme_classic() +
  scale_fill_manual(values = Paleta[c(1:3)], guide="none") +
  theme(plot.title = element_text(size=8))

p2 <- ggplot(YearRainfall, aes(x = Rainfall, group = Regime)) +
  geom_density(aes(fill = Regime), position = "fill") +
  ylab("Empirical probabilities") +
  xlab("Annual precipitation (mm)") +
  theme_classic() +
  scale_fill_manual(values = Paleta[c(1:3)], guide="none") +
  theme(plot.title = element_text(size=8))

p2

p_data <- layer_data(p2)

# These are the columns of interest    
head(p_data)
p_data$data[[1]]$density
p_data$data[[1]]$PANEL

Density <- data.frame("Grupo1" = p_data %>% filter(group == 1) %>% pull("density"),
                      "Grupo2" = p_data %>% filter(group == 2) %>% pull("density"),
                      "Grupo3" = p_data %>% filter(group == 3) %>% pull("density"),
                      "Prec" = p_data %>% filter(group == 1) %>% pull("x"),
                      "density" = p_data %>% filter(group == 1) %>% pull("ymin"),
                      "density" = p_data %>% filter(group == 2) %>% pull("ymin"))
DensityMax <- apply(Density %>% dplyr::select("Grupo1","Grupo2","Grupo3"), 1, which.max)
Density$Prec[max(which(DensityMax==1))]
Density$Prec[min(which(DensityMax==3))]

Density$density[max(which(DensityMax==1))]
Density$density[min(which(DensityMax==3))]

head(p_data[p_data$x >= 1108,])
P <- 1350

t.test(YearRainfall %>% filter(Regime == 1)  %>% pull("Rainfall"),
       YearRainfall %>% filter(Regime == 2)  %>% pull("Rainfall"))

t.test(YearRainfall %>% filter(Regime == 2)  %>% pull("Rainfall"),
       YearRainfall %>% filter(Regime == 3)  %>% pull("Rainfall"))


g1 <- ggplot(YearRainfall %>% filter(Regime == 1), aes(x = Rainfall)) +
  geom_density(aes(fill = Regime))
g1 <- ggplot_build(g1)[[1]][[1]]
g1 <- g1[which.min(abs(g1$x-P)),]$density

g2 <- ggplot(YearRainfall %>% filter(Regime == 2), aes(x = Rainfall)) +
  geom_density(aes(fill = Regime))
g2 <- ggplot_build(g2)[[1]][[1]]
g2 <- g2[which.min(abs(g2$x-P)),]$density

g3 <- ggplot(YearRainfall %>% filter(Regime == 3), aes(x = Rainfall)) +
  geom_density(aes(fill = Regime))
g3 <- ggplot_build(g3)[[1]][[1]]
g3 <- g3[which.min(abs(g3$x-P)),]$density

g1/(g1+g2+g3)
g2/(g1+g2+g3)
g3/(g1+g2+g3)

plot_row <- plot_grid(p2, p1,
                      labels = c("(a)", "(b)"),
                      label_size = 10,
                      ncol = 2, nrow = 1)


title <- ggdraw() + 
  draw_label(
    "Precipitation in different Low-water level regimes in the Paraguay river basin",
    fontface = 'bold',
    x = 0,
    hjust = 0
  )

subtitle <- ggdraw() + 
  draw_label(
    "Precipitation estimated in the entire catchment area with CHIRPS (Years >= 1981) and GPCC (Years < 1981)",
    fontface = 'italic',
    size = 10,
    x = 0,
    hjust = 0
  )


plot_grid(
  title, subtitle, plot_row,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 0.05, 1)
)


png(paste0("plots/rainfall_Classes.png"),
    height = 500 * 3, width = 800 * 6, res = 75 * 6)
plot_grid(
  #title, subtitle,
  plot_row,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 0.05, 1)
)
dev.off()

# Mean rainfall in each drought class
YearRainfall %>%
  group_by(Regime) %>%
  dplyr::summarize(Media = mean(Rainfall),
                   DP = sd(Rainfall)) %>%
  ungroup()


# Analisa chuva
YearRainfall$Scaled <- scale(YearRainfall$Rainfall)

N <- 7
AIC <- array(NA, length(N))
BIC <- array(NA, length(N))
for(i in 1:N){
  set.seed(i)
  mod <- depmix(Scaled ~ 1,
                data = YearRainfall,
                nstates = i,
                na.allow = TRUE,
                family = gaussian(),
                instart = runif(i))
  
  fm_total <- fit(mod, verbose = FALSE, emc=em.control(rand=FALSE))
  AIC[i] <- AIC(fm_total)
  BIC[i] <- BIC(fm_total)
}


CompareFits <- data.frame("Classes" = c(1:N),
                          "BIC" = BIC,
                          "AIC" = AIC)
CompareFits <- CompareFits %>%
  pivot_longer(!Classes,
               names_to = "InfoCriteria", values_to = "Value")

p <- ggplot(CompareFits, aes(x = Classes, y = Value, colour = InfoCriteria)) +
  geom_line(size = 3) +
  theme_bw() +
  ggtitle("Bayesian and Akaike information criteria with different HMM fits applied to rainfall data",
          subtitle = "Information criteria evolution as the number of classes increases") +
  xlab("Number of drought regimes") +
  scale_colour_discrete(name = "Information Criteria") +
  scale_color_manual(values=c("#e76f51", "#023047"))

png(paste0("plots/InfoCriteriaRainfall.png"),
    height = 500 * 3, width = 800 * 3, res = 75 * 5)
print(p)
dev.off()


mod <- depmix(Scaled ~ 1,
              data = YearRainfall,
              nstates = 3,
              na.allow = TRUE,
              family = gaussian(),
              instart = runif(3))

fm_total <- fit(mod, verbose = FALSE, emc=em.control(rand=FALSE))

prstates_total <- apply(posterior(fm_total, type = "viterbi")[,c("S1", "S2","S3")], 1, which.max)

Medias <- tapply(YearRainfall$Rainfall, prstates_total, mean)
SD <- tapply(YearRainfall$Rainfall, prstates_total, sd)

prstates_valores <- as.numeric(prstates_total)
prstates_valores[prstates_valores == 1] <- Medias[1]
prstates_valores[prstates_valores == 2] <- Medias[2]
prstates_valores[prstates_valores == 3] <- Medias[3]

# Faz o gráfico
Paleta <- c(rgb(0.2,0.2,0.8,0.8),
            rgb(1,0.5,0.2,0.8),
            rgb(0.2,0.8,0.2,0.8),
            rgb(0.8,0.2,0.8,0.8))

#Paleta <- Paleta[c(2,3,1,4)]

fit_comparison <- lm(Rainfall ~ as.factor(Regime)-1,
                     data = YearRainfall)
fit_comparison %>% summary()


YearRainfall$class <- prstates_total

png(paste0("plots/RainfallRegimes.png"),
    height = 500 * 5, width = 800 * 5, res = 75 * 5)
par(mar = c(4.5,4.5,0.5,7), mfrow = c(1,1))
plot(YearRainfall$hydroYear, YearRainfall$Rainfall, type = "l", bty = "n", col = "grey",
     xlab = "Year", ylab = "Precipitation (cm)",
     main = "")

abline(h = seq(0,2000,200), col = "grey", lty = 2)
abline(v = seq(1900,2020,20), col = "grey", lty = 2)

points(YearRainfall$hydroYear, YearRainfall$Rainfall, pch = 20, cex = 2.5)
points(YearRainfall$hydroYear, YearRainfall$Rainfall, pch = 20, cex = 2,
       col = Paleta[prstates_total])

lines(YearRainfall$hydroYear, prstates_valores)
points(YearRainfall$hydroYear, prstates_valores,
       col = Paleta[prstates_total], pch = 20)

add_legend("right",
           legend = paste0(c("Dry Years: \n", "Normal Years: \n", "Wet Years: \n"),
                           "Mean = ", round(Medias[c(2,3,1)]), " mm \n",
                           "StdDev = ", round(SD[c(2,3,1)]), " mm \n"),
           pch=20, 
           pt.cex = 2,
           col=Paleta[c(2,3,1)],
           horiz=FALSE, bty='n', cex=0.8)


dev.off()

png(paste0("plots/ACF_Rainfall.png"),
    height = 500 * 5, width = 800 * 5, res = 75 * 7)
par(mar = c(4.5,4.5,2.5,1))
acf(YearRainfall$Rainfall, main = "", bty = "n")
mtext("Autocorrelation function of annual precipitation in the PRB",
      side = 3, line = 1, font = 2, cex = 1.3, adj = 0)
mtext("Precipitation is estimated using CHIRPS and GPCC data",
      side = 3, line = 0, font = 3, cex = 0.8, adj = 0)
dev.off()

hist(YearRainfall$Rainfall, freq = FALSE, nclass = 10)
A <- c(800:1500)
lines(A, dnorm(A, mean(YearRainfall$Rainfall), sd(YearRainfall$Rainfall)),
      col = 2)
kurtosis(YearRainfall$Rainfall)
########################################################################
########################################################################
##
## 6: Prediction with Non-homogeneous Multinomial regression
##
########################################################################
########################################################################

YearRainfall$Regime_Prev <- lag(YearRainfall$Regime)
addmargins(table(YearRainfall$Regime_Prev))

# Calibrate with 7/9 of the data and validate with 2/9
set.seed(10)
TrainingLevels <- sample(18,14)
TestingLevels <- c(1:18)[!(c(1:18) %in% TrainingLevels)]

Set_1 <- YearRainfall %>% filter(Regime_Prev == 1)
Set_2 <- YearRainfall %>% filter(Regime_Prev == 2)
Set_3 <- YearRainfall %>% filter(Regime_Prev == 3)

png(paste0("plots/PrecipitationRegimes.png"),
    height = 500 * 5, width = 1000 * 5, res = 75 * 8)
par(mfrow = c(1,3), oma = c(2,2,1,0), mar = c(2.5,2.5,2.5,1))
with(Set_1, plot(Regime, Rainfall, ylim = c(900,1500), col = Paleta[c(2,3,1)],
                 ylab = "", xlab = "", names = c("Dry", "Normal", "Wet")))
mtext("Initial condition: DRY YEAR", line = 0.5, cex = 0.8)
title(main = list("(a)", cex = 1.5, font = 2), adj = 0, line = -0.5, outer = TRUE)

with(Set_2, plot(Regime, Rainfall, ylim = c(900,1500), col = Paleta[c(2,3,1)],
                 ylab = "", xlab = "", names = c("Dry", "Normal", "Wet")))
mtext("Initial condition: NORMAL YEAR", line = 0.5, cex = 0.8)
title(main = list("(b)", cex = 1.5, font = 2), adj = 0.33, line = -0.5, outer = TRUE)

with(Set_3, plot(Regime, Rainfall, ylim = c(900,1500), col = Paleta[c(2,3,1)],
                 ylab = "", xlab = "", names = c("Dry", "Normal", "Wet")))
mtext("Initial condition: WET YEAR", line = 0.5, cex = 0.8)
title(main = list("(c)", cex = 1.5, font = 2), adj = 0.66, line = -0.5, outer = TRUE)

mtext("Precipitation (mm)", side = 2, line = 0.5, cex = 0.8, outer = TRUE)
mtext("Next year condition", side = 1, line = 0.5, cex = 0.8, outer = TRUE)
#mtext("Precipitation in each low flow class and next year's class",
#      side = 3, line = 1.0, cex = 1.2, outer = TRUE, font = 2)
dev.off()

t.test(Set_1 %>% filter(Regime == 1)  %>% pull("Rainfall"),
       Set_1 %>% filter(Regime == 2)  %>% pull("Rainfall"))

t.test(Set_2 %>% filter(Regime == 1)  %>% pull("Rainfall"),
       Set_2 %>% filter(Regime == 2)  %>% pull("Rainfall"))

t.test(Set_2 %>% filter(Regime == 2)  %>% pull("Rainfall"),
       Set_2 %>% filter(Regime == 3)  %>% pull("Rainfall"))

t.test(Set_3 %>% filter(Regime == 2)  %>% pull("Rainfall"),
       Set_3 %>% filter(Regime == 3)  %>% pull("Rainfall"))


N <- 200
i <- 1
Acuracy <- matrix(NA, N, 3)
Sensitivity <- matrix(NA, N, 3)
Specificity <- matrix(NA, N, 3)

PrevRain <- seq(500,2500,10)
pp_class_1_1 <- matrix(NA, N, length(PrevRain))
pp_class_1_2 <- matrix(NA, N, length(PrevRain))
pp_class_1_3 <- matrix(0, N, length(PrevRain))
pp_class_2_1 <- matrix(NA, N, length(PrevRain))
pp_class_2_2 <- matrix(NA, N, length(PrevRain))
pp_class_2_3 <- matrix(NA, N, length(PrevRain))
pp_class_3_1 <- matrix(0, N, length(PrevRain))
pp_class_3_2 <- matrix(NA, N, length(PrevRain))
pp_class_3_3 <- matrix(NA, N, length(PrevRain))

Limite_1_1 <- array(NA, N)

Limite_2_1 <- array(NA, N)
Limite_2_3 <- array(NA, N)

Limite_3_2 <- array(NA, N)

for(i in 1:N){
  
  # Separating training and testing sets
  set.seed(i)
  index <- createDataPartition(Set_1$Regime, p = 7/9, list = FALSE)
  Training_Set_1 <- Set_1[index,] %>% arrange(hydroYear)
  Testing_Set_1 <- Set_1[-index,] %>% arrange(hydroYear)
  
  set.seed(i)
  index <- createDataPartition(Set_2$Regime, p = 7/9, list = FALSE)
  Training_Set_2 <- Set_2[index,] %>% arrange(hydroYear)
  Testing_Set_2 <- Set_2[-index,] %>% arrange(hydroYear)
  
  set.seed(i)
  index <- createDataPartition(Set_3$Regime, p = 7/9, list = FALSE)
  Training_Set_3 <- Set_3[index,] %>% arrange(hydroYear)
  Testing_Set_3 <- Set_3[-index,] %>% arrange(hydroYear)
  
  # Fit
  # https://www.datasciencecentral.com/alternatives-to-logistic-regression/
  multinom_model_1 <- multinom(Regime ~ Rainfall, data = Training_Set_1)
  multinom_model_2 <- multinom(Regime ~ Rainfall, data = Training_Set_2)
  multinom_model_3 <- multinom(Regime ~ Rainfall, data = Training_Set_3)
  
  # Summary
  #summary(multinom_model_1)
  #summary(multinom_model_2)
  #summary(multinom_model_3)
  
  #exp(coef(multinom_model_1))
  #exp(coef(multinom_model_2))
  #exp(coef(multinom_model_3))
  
  #round(fitted(multinom_model_1), 2)
  #round(fitted(multinom_model_2), 2)
  #round(fitted(multinom_model_3), 2)
  
  # Predictability
  Testing_Set_1$ClassPredicted <- predict(multinom_model_1, newdata = Testing_Set_1, "class")
  Testing_Set_2$ClassPredicted <- predict(multinom_model_2, newdata = Testing_Set_2, "class")
  Testing_Set_3$ClassPredicted <- predict(multinom_model_3, newdata = Testing_Set_3, "class")
  
  cM <- confusionMatrix(factor(Testing_Set_1$Regime, levels = c(1,2,3)),
                        factor(Testing_Set_1$ClassPredicted, levels = c(1,2,3)))
  
  Acuracy[i,] <- cM$overall[c(1,3,4)]
  Sensitivity[i,] <- cM$byClass[,1]
  Specificity[i,] <- cM$byClass[,2]
  
  # Predictability
  pp_class_1_1[i,] <- 1-predict(multinom_model_1,
                                newdata = data.frame("Rainfall" = PrevRain),
                                type = "probs", se = TRUE)
  pp_class_1_2[i,] <- predict(multinom_model_1,
                              newdata = data.frame("Rainfall" = PrevRain),
                              type = "probs", se = TRUE)
  
  Limite_1_1[i] <- PrevRain[which.min(abs(pp_class_1_2[i,]-0.5))]
  
  pp_class_2_1[i,] <- predict(multinom_model_2,
                              newdata = data.frame("Rainfall" = PrevRain),
                              type = "probs", se = TRUE)[,1]
  pp_class_2_2[i,] <- predict(multinom_model_2,
                              newdata = data.frame("Rainfall" = PrevRain),
                              type = "probs", se = TRUE)[,2]
  pp_class_2_3[i,] <- predict(multinom_model_2,
                              newdata = data.frame("Rainfall" = PrevRain),
                              type = "probs", se = TRUE)[,3]
  
  Limite_2_1[i] <- head(PrevRain[pp_class_2_1[i,] < pp_class_2_2[i,]],1)
  Limite_2_3[i] <- head(PrevRain[pp_class_2_2[i,] < pp_class_2_3[i,]],1)
  
  #par(mfrow = c(1,1))
  #plot(PrevRain,pp_class_2_1[i,], type = "l")
  #lines(PrevRain,pp_class_2_2[i,], col = 2)
  #lines(PrevRain,pp_class_2_3[i,], col = 3)
  
  pp_class_3_2[i,] <- 1-predict(multinom_model_3,
                              newdata = data.frame("Rainfall" = PrevRain),
                              type = "probs", se = TRUE)
  pp_class_3_3[i,] <- predict(multinom_model_3,
                                newdata = data.frame("Rainfall" = PrevRain),
                                type = "probs", se = TRUE)
  
  Limite_3_2[i] <- PrevRain[which.min(abs(pp_class_3_2[i,]-0.5))]
  
  #par(mfrow = c(1,1))
  #plot(PrevRain,pp_class_3_2[i,], type = "l", ylim = c(0,1))
  #lines(PrevRain,pp_class_3_3[i,], col = 2)
}

png(paste0("plots/ChangingProbabilities.png"),
    height = 500 * 6, width = 800 * 5, res = 75 * 5)
par(mfrow = c(3,1), oma = c(2.5,2.5,0,0), mar = c(2,2,1.5,1))
# From class 1
plot(PrevRain,pp_class_1_1[1,], type = "l", ylim = c(0,1),
     col = rgb(1,0.5,0,0.1), bty = "n",
     main = "Initial condition: DRY YEAR",
     ylab = "",
     xlab = "")
lines(PrevRain,pp_class_1_2[1,], col = rgb(0.1,0.8,0.1,0.1))
for(i in 2:N){
  lines(PrevRain,pp_class_1_1[i,], col = rgb(1,0.5,0,0.1))
  lines(PrevRain,pp_class_1_2[i,], col = rgb(0.1,0.8,0.1,0.1))
}

lines(PrevRain, colMeans(pp_class_1_1), col = rgb(1,0.5,0,1), lwd = 3)
lines(PrevRain, colMeans(pp_class_1_2), col = rgb(0.1,0.8,0.1,1), lwd = 3)

abline(v = mean(Limite_1_1), lty = 2)
polygon(c(qnorm(c(0.025, 0.975), mean(Limite_1_1), sd(Limite_1_1)),
          rev(qnorm(c(0.025, 0.975), mean(Limite_1_1), sd(Limite_1_1)))),
        c(0,0,1,1), col = rgb(0.5,0.5,0.5,0.3), border = NA)

legend("topright",
       lwd = 3, lty = 1,cex = 1.5,
       col = c(rgb(1,0.5,0,1), rgb(0.1,0.8,0.1,1), rgb(0.1,0.1,0.8,1)),
       legend = c("Dry", "Normal", "Wet"), box.col = "white")

title(main = list("(a)", cex = 1.5, font = 2), adj = 0, line = .5)


# From class 2
plot(PrevRain,pp_class_2_1[1,], type = "l", ylim = c(0,1),
     col = rgb(1,0.5,0,0.1), bty = "n",
     ylab = "",
     main = "Initial condition: NORMAL YEAR",
     xlab = "")
lines(PrevRain,pp_class_2_2[1,], col = rgb(0.1,0.8,0.1,0.1))
lines(PrevRain,pp_class_2_3[1,], col = rgb(0.1,0.1,0.8,0.1))
for(i in 2:N){
  lines(PrevRain,pp_class_2_1[i,], col = rgb(1,0.5,0,0.1))
  lines(PrevRain,pp_class_2_2[i,], col = rgb(0.1,0.8,0.1,0.1))
  lines(PrevRain,pp_class_2_3[i,], col = rgb(0.1,0.1,0.8,0.1))
}

lines(PrevRain, colMeans(pp_class_2_1), col = rgb(1,0.5,0,1), lwd = 3)
lines(PrevRain, colMeans(pp_class_2_2), col = rgb(0.1,0.8,0.1,1), lwd = 3)
lines(PrevRain, colMeans(pp_class_2_3), col = rgb(0.1,0.1,0.8,1), lwd = 3)

abline(v = mean(Limite_2_1), lty = 2)
abline(v = mean(Limite_2_3), lty = 2)

polygon(c(qnorm(c(0.025, 0.975), mean(Limite_2_1), sd(Limite_2_1)),
          rev(qnorm(c(0.025, 0.975), mean(Limite_2_1), sd(Limite_2_1)))),
        c(0,0,1,1), col = rgb(0.5,0.5,0.5,0.3), border = NA)

polygon(c(qnorm(c(0.025, 0.975), mean(Limite_2_3), sd(Limite_2_3)),
          rev(qnorm(c(0.025, 0.975), mean(Limite_2_3), sd(Limite_2_3)))),
        c(0,0,1,1), col = rgb(0.5,0.5,0.5,0.3), border = NA)

title(main = list("(b)", cex = 1.5, font = 2), adj = 0, line = .5)

# From class 3
plot(PrevRain,pp_class_3_2[1,], type = "l", ylim = c(0,1),
     col = rgb(0.1,0.8,0.1,0.1), bty = "n",
     main = "Initial condition: WET YEAR",
     ylab = "",
     xlab = "")
lines(PrevRain,pp_class_3_3[1,], col = rgb(0.1,0.1,0.8,0.1))
for(i in 2:N){
  lines(PrevRain,pp_class_3_2[i,], col = rgb(0.1,0.8,0.1,0.1))
  lines(PrevRain,pp_class_3_3[i,], col = rgb(0.1,0.1,0.8,0.1))
}

lines(PrevRain, colMeans(pp_class_3_2), col = rgb(0.1,0.8,0.1,1), lwd = 3)
lines(PrevRain, colMeans(pp_class_3_3), col = rgb(0.1,0.1,0.8,1), lwd = 3)

abline(v = mean(Limite_3_2), lty = 2)
polygon(c(qnorm(c(0.025, 0.975), mean(Limite_3_2), sd(Limite_3_2)),
          rev(qnorm(c(0.025, 0.975), mean(Limite_3_2), sd(Limite_3_2)))),
        c(0,0,1,1), col = rgb(0.5,0.5,0.5,0.3), border = NA)


mtext("Annual cummulative precipitation (mm)", side = 1, line = 0.5, outer = TRUE)
mtext("Transition probability", side = 2, line = 1, outer = TRUE)
#mtext("Probability of changing drought regime",
#      side = 3, line = 0.5, outer = TRUE, font = 2, cex = 1.2)

title(main = list("(c)", cex = 1.5, font = 2), adj = 0, line = .5)

dev.off()

########################################################################
########################################################################
##
## 5: Land Use
##
########################################################################
########################################################################

Files <- dir("tables")

Lu_Data <- vector(mode='list', length = 7)
for(i in 1:7){
  Lu_Data[[i]] <- read_csv2(paste0("tables/",Files[i]))
  Title <- strsplit(Files[i], "_Dias.csv")[[1]]
  
  png(paste0("plots/landuse/",Title,".png"), height = 500, width = 1000)
  
  compare <- merge(Lu_Data[[i]], plotData, by.x = "Time", by.y = "Year")
  
  par(mfrow = c(1,1))
  with(compare, plot(Class, LU, main = paste(Files[i])))
  
  par(mfrow = c(1,2), oma = c(0,0,2,0), mar = c(4.5,4.5,1,1), xpd=FALSE)
  with(compare,
       plot(Time, LU,
            main = "",
            type = "n",
            ylab = paste(Title, "coverage (%)"),
            bty = "n"))
  points(compare$Time, compare$LU, pch = 20, cex = 2,
         col = Paleta[compare$Class])
  
  legend("topleft",
         bty = "n",
         pch = 20, pt.cex = 3,
         col = Paleta[1:3],
         legend = c("Dry","Normal","Wet"))
  
  with(compare,
       plot(LU, Amf,
            col = Paleta[compare$Class],
            pch = 20, bty = "n", cex = 2,
            main = "",
            ylab = paste("Annual miminum flow (m³/s)"),
            xlab = paste(Title, "coverage (%)")))
  fit <- with(compare, lm(Amf ~ LU))
  summary(fit)
  abline(fit, col = 2, lwd = 3)
  
  mtext(paste("Relationship between",Title, ", Time and Annual miminum flow"),
        side = 3, outer = TRUE, cex = 1.2, font = 2)
  
  dev.off()
  
  
  p <- compare %>%
    ggplot(aes(x = factor(Class),
               y = LU,
               fill = factor(Class))) +
    geom_boxplot(
      width = 0.11,
      alpha = 0.5) +
    ggdist::stat_dots(
      side = "left",
      justification = 1.1) +
    scale_fill_manual(values=c("#FF8033CC", "#33CC33CC", "#3333CCCC")) +
    theme_tq() +
    labs(
      title = paste("Land use cover during different episodes of droughts - Product:", Title),
      subtitle = "Classification on different drought classes was performed using a multisite HMM, and land cover data was provided by Dias et al. (2016)",
      x = "Drought Classes",
      y = paste(Title, "coverage (%)"),
      fill = "Drought classes"
    ) +
    coord_flip()
  p
  
  png(paste0("plots/landuse/RainPlots_",Title,".png"),
      height = 500 * 5, width = 800 * 5, res = 75 * 6)
  print(p)
  dev.off()
  
}

names(Lu_Data) <- unlist(strsplit(Files, "_Dias.csv"))

for(i in 1:7){
  Lu_Data[[i]] <- read_csv2(paste0("tables/",Files[i]))
  Title <- strsplit(Files[i], "_Dias.csv")[[1]]
  
  compare <- merge(plotData, YearRainfall, by.x = "Year", by.y = "hydroYear")
  compare <- merge(Lu_Data[[i]], compare, by.x = "Time", by.y = "Year") %>%
    mutate(Scaled_LU = scale(LU))
  
  head(compare)
  fit1 <- loess(Amf ~ Scaled, data = compare)
  fit2 <- lm(Amf ~ Scaled_LU, data = compare)
  
  png(paste0("plots/landuse/Correlation_2_",Title,".png"), height = 500, width = 1000)
  
  par(mfrow = c(2,2), oma = c(0,0,2,0))
  
  plot(compare$Time, compare$Amf, pch = 20)
  lines(compare$Time, compare$Amf, col = "grey")
  
  plot(compare$Scaled, compare$Amf)
  lines(seq(-2,2,0.1), predict(fit1, newdata = data.frame("Scaled" = seq(-2,2,0.1))),
        lwd = 3, col = 2)
  
  plot(compare$Scaled_LU, compare$Amf)
  abline(fit2, col = 2, lwd = 3)
  
  plot(compare$Scaled_LU, fit1$residuals)
  fit3 <- lm(fit1$residuals ~ compare$Scaled_LU)
  abline(fit3, col = 2, lwd = 3)
  summary(fit3)
  
  mtext(paste(Title), outer = TRUE, font = 2, side = 3, line = 0.5)
  
  dev.off()
  
}

par(mfrow = c(1,1))
with(compare, plot(Class, LU, main = paste(Files[i])))
par(mfrow = c(1,2), oma = c(0,0,2,0), mar = c(4.5,4.5,1,1), xpd=FALSE)
with(compare,
     plot(Time, LU,
          main = "",
          type = "n",
          ylab = paste(Title, "coverage (%)"),
          bty = "n"))
points(compare$Time, compare$LU, pch = 20, cex = 2,
       col = Paleta[compare$Class])
legend("topleft",
         bty = "n",
         pch = 20, pt.cex = 3,
         col = Paleta[1:3],
         legend = c("Dry","Normal","Wet"))
  

with(compare,
       plot(LU, Amf,
            col = Paleta[compare$Class],
            pch = 20, bty = "n", cex = 2,
            main = "",
            ylab = paste("Annual miminum flow (m³/s)"),
            xlab = paste(Title, "coverage (%)")))
fit <- with(compare, lm(Amf ~ LU))
summary(fit)
abline(fit, col = 2, lwd = 3)
  
mtext(paste("Relationship between",Title, ", Time and Annual miminum flow"),
        side = 3, outer = TRUE, cex = 1.2, font = 2)
  
p <- compare %>%
  ggplot(aes(x = factor(Class),
             y = LU,
             fill = factor(Class))) +
  geom_boxplot(
    width = 0.11,
    alpha = 0.5) +
  ggdist::stat_dots(
    side = "left",
    justification = 1.1) +
  scale_fill_manual(values=c("#FF8033CC", "#33CC33CC", "#3333CCCC")) +
  theme_tq() +
  labs(
    title = paste("Land use cover during different episodes of droughts - Product:", Title),
      subtitle = "Classification on different drought classes was performed using a multisite HMM, and land cover data was provided by Dias et al. (2016)",
      x = "Drought Classes",
      y = paste(Title, "coverage (%)"),
      fill = "Drought classes"
    ) +
  coord_flip()

png(paste0("plots/landuse/RainPlots_",Title,".png"),
    height = 500 * 5, width = 800 * 5, res = 75 * 6)
print(p)
dev.off()
  
########################################################################
########################################################################
##
## 6: Teleconnections
##
########################################################################
########################################################################

source("scripts/climateIndices.R")

climaClasse <- merge(plotData, nino, by.x = "Year", by.y = "HydroYear")
climaClasse$Year.y <- NULL
names(climaClasse)[4] <- "nino"

climaClasse <- merge(climaClasse, pdo, by.x = "Year", by.y = "HydroYear")
climaClasse$Year.y <- NULL
names(climaClasse)[5] <- "pdo"

climaClasse <- merge(climaClasse, amo, by.x = "Year", by.y = "HydroYear")
climaClasse$Year.y <- NULL
names(climaClasse)[6] <- "amo"

climaClasse <- merge(climaClasse, YearRainfall %>% dplyr::select("hydroYear", "Rainfall"),
                     by.x = "Year", by.y = "hydroYear") %>%
  mutate(Rain_TH_1 = ifelse(Rainfall > 1221, 1, 0),
         Rain_TH_2 = ifelse(Rainfall > 1037, 1, 0))

climaClasse$Rain_TH_1 <- as.factor(climaClasse$Rain_TH_1)
climaClasse$Rain_TH_2 <- as.factor(climaClasse$Rain_TH_2)

climaClasse$Class <- as.factor(climaClasse$Class)
climaClasse$Class_2 <- lag(climaClasse$Class)
climaClasse$Class_2[1] <- 1

model_1 <- glm(Rain_TH_1 ~ nino + pdo + amo, family=binomial(link='logit'),data=climaClasse)
summary(model_1)
model_2 <- glm(Rain_TH_2 ~ nino + pdo + amo, family=binomial(link='logit'),data=climaClasse)
summary(model_2)

p_nino <- climaClasse %>%
  ggplot(aes(x = factor(Rain_TH_1),
             y = nino,
             fill = factor(Rain_TH_1))) +
  geom_boxplot(
    width = 0.11,
    alpha = 0.5) +
  ggdist::stat_halfeye(
    adjust = 0.5,
    justification = -.1,
    .width = 0,
    point_colour = NA) +
  ggdist::stat_dots(
    side = "left",
    justification = 1.1) +
  scale_fill_manual(values=c("#FF8033CC", "#33CC33CC", "#3333CCCC")) +
  theme_tq() +
  labs(
    title = paste("Land use cover during different episodes of droughts - Product:"),
    subtitle = "Classification on different drought classes was performed using a multisite HMM, and land cover data was provided by Dias et al. (2016)",
    x = "Drought Classes",
    y = paste("Nino"),
    fill = "Drought classes"
  ) +
  coord_flip()
p_nino


p_pdo <- climaClasse %>%
  ggplot(aes(x = factor(Rain_TH_1),
             y = pdo,
             fill = factor(Rain_TH_1))) +
  geom_boxplot(
    width = 0.11,
    alpha = 0.5) +
  ggdist::stat_halfeye(
    adjust = 0.5,
    justification = -.1,
    .width = 0,
    point_colour = NA) +
  ggdist::stat_dots(
    side = "left",
    justification = 1.1) +
  scale_fill_manual(values=c("#FF8033CC", "#33CC33CC", "#3333CCCC")) +
  theme_tq() +
  labs(
    title = paste("Land use cover during different episodes of droughts - Product:"),
    subtitle = "Classification on different drought classes was performed using a multisite HMM, and land cover data was provided by Dias et al. (2016)",
    x = "Drought Classes",
    y = paste("PDO"),
    fill = "Drought classes"
  ) +
  coord_flip()
p_pdo

p_amo <- climaClasse %>%
  ggplot(aes(x = factor(Rain_TH_1),
             y = amo,
             fill = factor(Rain_TH_1))) +
  geom_boxplot(
    width = 0.11,
    alpha = 0.5) +
  ggdist::stat_halfeye(
    adjust = 0.5,
    justification = -.1,
    .width = 0,
    point_colour = NA) +
  ggdist::stat_dots(
    side = "left",
    justification = 1.1) +
  scale_fill_manual(values=c("#FF8033CC", "#33CC33CC", "#3333CCCC")) +
  theme_tq() +
  labs(
    title = paste("Land use cover during different episodes of droughts - Product:"),
    subtitle = "Classification on different drought classes was performed using a multisite HMM, and land cover data was provided by Dias et al. (2016)",
    x = "Drought Classes",
    y = paste("AMO"),
    fill = "Drought classes"
  ) +
  coord_flip()
p_amo

str(climaClasse)
y1 <- climaClasse %>% filter(nino < -0.5 & pdo < 0) %>% pull(Year)
y2 <- climaClasse %>% filter(!(Year %in% c(y1))) %>% pull(Year)

c1<- climaClasse %>% filter(Year %in% y1) %>% pull(Rain_TH_1)
c2<- climaClasse %>% filter(Year %in% y2) %>% pull(Rain_TH_1)

sum(c1 == 1)/length(c1)
sum(c2 == 1)/length(c2)

y1 <- climaClasse %>% filter(nino < -0.5 & pdo < 0) %>% pull(Year)
y2 <- climaClasse %>% filter(!(Year %in% c(y1))) %>% pull(Year)

c1<- climaClasse %>% filter(Year %in% y1) %>% pull(Rain_TH_2)
c2<- climaClasse %>% filter(Year %in% y2) %>% pull(Rain_TH_2)

sum(c1 == 1)/length(c1)
sum(c2 == 1)/length(c2)

y1 <- climaClasse %>% filter(nino < -0.5 & pdo < -0.5) %>% pull(Year)
y2 <- climaClasse %>% filter(!(Year %in% c(y1))) %>% pull(Year)

c1<- climaClasse %>% filter(Year %in% y1) %>% pull(Rainfall)
c2<- climaClasse %>% filter(Year %in% y2) %>% pull(Rainfall)

t.test(climaClasse %>% filter(Year %in% y1)  %>% pull(Rainfall),
       climaClasse %>% filter(Year %in% y2)  %>% pull(Rainfall))

p_nino <- climaClasse %>%
  ggplot(aes(x = factor(Class),
             y = nino,
             fill = factor(Class))) +
  geom_boxplot(
    width = 0.11,
    alpha = 0.5) +
  ggdist::stat_halfeye(
    adjust = 0.5,
    justification = -.1,
    .width = 0,
    point_colour = NA) +
  ggdist::stat_dots(
    side = "left",
    justification = 1.1) +
  scale_fill_manual(values=c("#FF8033CC", "#33CC33CC", "#3333CCCC")) +
  theme_tq() +
  labs(
    title = paste("ENSO34 values during different drought episodes"),
    subtitle = "Classification on different drought classes was performed using a multisite HMM",
    x = "Drought Classes",
    y = paste("Nino"),
    fill = "Drought classes"
  ) +
  coord_flip()
p_nino

t.test(climaClasse %>% filter(Class == 1)  %>% pull("nino"),
       climaClasse %>% filter(Class == 2)  %>% pull("nino"))

t.test(climaClasse %>% filter(Class == 1)  %>% pull("nino"),
       climaClasse %>% filter(Class == 3)  %>% pull("nino"))

p_pdo <- climaClasse %>%
  ggplot(aes(x = factor(Class),
             y = pdo,
             fill = factor(Class))) +
  geom_boxplot(
    width = 0.11,
    alpha = 0.5) +
  ggdist::stat_halfeye(
    adjust = 0.5,
    justification = -.1,
    .width = 0,
    point_colour = NA) +
  ggdist::stat_dots(
    side = "left",
    justification = 1.1) +
  scale_fill_manual(values=c("#FF8033CC", "#33CC33CC", "#3333CCCC")) +
  theme_tq() +
  labs(
    title = paste("PDO values during different drought episodes"),
    subtitle = "Classification on different drought classes was performed using a multisite HMM",
    x = "Drought Classes",
    y = paste("PDO"),
    fill = "Drought classes"
  ) +
  coord_flip()
p_pdo

t.test(climaClasse %>% filter(Class == 1)  %>% pull("pdo"),
       climaClasse %>% filter(Class == 2)  %>% pull("pdo"))

t.test(climaClasse %>% filter(Class == 1)  %>% pull("pdo"),
       climaClasse %>% filter(Class == 3)  %>% pull("pdo"))


p_amo <- climaClasse %>%
  ggplot(aes(x = factor(Class),
             y = amo,
             fill = factor(Class))) +
  geom_boxplot(
    width = 0.11,
    alpha = 0.5) +
  ggdist::stat_halfeye(
    adjust = 0.5,
    justification = -.1,
    .width = 0,
    point_colour = NA) +
  ggdist::stat_dots(
    side = "left",
    justification = 1.1) +
  scale_fill_manual(values=c("#FF8033CC", "#33CC33CC", "#3333CCCC")) +
  theme_tq() +
  labs(
    title = paste("AMO values  during different drought episodes"),
    subtitle = "Classification on different drought classes was performed using a multisite HMM",
    x = "Drought Classes",
    y = paste("AMO"),
    fill = "Drought classes"
  ) +
  coord_flip()
p_amo


t.test(climaClasse %>% filter(Class == 1)  %>% pull("amo"),
       climaClasse %>% filter(Class != 1)  %>% pull("amo"))

t.test(climaClasse %>% filter(Class == 1)  %>% pull("amo"),
       climaClasse %>% filter(Class == 3)  %>% pull("amo"))


p <- plot_grid(p_nino, p_pdo, p_amo,
              labels = c("(a)", "(b)", "(c)"),
              label_size = 10,
              ncol = 1, nrow = 3)

png(paste0("plots/climateIndices.png"),
    height = 800 * 5, width = 500 * 7, res = 75 * 6)
print(p)
dev.off()

