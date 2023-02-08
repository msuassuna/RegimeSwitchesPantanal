##########################################################
##########################################################
##
## Marcus Suassuna Santos, June 2022
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

# Load packages
library(tidyverse)
library(lubridate)
library(ForestFit)
library(caret)
library(cowplot)
library(markovchain)
library(e1071)
library(Gmisc)
library(depmixS4)
library(plyr)
library(cowplot)
library(ggpubr)
library(diagram)
library(zoo) # Windowed standard deviation
library(ggplot2);cat('\014')

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
sites <- dir("data/HidroWeb")

# Add the metadata
# Read the ANA inventory of the Brazilian stream gauges
# Metadata is saved as a .mdb file and requires RODBC package to read
arquivo <- RODBC::odbcConnectAccess2007("data/hidroWebMetaData/HIDRO.mdb")
streamGauges <- RODBC::sqlFetch(arquivo,"Estacao",as.is=TRUE)
streamGauges <- streamGauges %>%
  dplyr::select("Codigo", "Nome","Latitude","Longitude","Altitude","AreaDrenagem", "RioCodigo", "EstadoCodigo") %>%
  filter(Codigo %in% sites) %>%
  arrange(Codigo) %>%
  dplyr::select(Codigo, Nome, Latitude, Longitude, Altitude, AreaDrenagem)
rm(arquivo)

names(streamGauges) <- c("Site", "Name", "lat", "lon", "elevation", "DA")
streamGauges <- streamGauges %>%
  arrange(factor(Site, levels = sites))

# In this loop the data is read. Also, daily plots are built with
# annual miminum water levels is included
# In the stream gauge 66870000, stream flow data was used
# In all of the others, water level time series were used
df <- NULL
for(i in 1:length(sites)){
  
  print(i)
  # Beginning of the dry-hydro-year is the month of March
  breaks <- seq(as.Date("1900-03-01"), length=125, by="year")
  
  if(sites[i] == 66870000){
    
    # Read HidroWeb data
    Est <- read.csv2(unz(paste0("data/HidroWeb/",sites[i],"/vazoes_T_",sites[i],".zip"),
                         paste0("vazoes_T_",sites[i],".txt")),
                     skip = 12, sep = ";", row.names=NULL)
    
    # It requires adjustments on the header
    RepNomes <- names(Est)[2:48]
    Est <- Est[1:47]
    names(Est) <- RepNomes
    
    # From wide to long, and other adjustments
    Est <- FlowsHidroWeb(Est)
    
    # Create plot
    png(paste0("plots/daily/",sites[i],".png"),
        height = 500 * 5, width = 800 * 5, res = 75 * 5)
    
    with(Est,
         plot(Data, Vazao,
              type = "l", col = "grey", bty = "n",
              main = paste("Daily flows at",
                           streamGauges$Name[i], "stream gauge"),
              xlab = "Year", ylab = "Discharge (m³/s)",
              xlim = c(as.Date("1900-01-01"),as.Date("2021-01-01"))))
    
    # Summarizes annual minima data
    Est <- Est %>%
      filter(complete.cases(Vazao)) %>%
      mutate(hydroYear = cut(Data, breaks, labels=1900:2023)) %>%
      group_by(hydroYear) %>%
      dplyr::summarize(Minima = min(Vazao, na.rm = TRUE),
                       Dias = length(Vazao),
                       Data = Data[which.min(Vazao)]) %>%
      ungroup() %>%
      filter(Dias >= 210) %>%
      filter(as.numeric(as.character(hydroYear)) < 2022) %>%
      mutate(Site = sites[i])
    
    # Include annual minima data in the daily plots
    with(Est, points(Data, Minima, pch = 20))
    with(Est, lines(Data, Minima, lty = 2))
    
    dev.off()
    
  } else {
    
    # Read HidroWeb data
    Est <- read.csv2(unz(paste0("data/HidroWeb/",sites[i],"/cotas_T_",sites[i],".zip"),
                         paste0("cotas_T_",sites[i],".txt")),
                     skip = 12, sep = ";", row.names=NULL)
    
    # It requires adjustments on the header
    RepNomes <- names(Est)[2:48]
    Est <- Est[1:47]
    names(Est) <- RepNomes
    
    # From wide to long, and other adjustments
    Est <- WaterLevelsHidroWeb(Est)
    
    # Clear inconsistent data manually
    if(sites[i] == 66090000){
      Est[Est$Data > as.Date("2014-07-01") & Est$Data < as.Date("2014-09-01"),]$Cota <- NA
      Est[Est$Data > as.Date("2015-10-31") & Est$Data < as.Date("2015-11-15"),]$Cota <- NA
    } else if(sites[i] == 66125000) {
      Est[Est$Data > as.Date("1977-01-01") & Est$Data < as.Date("1977-02-01"),]$Cota <- NA
    } else if(sites[i] == 66750000) {
      Est[Est$Data > as.Date("1991-05-01") & Est$Data < as.Date("1991-07-01"),]$Cota <- NA
    } else if(sites[i] == 66910000) {
      Est[Est$Data > as.Date("1992-07-01") & Est$Data < as.Date("1992-09-01"),]$Cota <- NA
    }
    
    # Create plot
    png(paste0("plots/daily/",sites[i],".png"),
        height = 500 * 5, width = 800 * 5, res = 75 * 5)
    with(Est,
         plot(Data, Cota,
              type = "l", col = "grey", bty = "n",
              main = paste("Daily water levels at",
                           streamGauges$Name[i], "stream gauge"),
              xlab = "Year", ylab = "Water level (cm)",
              xlim = c(as.Date("1900-01-01"),as.Date("2021-01-01"))))
    
    # Summarizes annual minima data
    Est <- Est %>%
      filter(complete.cases(Cota)) %>%
      mutate(hydroYear = cut(Data, breaks, labels=1900:2023)) %>%
      group_by(hydroYear) %>%
      dplyr::summarize(Minima = min(Cota, na.rm = TRUE),
                       Dias = length(Cota),
                       Data = Data[which.min(Cota)]) %>%
      ungroup() %>%
      filter(complete.cases(hydroYear)) %>%
      filter(Dias >= 210) %>%
      filter(as.numeric(as.character(hydroYear)) < 2021) %>%
      mutate(Site = sites[i])
    
    # Include 2021 data manually in the sites were HidroWeb data is not updated
    # to 2021
    if(sites[i] == 67100000){
      
      Est <- rbind.data.frame(Est,
                              data.frame("hydroYear" = c(2021),
                                         "Minima" = c(77),
                                         "Dias" = c(365),
                                         "Data" = c(as.Date("2021-10-12")),
                                         "Site" = c(67100000)))
      
    } else if(sites[i] == 66825000){
      
      Est <- rbind.data.frame(Est,
                              data.frame("hydroYear" = c(2021),
                                         "Minima" = c(-60),
                                         "Dias" = c(365),
                                         "Data" = c(as.Date("2021-10-15")),
                                         "Site" = c(66825000)))
      
    } else if(sites[i] == 66070004){
      
      Est <- rbind.data.frame(Est,
                              data.frame("hydroYear" = c(2021),
                                         "Minima" = c(26),
                                         "Dias" = c(365),
                                         "Data" = c(as.Date("2021-09-26")),
                                         "Site" = c(66070004)))
      
    } else if(sites[i] == 66260001){
      
      Est <- rbind.data.frame(Est,
                              data.frame("hydroYear" = c(2020,2021),
                                         "Minima" = c(15,4),
                                         "Dias" = c(365,365),
                                         "Data" = c(as.Date("2020-09-30"),as.Date("2021-09-26")),
                                         "Site" = c(66260001,66260001)))
      
    } else if(sites[i] == 66470000){
      
      Est <- rbind.data.frame(Est,
                              data.frame("hydroYear" = c(2021),
                                         "Minima" = c(36),
                                         "Dias" = c(365),
                                         "Data" = c(as.Date("2021-09-29")),
                                         "Site" = c(66470000)))
      
    } else if(sites[i] == 66910000){
      Est <- Est
      Est <- rbind.data.frame(Est,
                              data.frame("hydroYear" = c(2021),
                                         "Minima" = c(124),
                                         "Dias" = c(365),
                                         "Data" = c(as.Date("2021-10-03")),
                                         "Site" = c(66910000)))
      
    } else if(sites[i] == 66945000){
      
      Est <- rbind.data.frame(Est,
                              data.frame("hydroYear" = c(2021),
                                         "Minima" = c(166),
                                         "Dias" = c(365),
                                         "Data" = c(as.Date("2021-09-26")),
                                         "Site" = c(66945000)))
      
    } 
    
    # Include annual minima data in the daily plots
    with(Est, points(Data, Minima, pch = 20))
    with(Est, lines(Data, Minima, lty = 2))
    
    dev.off()
    
  }
  
  df <- rbind(df, Est)

}

df <- df %>%
  mutate(hydroYear = as.numeric(as.character(hydroYear))) %>%
  arrange(hydroYear)

# Useful plot for data visualization
ggplot(df, aes(x = hydroYear, y = Minima)) +
  geom_line() + facet_wrap(~Site, scales = "free_y")

# Add the number length of the times series
siteN <- ddply(df,.(Site), summarise, N=length(Site))
streamGauges <- merge(streamGauges, siteN, by = "Site")

sum(streamGauges$N.y)
# 1373 minimum flows in the study area

# Convert the data from long to wide
df_wide <- df[,-c(3,4)] %>% 
  pivot_wider(names_from = Site, values_from = Minima) %>%
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
  set.seed(2)
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
  xlab("Number of drought regimes") +
  scale_colour_discrete(name = "Information Criteria") +
  scale_color_manual(values=c("#023047"))
p
png(paste0("plots/InfoCriteria.png"),
    height = 500 * 3, width = 800 * 3, res = 75 * 5)
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

# Check the number of stream gauges considered in each year to the multisite HMM
rowSums(!is.na(df_wide_scaled[,-1]))

# Check which class of drought was observed each year in the past
prstates_total <- apply(posterior(fm_total, type = "viterbi")[,c("S1", "S2","S3")], 1, which.max)

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
par(mar = c(4.5,4.5,3.5,7))
plot(df_wide$hydroYear, df_wide$`66825000`, type = "l", bty = "n", col = "grey",
     xlab = "Year", ylab = "Water level (cm)",
     main = "Annual minimum water levels at Ladário stream gauge")

abline(h = seq(-100,400,100), col = "grey", lty = 2)
abline(v = seq(1900,2020,20), col = "grey", lty = 2)

points(df_wide$hydroYear, df_wide$`66825000`, pch = 20, cex = 2.5)
points(df_wide$hydroYear, df_wide$`66825000`, pch = 20, cex = 2,
       col = Paleta[prstates_total])

lines(df_wide$hydroYear, prstates_valores)
points(df_wide$hydroYear, prstates_valores,
       col = Paleta[prstates_total], pch = 20)

add_legend("right",
           legend = paste0("Mean = ", round(Medias$Medias,1), " cm \n",
                           "StdDev = ", round(Medias$Sd,1), " cm \n"),
           pch=20, 
           pt.cex = 2,
           col=Paleta[1:3],
           horiz=FALSE, bty='n', cex=0.8)
dev.off()

# Plot mixed distribution
A <- c(-100: 300)
plot(A, dnorm(A, Medias$Medias[1], Medias$Sd[1]),
     type = "l", bty = "n", ylim = c(0,0.015))
lines(A, dnorm(A, Medias$Medias[2],
               Medias$Sd[2]), type = "l", col = 2)
lines(A, dnorm(A, Medias$Medias[3],
               Medias$Sd[3]), type = "l", col = 3)
Q <- A
fixDistLad <- function(Q){
  D <- sum(prstates_total == 1)/length(prstates_total) *
    dnorm(Q, Medias$Medias[1], Medias$Sd[1]) +
    sum(prstates_total == 2)/length(prstates_total) *
    dnorm(Q, Medias$Medias[2], Medias$Sd[2]) +
    sum(prstates_total == 3)/length(prstates_total) *
    dnorm(Q, Medias$Medias[3], Medias$Sd[3])
  return(D)
}

lines(A, fixDistLad(A), type = "l", col = 4)

# Violin plots of annual minima separated by drought class
df_prstates_total <- data.frame(hydroYear = df_wide$hydroYear,class = prstates_total)
df2 <- merge(df,df_prstates_total, by = "hydroYear", all.x = TRUE)
df2$class <- as.factor(df2$class)

for(i in 1:length(sites)){
  
  if(sites[i] == "66870000"){
    
    p <- ggplot(df2 %>% filter(Site == sites[i]),
                aes(x = as.factor(class), y = Minima, fill = as.factor(class))) +
      geom_violin() + theme_bw() +
      ylab("Discharge (m³/s)") + xlab("Drought classes") +
      ggtitle("Violin plots of annual minimum flows",
              subtitle = paste(streamGauges$Name[i], "stream gauge")) +
      scale_fill_discrete(name = "Drought class")
    
    png(paste0("plots/minimaViolin/",sites[i],".png"),
        height = 500 * 3, width = 800 * 3, res = 75 * 4)
    print(p)
    dev.off()
    
  } else {
    
    p <- ggplot(df2 %>% filter(Site == sites[i]),
                aes(x = as.factor(class), y = Minima, fill = as.factor(class))) +
      geom_violin() + theme_bw() +
      ylab("Water levels (cm)") + xlab("Drought classes") +
      ggtitle("Violin plots of annual minimum water levels",
              subtitle = paste(streamGauges$Name[i], "stream gauge")) +
      scale_fill_discrete(name = "Drought class")
    
    png(paste0("plots/minimaViolin/",sites[i],".png"),
        height = 500 * 3, width = 800 * 3, res = 75 * 4)
    print(p)
    dev.off()
  }
}

# Evaluate if differences between classes are meaningful
Estimates <- list()
p_values <- matrix(NA, nrow = length(sites), ncol = 3)
for(i in 1:length(sites)){
  
  if(sites[i] == 66260002) next
  
  fit_comparison <- lm(Minima ~ as.factor(class),
                       data = df2 %>% filter(Site == sites[i]))
  
  Estimates[[i]] <- fit_comparison %>% summary()
  p_values[i,] <- round(Estimates[[i]]$coefficients[,4], 4)
  
}

names(Estimates) <- sites
Estimates[["66810000"]]

streamGauges$Regimes <- rowSums(ifelse(p_values > 0.05,0,1))
names(streamGauges)[7] <- "N"

df2$Site <- ifelse(df2$Site %in% c("66460000","6470000","66710000",
                                   "66260001","66270000","66340000",
                                   "66900000","66910000","66945000"),
       paste0(df2$Site,"*"), df2$Site)

p <- ggplot(df2 %>% filter(Site != 66260002),
            aes(x = class, y = Minima, fill = class)) +
  geom_violin(colour = "white") + facet_wrap(~Site, scale = "free") +
  theme(panel.background = element_rect(fill = "white")) +
  scale_fill_manual(values = Paleta[1:3],
                    name = "Drought classes") +
  ylab("Annual minima (cm)") +
  xlab("Drought classes") +
  ggtitle("Violin plots of annual minimum water levels in all stream gauges for each drought class")

png(paste0("plots/Regimes.png"),
    height = 500 * 5, width = 800 * 5, res = 75 * 5)
print(p)
dev.off()

streamGauges$Regimes2 <- ifelse(streamGauges$Site %in% c("66460000","6470000","66710000",
                       "66260001","66270000","66340000",
                       "66900000","66910000","66945000"),
       "Regional", "Different")

write.csv(streamGauges, "outputTables/streamGauges.txt", row.names=FALSE)

########################################################################
########################################################################
##
## 3: Evaluate Markov Chains
##
########################################################################
########################################################################

df_prstates_total$hydroYear <- as.numeric(as.character(df_prstates_total$hydroYear))
set.seed(1)
Drought_Regional_Total <- markovchainFit(data = df_prstates_total$class,
                                         method = "bootstrap", nboot = 50)

set.seed(1)
Drought_Regional_1900_1960 <- markovchainFit(data = df_prstates_total$class[df_prstates_total$hydroYear <= 1960],
                                             method = "bootstrap", nboot = 50)

set.seed(1)
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
transitionProbability(Drought_Regional_Total$estimate, "1", "2")
transitionProbability(Drought_Regional_1900_1960$estimate, "1", "2")
transitionProbability(Drought_Regional_1960_2020$estimate, "1", "2")

# Build the transition diagrams
M1 <- round(Drought_Regional_Total$estimate[], 2)
M2 <- round(Drought_Regional_1900_1960$estimate[], 2)
M3 <- round(Drought_Regional_1960_2020$estimate[], 2)

Col <- M1
Col[] <- "darkred"
diag(Col) <- "darkred"

png(paste0("plots/Transition.png"),
    height = 500 * 5, width = 800 * 10, res = 75 * 10)
par(mar = c(3, 1, 2, 1), mfrow = c(1, 3))
names <- paste0("Regime ",c(1,2,3))
plotmat(M1, pos = c(2, 1), name = names, lwd = 1,
        box.lwd = 2, cex.txt = 0.8, box.size = 0.1,
        box.type = "circle", box.prop = 0.5,
        box.col = Paleta[c(1:3)],
        self.lwd = 2,
        self.shiftx = c(-0.1, 0.1, -0.1), self.cex = 0.5,
        arr.col = Col)
mtext("Transition matrix 1900 - 2021")

plotmat(M2, pos = c(2, 1), name = names, lwd = 1,
        box.lwd = 2, cex.txt = 0.8, box.size = 0.1,
        box.type = "circle", box.prop = 0.5,
        box.col = Paleta[c(1:3)],
        self.lwd = 2,
        self.shiftx = c(-0.1, 0.1, -0.1), self.cex = 0.5,
        arr.col = Col)
mtext("Transition matrix 1900 - 1960")

plotmat(M3, pos = c(2, 1), name = names, lwd = 1,
        box.lwd = 2, cex.txt = 0.8, box.size = 0.1,
        box.type = "circle", box.prop = 0.5,
        box.col = Paleta[c(1:3)],
        self.lwd = 2,
        self.shiftx = c(-0.1, 0.1, -0.1), self.cex = 0.5,
        arr.col = Col)
mtext("Transition matrix 1961 - 2021")

par(mfrow = c(1,1))
dev.off()

# Explore other aspects of hidden states
initialState <- c(1,0,0)
after7Years <- initialState * (Drought_Regional_Total$estimate ^ 7)
after7Years

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

fvals(mchain = Drought_Regional_1900_1960$estimate, initialstate = c(0,100,0), n=20)

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

meanFirstPassageTime(Drought_Regional_Total$estimate)
meanRecurrenceTime(Drought_Regional_Total$estimate)

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
## 3: Simulate with the Markov Chains
##
########################################################################
########################################################################

par(mar = c(4.5,4.5,2.5,1))
with(df_prstates_total,
     plot(hydroYear, class, bty = "n", pch = 20, cex = 2.5,
          xlim = c(1970, 2021)))

with(df_prstates_total,
     points(hydroYear, class, pch = 20, col = Paleta[class], cex = 2))
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
my_breaks <- H$breaks
my_colors <- rep(Paleta[3], length(my_breaks))
my_colors[1] <- rgb(0.8,0.2,0.2)

png(paste0("plots/DroughtYearsCount.png"),
    height = 500 * 3, width = 800 * 6, res = 75 * 5)
par(mfrow = c(1,3),
    mar = c(4.5,2.5,1.5,0),
    oma = c(0,2.5,2.5,1))
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

mtext("Histogram of simulated drought-regime years in a 46 years period for each transition matrix",
      outer = TRUE, line = 1, font = 2, cex = 1.2)
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
## 4: Precipitation
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

head(RainfallCHIRPS)
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

png(paste0("plots/RainfallComparison.png"),
    height = 500 * 3, width = 800 * 6, res = 75 * 5)
par(mfrow = c(1,2),
    mar = c(4.5,4.5,2.5,1),
    oma = c(0,0,0,0))
with(Compara,
     plot(Rain_GPCC, Rain_CHIRPS, pch = 20, bty = "n", col = Paleta[Season],
          xlab = "Precipitation with GPCC (mm)",
          ylab = "Precipitation with CHIRPS (mm)",
          main = "Scatter plot of precipiation with CHIRPS against GPCC"))
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




YearRainfall




# Faz o gráfico
Paleta <- c(rgb(0.2,0.2,0.8,0.8),
            rgb(1,0.5,0.2,0.8),
            rgb(0.2,0.8,0.2,0.8),
            rgb(0.8,0.2,0.8,0.8))

Paleta <- Paleta[c(2,3,1,4)]

p1 <- ggplot(YearRainfall, aes(x = Regime, y = Rainfall, group = Regime)) +
  geom_violin(aes(fill = Regime)) +
  theme_bw() +
  ggtitle("Violin plots of annual cummulative precipitation in each drought regime") +
  xlab("Drought classes") +
  scale_fill_manual(values = Paleta[c(1:3)]) +
  theme(plot.title = element_text(size=8))

p2 <- ggplot(YearRainfall, aes(x = Rainfall, group = Regime)) +
  geom_density(aes(fill = Regime), position = "fill") +
  ggtitle("Evolution of Density Kernels of different drought classes with annual precipitation") +
  ylab("Empirical probabilities") +
  xlab("Annual precipitation (mm)") +
  theme_classic() +
  scale_fill_manual(values = Paleta[c(1:3)], guide="none") +
  theme(plot.title = element_text(size=8)) +
  geom_vline(xintercept = 1200) +
  geom_hline(yintercept = 0.87) +
  geom_hline(yintercept = 0.45)
p2

P <- 1350

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
    "Precipitation in different drought regimes in the Paraguay river basin",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

subtitle <- ggdraw() + 
  draw_label(
    "Precipitation estimated in the entire catchment area with CHIRPS (Years >= 1981) and GPCC (Years < 1981)",
    fontface = 'italic',
    size = 10,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )


plot_grid(
  title, subtitle, plot_row,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 0.05, 1)
)


png(paste0("plots/rainfall_Classes.png"),
    height = 500 * 3, width = 800 * 6, res = 75 * 5)
plot_grid(
  title, subtitle, plot_row,
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

Paleta <- Paleta[c(2,3,1,4)]

fit_comparison <- lm(Rainfall ~ as.factor(Regime)-1,
                     data = YearRainfall)
fit_comparison %>% summary()
p_values[i,] <- round(Estimates[[i]]$coefficients[,4], 4)

YearRainfall$class <- prstates_total

png(paste0("plots/RainfallRegimes.png"),
    height = 500 * 5, width = 800 * 5, res = 75 * 5)
par(mar = c(4.5,4.5,3.5,7), mfrow = c(1,1))
plot(YearRainfall$hydroYear, YearRainfall$Rainfall, type = "l", bty = "n", col = "grey",
     xlab = "Year", ylab = "Precipitation (cm)",
     main = "Annual cummulative precipitation over the Paraguay river basin")

abline(h = seq(0,2000,200), col = "grey", lty = 2)
abline(v = seq(1900,2020,20), col = "grey", lty = 2)

points(YearRainfall$hydroYear, YearRainfall$Rainfall, pch = 20, cex = 2.5)
points(YearRainfall$hydroYear, YearRainfall$Rainfall, pch = 20, cex = 2,
       col = Paleta[prstates_total])

lines(YearRainfall$hydroYear, prstates_valores)
points(YearRainfall$hydroYear, prstates_valores,
       col = Paleta[prstates_total], pch = 20)

add_legend("right",
           legend = paste0("Mean = ", round(Medias), " mm \n",
                           "StdDev = ", round(SD), " mm \n"),
           pch=20, 
           pt.cex = 2,
           col=Paleta[1:3],
           horiz=FALSE, bty='n', cex=0.8)
dev.off()
