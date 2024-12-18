nino <- "https://psl.noaa.gov/gcos_wgsp/Timeseries/Data/nino34.long.anom.data"
nino <- read.table(nino, skip = 1, nrows = 155, col.names = c("Year", month.abb), na.strings = "-99.99") %>%
  dplyr::select(Year, Dec)
nino$HydroYear <- nino$Year + 1

pdo <- "https://www.ncei.noaa.gov/pub/data/cmb/ersst/v5/index/ersst.v5.pdo.dat"
pdo <- read.table(pdo, skip = 2, nrows = 170, col.names = c("Year", month.abb), na.strings = "-99.99") %>%
  dplyr::select(Year, Dec)
pdo$HydroYear <- pdo$Year + 1
tail(pdo)

amo <- "https://psl.noaa.gov/data/correlation/amon.us.long.data"
amo <- read.table(amo, skip = 1, nrows = 167, col.names = c("Year", month.abb), na.strings = "-99.99") %>%
  dplyr::select(Year, Dec)
amo$HydroYear <- amo$Year + 1


par(mfrow = c(3,1))
with(nino, plot(Year, Dec, bty = "n", pch = 20, xlim = c(1900, 2030)))
abline(h = 0, lty = 2, col = "grey")
abline(v = seq(1800,2100,10), lty = 2, col = "grey")

fit <- loess(Dec ~ Year, nino, span = 0.2)
lines(nino$Year, predict(fit, nino), col = "orange")
mtext("ENSO 3.4", font = 2)

with(pdo, plot(Year, Dec, bty = "n", pch = 20, xlim = c(1900, 2030)))
abline(h = 0, lty = 2, col = "grey")
abline(v = seq(1800,2100,10), lty = 2, col = "grey")
fit <- loess(Dec ~ Year, pdo, span = 0.2)
lines(pdo$Year, predict(fit, pdo), col = "orange")
mtext("PDO", font = 2)

with(amo, plot(Year, Dec, bty = "n", pch = 20, xlim = c(1900, 2030)))
abline(h = 0, lty = 2, col = "grey")
abline(v = seq(1800,2100,10), lty = 2, col = "grey")
fit <- loess(Dec ~ Year, amo, span = 0.5)
lines(amo$Year, predict(fit, amo), col = "orange")
mtext("AMO", font = 2)
