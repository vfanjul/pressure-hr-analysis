## Set parameters
project = "Heart progerin" # mir29   Heart progerin   Heart lamin  ISO Challenge
technique = "Pressure.HR" # Pressure.HR
daysremoved = 0 # Remove training days
daymedians = F # Calculate medians for each subject per day
outpaired = F # T to exclude outliers in all pressure variables and not only the affected

DPlim = c(30, 140) # Empirical 30 133
SPlim = c(50, 190) # Empirical 51 179
factors = c("Id", "Date")


## Set directories
if (.Platform$OS.type == "unix") setwd("/Volumes/Victor/") else setwd("S:/LAB_VA/LAB/Victor/")
baseroute = paste0(project, " project/Raw data/")
route = paste0(baseroute, technique, "/")

## Functions and libraries
packages = c("beeswarm", "zoo")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) install.packages(setdiff(packages, rownames(installed.packages())))  
for (i in packages) library(i, character.only = T)

datefun = function (x) {
  x = gsub("/", "-", x)
  date1 = sub("-.*", "", x)
  date2 = sub("-.*", "", sub("[[:digit:]]*-", "", x))
  date3 = sub(".*-", "", x)
  if (max(nchar(date1)) == 4) {
    as.Date(x)
  } else if (max(nchar(date3)) == 4) {
    if (max(as.numeric(date2), na.rm = T) > 12 & max(as.numeric(date1), na.rm = T) <= 12) {
      as.Date(paste0(date3, "/", date1, "/", date2))
    } else if (max(as.numeric(date1), na.rm = T) > 12 & max(as.numeric(date2)) <= 12){
      as.Date(paste0(date3, "/", date2, "/", date1))
    } else as.Date(paste0(date3, "/", date2, "/", date1)) # Might mistake days for months
  } else if (max(as.numeric(date1), na.rm = T) > 31 | max(as.numeric(date3), na.rm = T) > as.integer(format(Sys.Date(),"%Y"))-2000) {
    as.Date(paste0(20,x))
  } else if (max(as.numeric(date3), na.rm = T) > 31 | max(as.numeric(date1), na.rm = T) > as.integer(format(Sys.Date(),"%Y"))-2000) {
    if (max(as.numeric(date2), na.rm = T) > 12 & max(as.numeric(date1), na.rm = T) <= 12) {
      as.Date(paste0(20, date3, "/", date1, "/", date2))
    } else if (max(as.numeric(date1), na.rm = T) > 12 & max(as.numeric(date2), na.rm = T) <= 12){
      as.Date(paste0(20, date3, "/", date2, "/", date1))
    } else as.Date(paste0(20, date3, "/", date2, "/", date1)) # Might mistake days for months
  } else if (max(as.numeric(date2), na.rm = T) > 12) {
    as.Date(paste0(20, date3, "/", date1, "/", date2))
  } else as.Date(paste0(20, date3, "/", date2, "/", date1)) # Might mistake between days, months and years
}


## Import and merge files
rawdata = read.delim(paste0(route, technique, ".data.txt"), encoding = "latin1", row.names = NULL, stringsAsFactors = F, na.strings = c("N/A", "NA", "#N/A", "", "nan"))
ylabels = as.character(read.delim(paste0(route, technique, ".data.txt"), encoding = "latin1", header = F, row.names = NULL, stringsAsFactors = F, dec = ".", na.strings = "N/A")[1,])
ylabels = setdiff(ylabels, factors)
names(rawdata)[grep("ystolic", names(rawdata))] = "SP"
names(rawdata)[grep("iastolic", names(rawdata))] = "DP"
names(rawdata)[grep("eart", names(rawdata))] = "HR"
vars = setdiff(names(rawdata), factors)
rawdata = rawdata[sapply(1:nrow(rawdata), function (x) !all(is.na(rawdata[x,vars]))),]
for (var in vars) if (nrow(rawdata[!is.na(rawdata[,var]),]) == 0) {
  ylabels = ylabels[which(vars != var)]
  rawdata = rawdata[,-(length(factors) + which(vars == var))]
  vars = vars[vars != var]
}
rawdata$Id = toupper(rawdata$Id)
rawdata$Date = datefun(rawdata$Date)


## Remove global outliers
if (length(intersect(vars, c("SP", "DP"))) > 1) for (pres in c("SP", "DP")) {
  if (outpaired) colpress = c("SP", "DP") else colpress = pres
  rawdata[!is.na(rawdata[,pres]) & (rawdata[,pres] > get(paste0(pres, "lim"))[2] | rawdata[,pres] < get(paste0(pres, "lim"))[1]), colpress] = NA
}


## Remove outliers within subjects by day
for (pres in rev(vars)) {
  bxp = boxplot(get(pres) ~ Id + Date, rawdata, plot = F)
  outs = cbind.data.frame("Id" = gsub("\\.[[:graph:]]*", "", bxp$names[bxp$group]), "Date" = gsub("[[:graph:]]*\\.", "", bxp$names[bxp$group]), bxp$out, 1)
  names(outs) [3:4] = c(pres, paste0("out.", pres))
  rawdata = merge(rawdata, outs, by = names(outs)[-4], all = T)
  rawdata[!is.na(rawdata[,names(outs)[4]]), names(outs)[3]] = NA
}
if (outpaired) for (pres in c("SP", "DP")) rawdata[!is.na(rawdata[,paste0("out.", pres)]), c("SP", "DP")] = NA
rawdata = rawdata[-grep("out.", names(rawdata))]


## Calculate other variables
if (length(intersect(vars, c("SP", "DP"))) > 1) {
  rawdata$PP = rawdata$SP - rawdata$DP
  rawdata$MAP = rawdata$DP + rawdata$PP/3
  ylabels = c(ylabels, "Pulse pressure (mmHg)", "Mean arterial pressure (mmHg)")
  vars = setdiff(names(rawdata), factors)
}

## Days removed
if (daysremoved > 0) {
  followups = c()
  for (i in (unique(rawdata$Id))) followups = rbind.data.frame(followups, data.frame("Id" = i, "Date" = sort(unique(rawdata$Date[rawdata$Id == i])), "Followup" = 1:length(unique(rawdata$Date[rawdata$Id == i]))))
  rawdata = merge(rawdata, followups, by = factors, all = T)
  rawdata = rawdata[rawdata$Followup > daysremoved,]
  rawdata = rawdata[,-ncol(rawdata)]
}


## Day medians
daymeandata = rawdata[!duplicated(rawdata[,c("Id", "Date")]),]
rawdata$Dateid = paste0(rawdata$Id, rawdata$Date)
if (length(vars) > 1) {
  daymeandata[,vars] = t(sapply(unique(rawdata$Dateid), function (x) colMeans(rawdata[rawdata$Dateid == x,vars], na.rm = T)))
} else  daymeandata[,vars] = sapply(unique(rawdata$Dateid), function (x) colMeans(data.frame(rawdata[rawdata$Dateid == x,vars]), na.rm = T))
rawdata = rawdata[,-ncol(rawdata)]
if (daymedians) rawdata = daymeandata


## Export files
plots_per_col = ifelse(length(unique(rawdata$Id)) == 3, 1, round(sqrt(length(unique(rawdata$Id))),0))
plots_per_row = ifelse(length(unique(rawdata$Id)) == 3, 3, ceiling(length(unique(rawdata$Id))/plots_per_col))
if (length(intersect(vars, c("SP", "DP"))) > 1) {
  correlations = list(c("SP", "DP"), c("HR", "MAP"))
  for (corr in correlations) {
    pdf(paste0(route, paste0(corr, collapse = " "), " correlation.pdf"), plots_per_row*2, plots_per_col*2, pointsize = 12)
    par(mfrow = c(plots_per_col,plots_per_row), bty = "l", pch = 19, cex.axis = 0.7, cex.main = 0.7, mar = c(2,2,0.5,0.5), mgp = c(2,1,0), oma = c(1,1,1,1))
    for (i in unique(rawdata$Id)) {
      plot(get(corr[1]) ~ get(corr[2]), rawdata[rawdata$Id == i,], ylim = summary(rawdata[,corr[1]])[c(1,6)], xlim = summary(rawdata[,corr[2]])[c(1,6)], col = rainbow(nrow(rawdata[rawdata$Id == i,]), start = 0, end = 2/6))
      abline(h = boxplot(rawdata[rawdata$Id == i,corr[1]], plot = F, na.rm = T)$stats, col = 8, lty = c(3,2,1,2,3))
      abline(v = boxplot(rawdata[rawdata$Id == i,corr[2]], plot = F, na.rm = T)$stats, col = 8, lty = c(3,2,1,2,3))
      abline(h = boxplot(rawdata[,corr[1]], plot = F, na.rm = T)$stats[c(1,5)], col = 1, lty = 3)
      abline(v = boxplot(rawdata[,corr[2]], plot = F, na.rm = T)$stats[c(1,5)], col = 1, lty = 3)
      legend("top", legend = i, bty = "n")
    }
    title(paste(ylabels[which(vars == corr[1])], "vs", ylabels[which(vars == corr[2])]), outer = T, cex.main = 1)
    dev.off()
  }
}

for (var in vars) {
  pdf(paste0(route, var, " plot.pdf"), plots_per_row*2, plots_per_col*2, pointsize = 12)
  par(mfrow = c(plots_per_col,plots_per_row), bty = "l", pch = 19, cex.axis = 0.7, cex.main = 0.7, mar = c(2,2,0.5,0.5), mgp = c(2,1,0), oma = c(1,1,1,1))
  for (i in unique(rawdata$Id)) {
    boxplot(get(var) ~ Date, rawdata[rawdata$Id == i,], ylim = summary(rawdata[,var])[c(1,6)], add = F, border = rainbow(length(unique(rawdata$Date[rawdata$Id == i]))), boxlwd = 1)
    beeswarm(get(var) ~ Date, rawdata[rawdata$Id == i,], ylim = summary(rawdata[,var])[c(1,6)], add = T, col = adjustcolor(rainbow(length(unique(rawdata$Date[rawdata$Id == i]))), alpha.f = 0.2), corral = "wrap")
    abline(h = boxplot(rawdata[rawdata$Id == i, var], plot = F, na.rm = T)$stats, col = 8, lty = c(3,2,1,2,3))
    abline(h = boxplot(rawdata[,var], plot = F, na.rm = T)$stats[c(1,5)], col = 1, lty = 3)
    legend("top", legend = i, bty = "n")
  }
  title(ylabels[which(vars == var)], outer = T, cex.main = 1)
  dev.off()
}
graphics.off()

names(rawdata) = c(factors, ylabels)
names(daymeandata) = c(factors, ylabels)
write.table(rawdata, file = paste0(route, "rawdata.txt"), row.names = F, col.names = T, sep = "\t", append = F)
if (!daymedians) write.table(daymeandata, file = paste0(route, "daymeandata.txt"), row.names = F, col.names = T, sep = "\t", append = F)
