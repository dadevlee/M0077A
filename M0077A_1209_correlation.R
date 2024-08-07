library(astrochron)
library(readxl)
library(viridis)
library(quantmod)
library(jpeg)

setwd("/Users/daf_d/OneDrive - Universität Münster/Documents/Science/Exp 364/Paleocene/")
rm(list = ls())
my_colorbar = viridis(n = 3, direction = 1)
#source("/Users/pdoc3/Documents/Cenozoic_stack/MegaSplice2/Nat_Geo_GP/20160408_standardizationchange/filled.contour2.R")

La11 = getLaskar(sol = "la11", verbose = F)
La11 = La11[c(50000:70000),]
Cenogrid = read_excel("Cenogrid.xlsx")
Cenogrid = Cenogrid[,c(5,8,9)]
#write.csv(La11,"QAS/La11.csv", row.names = F)
#La11 = read.csv("QAS/La11.csv")

#####################
# Select the L* data of the studied interval
##################### ----
M0077A_Lstar=read.csv("364-M0077A_CRF.tab.csv")
M0077A_Lstar = M0077A_Lstar[,c(2,6)]
idx = which(M0077A_Lstar[,1]>= 607 & M0077A_Lstar[,1]<= 617)
idx2 = which(M0077A_Lstar[,1]>= 607.3 & M0077A_Lstar[,1]<= 616.5)
M0077A_Lstar_Paleocene = M0077A_Lstar[idx,]
M0077A_Lstar_cyclo = M0077A_Lstar[idx2,]
write.csv(M0077A_Lstar_cyclo,"QAS/M007A_Lstar_cyclo.csv",row.names = F)

dev.off()
plot(M0077A_Lstar_Paleocene, type = "l", lwd = 1.5)
lines(M0077A_Lstar_cyclo, col = "orange", lwd = 1.5)

#######################
# MTM spectral analysis in the depth domain
#######################
M0077A_Lstar_cyclo = linterp(M0077A_Lstar_cyclo)
M0077A_Lstar_cyclo = bandpass(M0077A_Lstar_cyclo, flow = 0.75,fhigh = 35, xmax = 2)
mtm_depth = mtm(M0077A_Lstar_cyclo, ntap = 3, output = 1, detrend = T, genplot = F)

dev.off()
plot(mtm_depth$Frequency, mtm_depth$Power, type = "l", xlim = c(0,25), xaxs = "i",
     ylim = c(0,0.2), yaxs = "i", axes = F, xlab = "", ylab = "")
lines(mtm_depth$Frequency, mtm_depth$AR1_95_power, col = "gray50", lwd = 2)
rect(1,5e-6,1.5,5e-1, col = "orange", border = NA)
rect(4,5e-6,5.6,5e-1, col = "orange", border = NA)
lines(mtm_depth$Frequency, mtm_depth$Power, lwd = 1.5)
lines(mtm_depth$Frequency, mtm_depth$AR1_95_power, col = "gray50", lwd = 2)
text(2, 0.21, '405-kyr', cex = 1.5, xpd = NA, srt = 45)
text(5, 0.21, '100-kyr', cex = 1.5, xpd = NA, srt = 45)
text(12, 0.11, 'Obliquity?', cex = 1.5, xpd = NA, srt = 45)
axis(1, tck = -0.01)
axis(2, tck = 0.01, at = c(0,0.2))
mtext("Spectral Power", side = 2, cex = 1.5, line = 2)
mtext("Frequency (cycles/m)", side = 1, cex = 1.5, line = 2)
#text(17, 0.11, 'Precession?', cex = 1.5, xpd = NA, srt = 45)
#text(1, 0.21, 'ecc.', cex = 1.5, xpd = NA)

#################
# Bandpass filtering eccentricity in the depth domain
#################
M007A_Lastar_cyclo_405 = bandpass(M0077A_Lstar_cyclo, flow = 1, fhigh = 1.5, xmax = 10)
M007A_Lastar_cyclo_100 = bandpass(M0077A_Lstar_cyclo, flow = 4, fhigh = 5.6, xmax = 10)
M007A_Lastar_cyclo_ecc = M007A_Lastar_cyclo_405
M007A_Lastar_cyclo_ecc[,2] = M007A_Lastar_cyclo_405[,2]+M007A_Lastar_cyclo_100[,2]

###################
# Tuning was performed manually in QAnalySeries
# Here, we read in the tie points (Table 1 in paper)
###################
ties = read.delim("QAS/QAS6.txt", sep = "", header = F)
ties[,3] = (ties[,1]-607)/10*(66.2-61.3)+61.3
# Age-depth plot 
plot(ties[,2], ties[,1], pch = 19, ylim = c(618,607), col = "darkorange")
lines(ties[,2], ties[,1], col = "darkorange", lwd = 2)

####################
# Tuning plot (Figure 6 in paper)
####################
dev.off()
layout(matrix(c(1,2,3,4),4,1, byrow = T), heights = c(1,0.5,0.25,1.25))
par(mar = c(0,5,5,2))
plot(Cenogrid$`Tuned time [Ma]`, Cenogrid$`Foram bent Î´13C [â€° PDB] (VPDB)`,
     xlim = c(61.3,66.2), type = "l", ylim = c(2.4,0.4),
     axes = F, xlab = "", ylab = "",
     xaxs = "i")
abline(v = ties[,2]/1000, col = "grey50", lwd = 2)
axis(2, at = c(0.5,2.0), cex.axis = 1.5,tck = 0.01)
mtext(expression(paste(delta^13,"C"[benthic])), side = 2, cex = 1.4, line = 2.5)
axis(3, at = seq(61.5,66,0.5), labels = seq(61.5,66,0.5), cex.axis = 1.5, tck = 0.01)
mtext("Time (Ma)", side = 3, cex = 1.4, line = 2.5)

par(mar = c(0,5,0,2))
plot(La11$Time_ka, La11$ecc_LA11, xlim = c(61300,66200),type = "l", 
     axes = F, xlab = "", ylab = "",
     xaxs = "i")
axis(2, at = c(0,0.05), cex.axis = 1.5,tck = 0.01)
mtext("Eccentricity", side = 2, cex = 1.4, line = 2.5)
abline(v = ties[,2], col = "grey50", lwd = 2)

par(mar = c(0,5,0,2))
plot(c(0), xlim = c(61.3,66.2), xaxs = "i", ylim = c(0,1), axes = F, xlab = "", ylab = "")
segments(ties[,2]/1000,1,ties[,3],0, col = "grey50", lwd = 2)

par(mar = c(5,5,0,2))
plot(M0077A_Lstar_cyclo, type = "l", 
     xlim = c(607,617), ylim = c(70,50),
     axes = F, xlab = "", ylab = "", xaxs = "i")
axis(1, at = seq(607,617,1), cex.axis = 1.5, tck = 0.01)
axis(2, at = c(50,70), cex.axis = 1.5,tck = 0.01)
mtext("Depth (mbsf)", side = 1, cex = 1.6, line = 2.5)
mtext("L*", side = 2, cex = 1.6, line = 1)
par(new = T)
plot(M007A_Lastar_cyclo_ecc, col = "red", type = "l", lwd = 2,
     xlim = c(607,617), ylim = c(130,110),
     axes = F, xlab = "", ylab = "", xaxs = "i")
abline(v = ties[,1], col = "grey50", lwd = 2)

################
# Figure 7: Comparison with CENOGRID ----
################
M0077A_Lstar_tuned = tune(M0077A_Lstar,ties[,c(1,2)])
M0077A_Lstar_tuned_i = tune(M0077A_Lstar_cyclo,ties[,c(1,2)])
M0077A_Lstar_tuned_i = linterp(M0077A_Lstar_tuned_i, dt = 6)
MTM_tuned = mtm(M0077A_Lstar_tuned_i, ntap = 3, output = 1)

idx = which(Cenogrid$`Tuned time [Ma]` > 61.3 & Cenogrid$`Tuned time [Ma]`<66)
Cenogrid_i = linterp(Cenogrid[idx,c(1,2)], dt = 0.005)
MTM_Cenogrid = mtm(Cenogrid_i, ntap = 3, output = 1)

dev.off()
layout(matrix(c(1,2),1,2, byrow = T), heights = c(1), widths = c(1,2))
par(mar = c(5,5,2,2))
plot(MTM_tuned$Frequency, MTM_tuned$Power, xlim = c(0,0.05), ylim = c(0,0.2), type = "l", 
     col = "darkorange", xaxs = "i", yaxs = "i", axes = F, 
     xlab = "", ylab = "")
axis(1, at = c(0,0.01,0.025,0.05), labels = c("0", "0.01", "0.025", "0.05"), tck = -0.01)
axis(2, at = c(0,0.2), tck = 0.01)
mtext("Spectral Power", side = 2, cex = 1.4, line = 2)
mtext("Frequency (cycles/kyr)", side = 1, cex = 1.4, line = 2.5)
text(0.004, 0.19, '405-kyr', cex = 1.5, xpd = NA, srt = 45)
text(0.015, 0.14, '100-kyr', cex = 1.5, xpd = NA, srt = 45)
text(0.03, 0.07, 'Obliquity?', cex = 1.5, xpd = NA, srt = 45)
text(0.04,0.19,"A", cex = 2)


#par(new = T)
#plot(MTM_Cenogrid$Frequency/1000, MTM_Cenogrid$Power, xlim = c(0,0.05), type = "l", log = "y")
par(mar = c(5,5,2,5))
plot(Cenogrid$`Tuned time [Ma]`, Cenogrid$`Foram bent Î´13C [â€° PDB] (VPDB)`, 
     type = "l", xlim = c(61.3,66.000), xaxs = "i", ylim = c(2.5,0),
     axes = F, xlab = "", ylab = "")
axis(1)
axis(2, at = c(0.5, 2))
mtext(expression(paste(delta^13,"C"[benthic])), side = 2, cex = 1.4, line = 2.5)
mtext("Time (Ma)", side = 1, cex = 1.4, line = 2.5)
par(new = T)
plot(M0077A_Lstar_tuned, type = "l", col = "darkorange", ylim = c(80,40), xaxs = "i", 
     yaxs = "i", xlim = c(61300,66000), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
rect(59200,80,61700,78, col = rgb(253/255,180/255,98/255))
rect(61700,80,66000,78, col = rgb(254/255,151/255,101/255))
axis(4, at= c(40,60,80))
mtext("L*", side = 4, cex = 1.4, line = 2.5)
text(63500,79, "Danian", cex = 1.4)
text(61500,79, "Sel.", cex = 1.4)
text(61500,42,"B", cex = 2)

###############################
# Figure 8
##############################
La04 = getLaskar(sol = "la04")
La04 = La04[50000:70000,c(1,3)]
OblAM = hilbert(La04, demean = T, addmean = T)
OblAM = lowpass(OblAM, fcut = 1/300, xmax = 1/100)
idx = findValleys(OblAM[,2])
OblAM_ages = OblAM[idx,1]

SB_depths = matrix(c(616.58, 614.17, 611.64, 609.64, 607.74, 0,0,0,0,0), ncol = 2, nrow = 5)
SB_ages = tune(SB_depths, ties[,c(1,2)], extrapolate = T)

#read jpgs
section_depths = read.delim("Paper/Fig._8/364-M0077A_core_section.txt", skip = 18)
section_top_depths = cbind(section_depths[98:107,4],section_depths[98:107,4])
section_bottem_depths = cbind(section_depths[98:107,5],section_depths[98:107,5])
section_top_ages = tune(section_top_depths, ties[,c(1,2)], extrapolate = T)
section_bottem_ages = tune(section_bottem_depths, ties[,c(1,2)], extrapolate = T)

C37R1 = readJPEG("Paper/Fig._8/364-M0077A-037R-001_11_scan (1).jpg")
C37R2 = readJPEG("Paper/Fig._8/364-M0077A-037R-002_11_scan.jpg")
C37R3 = readJPEG("Paper/Fig._8/364-M0077A-037R-003_11_scan.jpg")
C38R1 = readJPEG("Paper/Fig._8/364-M0077A-038R-001_11_scan.jpg")
C38R2 = readJPEG("Paper/Fig._8/364-M0077A-038R-002_11_scan.jpg")
C38R3 = readJPEG("Paper/Fig._8/364-M0077A-038R-003_11_scan.jpg")
C39R1 = readJPEG("Paper/Fig._8/364-M0077A-039R-001_11_scan.jpg")
C39R1top = readJPEG("Paper/Fig._8/364-M0077A-039R-001_11_scan (1).jpg")
C39R1bottem = readJPEG("Paper/Fig._8/364-M0077A-039R-001_11_scan (2).jpg")
C39R2 = readJPEG("Paper/Fig._8/364-M0077A-039R-002_11_scan.jpg")
C39R3 = readJPEG("Paper/Fig._8/364-M0077A-039R-003_11_scan.jpg")
C40R1 = readJPEG("Paper/Fig._8/364-M0077A-040R-001_95_scan (1).jpg")

C37R1 <- aperm(C37R1, c(1, 2, 3))[nrow(C37R1):1, , ]
C37R2 <- aperm(C37R2, c(2, 1, 3))[ , 1:nrow(C37R2), ]
C37R3 <- aperm(C37R3, c(2, 1, 3))[ , 1:nrow(C37R3), ]
C38R1 <- aperm(C38R1, c(2, 1, 3))[ , 1:nrow(C38R1), ]
C38R2 <- aperm(C38R2, c(2, 1, 3))[ , 1:nrow(C38R2), ]
C38R3 <- aperm(C38R3, c(2, 1, 3))[ , 1:nrow(C38R3), ]
#C39R1 <- aperm(C39R1, c(2, 1, 3))[ , 1:nrow(C39R1), ]
C39R1top <- aperm(C39R1top, c(1, 2, 3))[nrow(C39R1top):1, , ]
C39R1bottem <- aperm(C39R1bottem, c(1, 2, 3))[nrow(C39R1bottem):1, , ]
C39R2 <- aperm(C39R2, c(2, 1, 3))[ , 1:nrow(C39R2), ]
C39R3 <- aperm(C39R3, c(2, 1, 3))[ , 1:nrow(C39R3), ]
C40R1 <- aperm(C40R1, c(1, 2, 3))[nrow(C40R1):1, ,]

# Read Sea Level reconstructions by other authors ----
SL_reconstructions = read_excel("Paper/Fig._8/miller-sealev-originalcontrib.xls", sheet = "Data")
SL_Haq1987 = cbind(SL_reconstructions$`Age (GTS12; Ma)
Haq et al. (1987)`, SL_reconstructions$`Sea level
Haq et al. (1987)
0-244 Ma`)
SL_Komminz2008 = cbind(SL_reconstructions$`Age (GTS12; Ma)
sea level estimated  (Kominz et al., 2008)`,SL_reconstructions$`Sea level (m) backstripped w/ estimated lowstand  (Kominz et al., 2008)`)
SL_Miller2005 = cbind(SL_reconstructions$`Age (GTS12; Ma)
sea level estimated  (Miller et al., 2005)`, SL_reconstructions$`Sea level (m) backstripped w/ estimated lowstand  (Miller et al., 2005)`)
SL_Miller2020 = read.delim("Paper/Fig._8/Cenozoic_sea_level_reconstruction_smoothed.tab", skip = 33)



# Start plotting ----
pdf("my_plot.pdf", width = 12, height = 9)
layout(matrix(c(1, 2), ncol = 1), heights = c(1, 3))
par(mar = c(1,5,1,1))
plot(OblAM[,1]/1000, OblAM[,2], xlim = c(61,66), type = "l", ylab = "", xlab = "", axes = F)
axis(1)
axis(2)
mtext("1.2-Myr Obliquity (rad)", side = 2, line = 2.5)
box()

par(mar = c(5,5,1,1))
plot(SL_Miller2020[,1]/1000, SL_Miller2020[,2], xlim = c(61,66), ylim = c(-20,100), 
     type = "l", yaxs = "i", xlab = "Time (Ma)", ylab = "Sea Level (m above present) ", col = '#228B22')
#lines(SL_Haq1987)
lines(SL_Komminz2008, col = '#DC143C')
lines(SL_Miller2005, col = '#4169E1')
rasterImage(C37R1,xleft = 61.53, xright = section_bottem_ages[1,1]/1000, ybottom = -22, ytop = 0)
rasterImage(C37R2,xleft = section_top_ages[2,1]/1000, xright = section_bottem_ages[2,1]/1000, ybottom = -22, ytop = 0)
rasterImage(C37R3,xleft = section_top_ages[3,1]/1000, xright = section_bottem_ages[3,1]/1000, ybottom = -22, ytop = 0)
rasterImage(C38R1,xleft = section_top_ages[4,1]/1000, xright = section_bottem_ages[4,1]/1000, ybottom = -22, ytop = 0)
rasterImage(C38R2,xleft = section_top_ages[5,1]/1000, xright = section_bottem_ages[5,1]/1000, ybottom = -22, ytop = 0)
rasterImage(C38R3,xleft = section_top_ages[6,1]/1000, xright = section_bottem_ages[6,1]/1000, ybottom = -22, ytop = 0)
#rasterImage(C39R1,xleft = section_top_ages[7,1]/1000, xright = section_bottem_ages[7,1]/1000, ybottom = -22, ytop = 0)
rasterImage(C39R1top,xleft = section_top_ages[7,1]/1000, xright = SB_ages[4,1]/1000, ybottom = -22, ytop = 0)
rasterImage(C39R1bottem,xleft = SB_ages[4,1]/1000, xright = section_bottem_ages[7,1]/1000, ybottom = -22, ytop = 0)
rasterImage(C39R2,xleft = section_top_ages[8,1]/1000, xright = section_bottem_ages[8,1]/1000, ybottom = -22, ytop = 0)
rasterImage(C39R3,xleft = section_top_ages[9,1]/1000, xright = section_bottem_ages[9,1]/1000, ybottom = -22, ytop = 0)
rasterImage(C40R1,xleft = section_top_ages[10,1]/1000, xright = 66.035, ybottom = -22, ytop = 0)

rect(SB_ages[1,1]/1000, 0, SB_ages[2,1]/1000, 10)
rect(SB_ages[2,1]/1000, 0, SB_ages[3,1]/1000, 10)
rect(SB_ages[3,1]/1000, 0, SB_ages[4,1]/1000, 10)
rect(SB_ages[4,1]/1000, 0, SB_ages[5,1]/1000, 10)
rect(67,0,60,10, col = rgb(253/255, 167/255, 95/255))
rect(67,0,66.035,10, col = rgb(166/255, 216/255, 74/255))
#rect(61.35,-20,59,0, col = "white")
abline(h = 10)
abline(h = 0)
text(61.2,5,"Seq. 5", cex = 2)
text(62,5,"Seq. 4", cex = 2)
text(63,5,"Seq. 3", cex = 2)
text(64.25,5,"Seq. 2", cex = 2)
text(65.5,5,"Seq. 1", cex = 2)
dev.off()

# Other stuff... !!! Not used in the paper!!! ----
idx = findValleys(M007A_Lastar_cyclo_100[,2])
depths = M007A_Lastar_cyclo_100[idx,1]
ages = seq()

age_model = matrix(c(609.7,611.7,614.2,616.6,62.35,63.38,64.72,66), nrow = 4, ncol = 2)
write.csv(age_model,"QAS/age_model.csv", row.names = F)

M0077A_Lstar_biostrat = tune(M0077A_Lstar, age_model, extrapolate = F)

dev.off()
plot(M0077A_Lstar_biostrat, type = "l")
M0077A_Lstar_biostrat = linterp(M0077A_Lstar_biostrat, dt = 0.01)

mtm_biostrat = mtm(M0077A_Lstar_biostrat, ntap = 3, output = 1)
dev.off()
plot(mtm_biostrat$Frequency, mtm_biostrat$Power, type = "l", log = "y")

bandpass(M0077A_Lstar_biostrat, flow = 1/0.150, fhigh = 1/0.080)

Fe_1209 = read.csv("1209_Fe.csv")
Fe_1209 = Fe_1209[,c(5,4)]
Fe_1209 = sortNave(Fe_1209)

Cenogrid = read_excel("Cenogrid.xlsx")
Cenogrid = Cenogrid[,c(5,8,9)]



dev.off()
par(mar = c(5,5,5,5))
plot(Fe_1209, log = "y", type = "l", xlim = c(61700,66000), xaxs = "i")
par(new = T)
plot(M0077A_Lstar, type = "l", col = "red", ylim = c(80,30), xaxs = "i", xlim = c(607.8,616.5), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
axis(3)
mtext("Depth M0077A (m)", 3, line = 2)
axis(4)
mtext("L*", 4, line = 2)

dev.off()
par(mar = c(5,5,5,5))
plot(La11, type = "l", xlim = c(59200,66000), xaxs = "i")
par(new = T)
plot(M0077A_Lstar, type = "l", col = "red", ylim = c(80,30), xaxs = "i", xlim = c(607.5,617.0), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
axis(3)
mtext("Depth M0077A (m)", 3, line = 2)
axis(4)
mtext("L*", 4, line = 2)

dev.off()
par(mar = c(5,5,5,5))
plot(Cenogrid$`Tuned time [Ma]`, Cenogrid$`Foram bent Î´13C [â€° PDB] (VPDB)`, type = "l", xlim = c(62,66.000), xaxs = "i", ylim = c(2.5,0))
par(new = T)
plot(M0077A_Lstar, type = "l", col = "red", ylim = c(80,30), xaxs = "i", xlim = c(607.8,616.5), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
axis(3)
mtext("Depth M0077A (m)", 3, line = 2)
axis(4)
mtext("L*", 4, line = 2)

dev.off()
par(mar = c(5,5,5,5))
plot(Cenogrid$`Tuned time [Ma]`, Cenogrid$`Foram bent Î´18O [â€° PDB] (VPDB)`, type = "l", xlim = c(59.200,66.000), xaxs = "i", ylim = c(1,-1))
par(new = T)
plot(M0077A_Lstar, type = "l", col = "red", ylim = c(80,30), xaxs = "i", xlim = c(607.5,617.0), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
axis(3)
mtext("Depth M0077A (m)", 3, line = 2)
axis(4)
mtext("L*", 4, line = 2)

# Let's try to look at what it looks like in the time domain.
agemodel = matrix(c(607.5,617,59200,66000), nrow = 2, ncol = 2)
M0077A_Lstar_time = tune(M0077A_Lstar, agemodel, extrapolate = F)
M0077A_Lstar_time = linterp(M0077A_Lstar_time)
M0077A_Lstar_time = detrend(M0077A_Lstar_time)

mtm_M0077A_Lstar = mtm(M0077A_Lstar_time, ntap = 3, output = 1)

# Calculate sedimentation rates
sedrate = c()
for (i in 1:13){sedrate[i] = (ties[i+1,1]-ties[i,1])*100/(ties[i+1,2]-ties[i,2])}
