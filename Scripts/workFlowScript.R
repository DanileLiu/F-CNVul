# 00. Package Loading ------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(readxl)
library(topsis)
library(mapdata)
library(maps)
library(maptools)
library(sf)
library(ggmap)


# 01. Functions -----------------------------------------------------------
## 定义归一化函数

Rescale = function(x, type=1) {
    # type=1正向指标, type=2负向指标
    rng = range(x, na.rm = TRUE)
    if (type == 1) {
        (0.996 - 0.002) * (x - rng[1]) / (rng[2] - rng[1]) + 0.002
    } else {
        (0.996 - 0.002) *  (rng[2] - x) / (rng[2] - rng[1]) + 0.002
    }
}

Entropy_Weight = function(X, index) {
    # 实现用熵权法计算各指标(列）的权重及各数据行的得分
    # X为指标数据, 一行代表一个样本, 每列对应一个指标
    # index指示向量，指示各列正向指标还是负向指标，1表示正向指标，2表示负向指标
    # s返回各行（样本）得分，w返回各列权重
    
    pos = which(index == 1)
    neg = which(index != 1)
    
    # 数据归一化
    X[,pos] = lapply(X[,pos], Rescale, type=1)
    X[,neg] = lapply(X[,neg], Rescale, type=2)
    
    # 计算第j个指标下，第i个样本占该指标的比重p(i,j)           
    P = data.frame(lapply(X, function(x) x / sum(x)))
    
    # 计算第j个指标的熵值e(j)
    e = sapply(P, function(x) sum(x * log(x)) *(-1/log(nrow(P))))
    
    d = 1 - e         # 计算信息熵冗余度
    w = d / sum(d)   # 计算权重向量
    
    # 计算样本得分
    s = as.vector(100 * as.matrix(X) %*% w)
    
    list(w=w, s=s)
}

plotTriAxes <- function(indusData) {
    propTab <- prop.table(x = as.matrix(indusData[-1]),
                          margin = 1)
    
    coordTab <- data.frame(year = indusData[,1],
                           x = (propTab[,1] - propTab[,3])/3,
                           y = (propTab[,2] - propTab[,3])/3)
    
    par(xaxs='i',yaxs='i',mar=c(1,1,3,3)+0.1)
    plot(coordTab$x, coordTab$y,axes=F,ann=F, 
         xlim = c(-.2, .2), ylim = c(-.2, .2), 
         type = 'l');
    
    arrows(x0 = c(0,0, 0),y0 = c(0,0, 0),
           x1 = c(0.2, -0.2, -0.2),
           y1 = c(0, sqrt(0.12),-sqrt(0.12)), .05); ## draw custom axes
    arrows(x0 = c(0,0, 0),y0 = c(0,0, 0),
           x1 = c(-.2, 0.2, 0.2),
           y1 = c(0, -sqrt(0.12),sqrt(0.12)), .05, lty = 'dotted'); 
    text(x = coordTab$x-.002, y = coordTab$y+.002, 
         labels = coordTab$年份)
    
    return(propTab)
    
}


# 02. Workflow ------------------------------------------------------------
# File location identification
chinaMap <- './Data/China_ProvinceMap/province.shp'
triAexesData <- './Data/3-axisTest.xlsx'
marineEcoData <- './Data/vulnerabilityAssessmentData.xlsx'

# Import Data
provinceList <- c(210000, 120000, 130000, 370000,
                  320000, 310000, 330000, 350000,
                  440000, 450000, 150303) 
# Province filter: Liaoning, Hebei, Tianjin, Shandong, Jiangsu, Shanghai,
# Zhejiang, Fujian, Guangdong, Guangxi

# Read in files && pre-process
strIndustry <- readxl::read_excel(triAexesData)
marineIndEco <- readxl::read_excel(marineEcoData, sheet = 1)
eastChinaMap <- st_read(chinaMap) %>% 
    filter(GB %in% provinceList)

# Tri-Axes plot
propIndustry <- plotTriAxes(strIndustry)

# Weight calculation by *Entropy Method*
indi <- c(1, 1, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2) 
# Defination of positive and negative indicators

cla <- c('Sensitivity', 'Response', 'Sensitivity', 'Sensitivity', 'Sensitivity',
           'Sensitivity', 'Response', 'Response', 'Response', 'Response', 'Response',
           'Response')
X <- marineIndEco[,5:dim(marineIndEco)[2]] # Select real indicators

srWeight <- Entropy_Weight(X, index = indi)$w
sensWeight <- srWeight[which(cla == 'Sensitivity')]
resWeight <- srWeight[which(cla != 'Sensitivity')]
Sens <- X[, which(cla == 'Sensitivity')]
Res <- X[, which(cla != 'Sensitivity')]

impacts = c('+', '+', '-', '-', '+', '-', '-', '-', '-', '-', '-', '-')
senScore <- topsis(decision = as.matrix(Sens), weights = sensWeight, 
       impacts = c('+', '-', '-', '+', '-'))$score
reScore <- topsis(decision = as.matrix(Res), weights = resWeight, 
                  impacts = c('+', '-', '-', '-', '-', '-', '-'))$score

summaryWeighTable <- tibble::tibble(Sensitivity = as.matrix(Sens) %*% sensWeight,
                                    Response = as.matrix(Res) %*% resWeight) 
systemWeight <- Entropy_Weight(summaryWeighTable,
               index = c(1, 2))$w

vulnerabilityCH <- (summaryWeighTable$Sensitivity * systemWeight[1])/
    (summaryWeighTable$Response * systemWeight[2])
vulnerabilityCH <- topsis(decision = vulnerabilityCH, weights = 1, 
       impacts = '-')$score

vulTableCH <- marineIndEco[,1:4] %>%
    bind_cols(vulnerabilityCH) %>%
    rename(vulValue = ...5)

srFinal <- marineIndEco[,1:4] %>% 
    mutate(sensitivity = senScore,
           response = reScore) %>%
    group_by(Year) %>%
    summarise(sensitivity = sum(sensitivity),
              response = sum(response)) %>%
    ungroup()
vulFinal <- vulTableCH %>% 
    group_by(Year) %>% 
    summarise(vulValue = sum(vulValue)) %>%
    ungroup() 

# Unfinished second axis needed
vulFinal %>%
    full_join(srFinal, by = 'Year') %>%
    pivot_longer(-Year) %>%
    ggplot(aes(Year, value, group = name)) +
    geom_point() +
    geom_line(aes(linetype = name))

# Rank needed
scaledIndicator <- scale(X, center = F, scale = T)
scaledIndicator %*% diag(srWeight)
prop.table(scaledIndicator %*% diag(srWeight), margin = 1)

ggplot(eastChinaMap) +
    geom_sf() + theme_nothing() +
    geom_sf_text(aes(label = NAME_PINGY))



