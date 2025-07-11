
#Baran KILINC, Medine YAZICI, Fitnat KOC

library(MVN)
library(dplyr)
library(forecast)
library(purrr)
library(corrplot)
library(ICSNP)
library(car)
library(GGally)
library (psych)
library(heplots)
library(gridExtra)
library(factoextra)
library(psych)
library(MASS)
library(klaR)
library(ggplot2)
library(GGally)
library(mlbench)
library(rstatix)
library(CCA)


##IMPORTING AND CLEANING THE DATA

df <- read.csv("NFH.csv")
dim(df)
sum(is.na(df)) #no NA values
sum(duplicated(df))#no duplicated obs.


df <- df[!grepl("\\*", df$Pregnant.women.age.15.49.years.who.are.anaemic...11.0.g.dl.22....), ] # "*" is NA values, so we extracted these from the data.
df$Pregnant.women.age.15.49.years.who.are.anaemic...11.0.g.dl.22.... <- gsub('[\\"()]', '', df$Pregnant.women.age.15.49.years.who.are.anaemic...11.0.g.dl.22....)
df$Pregnant.women.age.15.49.years.who.are.anaemic...11.0.g.dl.22.... <- as.numeric(df$Pregnant.women.age.15.49.years.who.are.anaemic...11.0.g.dl.22....)
class(df$Pregnant.women.age.15.49.years.who.are.anaemic...11.0.g.dl.22....)

df <- df %>%
  mutate(Category = case_when(
    State.UT %in% c("Kerala", "Tamil Nadu", "Karnataka", "Punjab", "Maharashtra", "NCT of Delhi",
                 "Gujarat", "Haryana", "Rajasthan", "Telangana", "Uttar Pradesh", "West Bengal",
                 "Maharastra") ~ "Urban", 
    State.UT %in% c("Andhra Pradesh", "Arunachal Pradesh", "Assam", "Bihar", "Chhattisgarh",
                 "Dadra and Nagar Haveli & Daman and Diu", "Himachal Pradesh", "Jammu & Kashmir",
                 "Jharkhand", "Ladakh", "Lakshadweep", "Madhya Pradesh", "Manipur", "Meghalaya",
                 "Mizoram", "Nagaland", "Odisha", "Puducherry", "Sikkim", "Tripura", "Uttarakhand",
                 " Lakshadweep ") ~ "Rural"
  )) # After examining the cities in terms of health facilities and development, a new column was added as rural and urban.

df <- df %>%
  mutate(Category_Numeric = ifelse(Category == "Urban", 1, 0))

df <- df %>%
  mutate(ObeseSituation = ifelse(df$Women..age.15.49.years..whose.Body.Mass.Index..BMI..is.below.normal..BMI..18.5.kg.m2.21.... > 25, 1, 0))

df$Category <- factor(df$Category)
df$ObeseSituation <- factor(df$ObeseSituation)
df$Category_Numeric <- factor(df$Category_Numeric)


dim(df)
sum(is.na(df)) #no NA values
sum(duplicated(df))#no duplicated obs.


##ASSUMPTIONS

#Normality Checking

numeric_df <- select_if(df, is.numeric)

# We are interested in anemia in women, so we proceed with the necessary variables.
numeric_df <- numeric_df[, c(35,36,37,38,39,40,41,42,47,48,49,56)]
dim(numeric_df)
colnames(numeric_df)

#We shorten variable names to see plots more easily.
names(numeric_df)[names(numeric_df) == "Women..age.15.49.years..whose.Body.Mass.Index..BMI..is.below.normal..BMI..18.5.kg.m2.21...."] <- "BMI"
names(numeric_df)[names(numeric_df) == "Women..age.15.49.years..who.are.overweight.or.obese..BMI..25.0.kg.m2.21...."] <- "Obese"
names(numeric_df)[names(numeric_df) == "Women..age.15.49.years..who.have.high.risk.waist.to.hip.ratio...0.85....."] <- "WaisToHipRatio"
names(numeric_df)[names(numeric_df) == "Non.pregnant.women.age.15.49.years.who.are.anaemic...12.0.g.dl.22...."] <- "NonPregnantAnaemic"
names(numeric_df)[names(numeric_df) == "Pregnant.women.age.15.49.years.who.are.anaemic...11.0.g.dl.22...."] <- "PregnantAnameic"
names(numeric_df)[names(numeric_df) == "All.women.age.15.49.years.who.are.anaemic22...."] <- "AllAnaemic"
names(numeric_df)[names(numeric_df) == "Women..age.15.years.and.above.with.high..141.160.mg.dl..Blood.sugar.level23...."] <- "HighBloodSugarLevel"
names(numeric_df)[names(numeric_df) == "Women.age.15.years.and.above.wih.very.high...160.mg.dl..Blood.sugar.level23...."] <- "VeryHighBloodSugarLevel"
names(numeric_df)[names(numeric_df) == "Women.age.15.years.and.above.wih.Mildly.elevated.blood.pressure..Systolic.140.159.mm.of.Hg.and.or.Diastolic.90.99.mm.of.Hg....."] <- "MildlyBloodPressure"
names(numeric_df)[names(numeric_df) == "Women.age.15.years.and.above.wih.Moderately.or.severely.elevated.blood.pressure..Systolic..160.mm.of.Hg.and.or.Diastolic..100.mm.of.Hg....."] <- "ModeratelyBloodPressure."
names(numeric_df)[names(numeric_df) == "Women.age.15.years.and.above.wih.Elevated.blood.pressure..Systolic..140.mm.of.Hg.and.or.Diastolic..90.mm.of.Hg..or.taking.medicine.to.control.blood.pressure...."] <- "ElevatedBloodPressure"
names(numeric_df)[names(numeric_df) == "Women.age.15.years.and.above.who.use.any.kind.of.tobacco...."] <- "Tobacco"


mvn(numeric_df, multivariatePlot= "qq")

par(mar=c(3, 3, 2, 1))

result <- mvn(data = numeric_df, mvnTest = "royston")
result <- mvn(data = numeric_df, mvnTest = "royston", univariatePlot = "qqplot")
result <- mvn(data = numeric_df, mvnTest = "royston", univariatePlot = "histogram")


result$multivariateNormality#MULTIVARIATE NORMALITY IS NOT SATISFIED.
result$univariateNormality#NOT ALL VARIABLES ARE DIST. NORMALLY(only Women..age.15.49.years..who.have.high.risk.waist.to.hip.ratio...0.85..... is normally) 


#This Yeo-Johnson transformation is meaningful for these variables 
#because the log(y+1) transformation, which is similar to the sqrt 
#transformation, is applied to the resulting lambda value. Log(y+1) if lambda=0 and y>0.
#Therefore, our lambda's values for these variables are equal to 0 and values>0.
#There is no problem with shifting meaning between variables.
#The sqrt() transformation, which provides normality and brings our other variables closer to normality, is applied.
Yeo_trans_non_preg_anaemic <- powerTransform(numeric_df$NonPregnantAnaemic, family = "yjPower")
summary(Yeo_trans_non_preg_anaemic)
transformed_data_all_anaemic <- yjPower(numeric_df$NonPregnantAnaemic, lambda = Yeo_trans_non_preg_anaemic$lambda)
numeric_df$NonPregnantAnaemic<- transformed_data_all_anaemic

Yeo_trans_preg_anaemic <- powerTransform(numeric_df$PregnantAnameic, family = "yjPower")
summary(Yeo_trans_preg_anaemic)
transformed_data_anaemic <- yjPower(numeric_df$PregnantAnameic, lambda = Yeo_trans_preg_anaemic$lambda)
numeric_df$PregnantAnameic<- transformed_data_anaemic

Yeo_trans_all_anaemic <- powerTransform(numeric_df$AllAnaemic, family = "yjPower")
summary(Yeo_trans_all_anaemic)
transformed_data_all_anaemic <- yjPower(numeric_df$AllAnaemic, lambda = Yeo_trans_all_anaemic$lambda)
numeric_df$AllAnaemic<- transformed_data_all_anaemic


#sqrt transformation
sqrt_transformed<- sqrt(numeric_df[,c(2,7,8,9,10,11,12)])

numeric_df_transformed <- cbind(numeric_df[, -c(2,7,8,9,10,11,12)], sqrt_transformed)

result <- mvn(data = numeric_df_transformed, mvnTest = "royston")
result <- mvn(data = numeric_df_transformed, mvnTest = "royston", univariatePlot = "qqplot")
result <- mvn(data = numeric_df_transformed, mvnTest = "royston", univariatePlot = "histogram")

result$multivariateNormality#MULTIVARIATE NORMALTY IS SATISFIED
result$univariateNormality#ALL VARIABLES (6 variables) ARE DISTRIBUTED NORMALLY(except 5 variables, we assume normally.
# Because its p-values are close to 0.05, and both histogram and QQ-plot show that
#this variable resembles normal distribution. Therefore, we may assume this
#variable as normally distributed in the next steps of the analysis.

colnames(numeric_df_transformed)

df_transformed <- cbind(df$Category,df$Category_Numeric,df$ObeseSituation, numeric_df_transformed)
names(df_transformed)[names(df_transformed) == "df$Category"] <- "Category"
names(df_transformed)[names(df_transformed) == "df$ObeseSituation"] <- "ObeseSituation"
names(df_transformed)[names(df_transformed) == "df$Category_Numeric"] <- "CategoricalBinary"

colnames(df_transformed)

#1-Identifying univariate outliers
for (i in 13:length(df_transformed)) {
  out <- df_transformed %>%
    identify_outliers(colnames(df_transformed)[i])%>%
    as.data.frame()
  if(count(out)!=0){
    print(colnames(df_transformed)[i])
    print(out)
  }
}

#We have some outliers but non of them are extreme value. For now,
#we can keep them in data.

out <- numeric_df_transformed %>% 
  mahalanobis_distance() %>%
  filter(is.outlier == T)%>%
  as.data.frame()
out
#not same observation

#2-Univariate normality
#We check this assumption above, and it satisfied.

#3-Multivariate normality
#We check this assumption above, and it satisfied.

#4-Check multicollinearity between the variables
testRes = cor.mtest(numeric_df_transformed, conf.level = 0.95)
a <- cor(numeric_df_transformed)
plot.new()
dev.off() 
corrplot(a, p.mat = testRes$p, method = 'circle', type = 'lower', insig='blank',
         addCoef.col ='black', number.cex = 0.8, order = 'AOE', diag=FALSE,tl.cex = 0.7)

#5. Check linearity
ggpairs(numeric_df_transformed)
#According to the plots of pairs, we assume linearity assumption in satisfied.


##Hypothesis Testing in Multivariate Analysis
#INFERENCE ABOUT MEAN
y<-df_transformed%>%select(PregnantAnameic, NonPregnantAnaemic)
y

colnames(df_transformed)
xbar = colMeans(y)
xbar

test<-mvn(y,mvnTest = "mardia")
test$multivariateNormality
test$univariateNormality

error.bars (y, ylab="Group Means", xlab=" Dependent Variables")

pregnant_anaemic <- lm(PregnantAnameic ~ 1, data = df_transformed)
confint(pregnant_anaemic)

non_preg_anaemic <- lm(NonPregnantAnaemic ~ 1, data = df_transformed)
confint(non_preg_anaemic)

mu0=c(50,50)

HotellingsT2(y,mu=mu0)

#Since we reject H0. Therefore, we don't have enough evidence to 
#conclude that the the mean vector equals to c(50,50).

#Two Independent Samples
#H0: variance-covariance matrices are equal for each combination formed by each group in the independent variable

subset_data <- df_transformed%>%select(PregnantAnameic,NonPregnantAnaemic,Category)
table(subset_data$Category)

boxM(Y = cbind(subset_data$PregnantAnameic,subset_data$NonPregnantAnaemic), group = factor(subset_data$Category))
#since p-value>0, variance-covariance matrices are equal for each dependent

HotellingsT2(cbind(subset_data$PregnantAnameic,subset_data$NonPregnantAnaemic) ~ subset_data$Category)

#We fail to reject H0. Therefore, we don't have enough evidence to prove that the mean of responses change with respect to category.
#however, p-value(0.062) close to 0.05. 


#ONE-WAY MANOVA
subset_data1 <- df_transformed %>% 
  select(PregnantAnameic, NonPregnantAnaemic, Category) %>% 
  mutate(
    PregnantAnaemic = PregnantAnameic,
    NonPregnant_Anaemic = NonPregnantAnaemic
  )

subset_data1 %>% head()

subset_data1 %>% group_by(Category) %>%  summarise(n = n(), 
                                               mean_pregnant = mean(PregnantAnameic), 
                                               sd_pregnant = sd(PregnantAnameic),
                                               mean_nonpregnant = mean(NonPregnantAnaemic),
                                               sd_nonpregnant = sd(NonPregnantAnaemic))

p1 <- ggplot(subset_data1, aes(x = Category, y = PregnantAnameic, fill = Category)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) + theme(legend.position="top")+theme_minimal()+
  labs(title = "The Box Plot of Pregnant Anaemic by Category.")
p2 <- ggplot(subset_data1, aes(x = Category, y = NonPregnantAnaemic, fill = Category)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) + theme(legend.position="top")+theme_minimal()+
  labs(title = "The Box Plot of Non Pregnant by Category.")
grid.arrange(p1, p2, ncol=2) 


#Null hypothesis: variance-covariance matrices are equal for each combination formed by each group in the independent variable

boxM(Y = cbind(subset_data1$PregnantAnaemic,subset_data1$NonPregnantAnaemic), group = factor(subset_data1$Category))

m1 <- manova(cbind(PregnantAnameic,NonPregnantAnaemic) ~ Category, data = subset_data1)
summary(m1)

#"Therefore, we are not 95% confident that at least one category type is significantly 
#different from the others, as the results do not provide sufficient evidence to 
#support a significant difference among the ruraal and urban.

subset_data1 <- df_transformed %>% 
select(PregnantAnameic, NonPregnantAnaemic, Category,ObeseSituation) %>% 
  mutate(
    PregnantAnaemic = PregnantAnameic,
    NonPregnant_Anaemic = NonPregnantAnaemic
  )


boxM(Y = cbind(df_transformed$PregnantAnameic, df_transformed$NonPregnantAnaemic), group = factor(subset_data1$ObeseSituation))

m2 <- manova(cbind(PregnantAnameic,NonPregnantAnaemic) ~ Category*ObeseSituation, data = subset_data1)
summary(m2)


#regression for PCA

mlm1 <- lm(cbind(PregnantAnameic, NonPregnantAnaemic) ~  Category + BMI + WaisToHipRatio + 
             Obese + HighBloodSugarLevel + VeryHighBloodSugarLevel + MildlyBloodPressure + 
             ModeratelyBloodPressure. + ElevatedBloodPressure, data = df_transformed)
summary(mlm1)

m1 <- lm(cbind(PregnantAnameic) ~Category + BMI + WaisToHipRatio + 
           Obese + HighBloodSugarLevel + VeryHighBloodSugarLevel + MildlyBloodPressure + 
           ModeratelyBloodPressure. + ElevatedBloodPressure, data = df_transformed)
summary(m1)

m2 <- lm(cbind(NonPregnantAnaemic) ~ Category + BMI + WaisToHipRatio + 
           Obese + HighBloodSugarLevel + VeryHighBloodSugarLevel + MildlyBloodPressure + 
           ModeratelyBloodPressure. + ElevatedBloodPressure, data = df_transformed)
summary(m2)

head(resid(mlm1))

head(fitted(mlm1))

coef(mlm1)

sigma(mlm1)

vif(m1)

# Remove PregnantAnaemic and NonPregnantAnaemic from df_transformed
predictor_data<- df_transformed[, !(colnames(df_transformed) %in% c("AllAnaemic", "NonPregnantAnaemic"))]
# Verify the removal
colnames(predictor_data)
step(lm(cbind(PregnantAnameic) ~ ., data = predictor_data, direction = "backward"))
##to fix adjr^2 low values we used step function to choose proper model for our dataset.
#After output we choosed lowest AIC value of model

#PCA
colnames(numeric_df_transformed)
scaled_df <- as.data.frame(scale(numeric_df_transformed))

cov(scaled_df)
#For PCA, the covariance matrix is critical.PCA identifies principal components by calculating 
#eigenvalues and eigenvectors of the covariance matrix.pca1 <- prcomp(scaled_df)


scatterplotMatrix(scaled_df,diagonal = "histogram")

numeric_scaled_data <- scaled_df[sapply(scaled_df, is.numeric)]
numeric_scaled_data$AllAnaemic <- NULL
res <- cor(numeric_scaled_data, method = "pearson")
corrplot::corrplot(res, method= "color", order = "hclust")

pca1 <- prcomp(numeric_scaled_data)
summary(pca1)

names(pca1)
pca1$rotation
pca1$x
pca1$sdev
fviz_eig(pca1,addlabels=TRUE) #represent the proportion values

pca<-pca1$x[,1:4]
head(pca)

res1 <- cor(pca, method="pearson")
corrplot::corrplot(res1, method= "color", order = "hclust")

biplot(pca1, col = c("gray", "black"))

fviz_pca_var(pca1, col.var = "contrib")

fviz_pca_var(pca1, select.var = list(contrib = 4))

fviz_pca_ind(pca1, col.ind = "#00AFBB")

fviz_contrib(pca1, choice = "ind", axes = 1:2) + coord_flip()

fviz_pca_ind(pca1, label="none", habillage=df_transformed$Category,
             addEllipses=TRUE, ellipse.level=0.95)


#Factor Analysis

cm <- cor(numeric_scaled_data, method="pearson")
cm

corrplot::corrplot(cm, method= "number", order = "hclust")

KMO(r=cm)

print(cortest.bartlett(cm,nrow(numeric_scaled_data)))

parallel <- fa.parallel(numeric_scaled_data, fm = "minres", fa = "fa")

factanal(numeric_scaled_data, factors = 4)$PVAL

factanal(numeric_scaled_data, factors = 6)$PVAL

f<-factanal(numeric_scaled_data, factors = 6)
f

load <- f$loadings[, 1:2]
plot(load, type = "n")  # Create an empty plot
text(load, labels = names(numeric_scaled_data), cex = 0.7)

names(f$loadings[,1])[abs(f$loadings[,1])>0.4]

f1<-numeric_scaled_data[,names(f$loadings[,1])[abs(f$loadings[,1])>0.4]]
summary(alpha(f1, check.keys=TRUE))

scores<-factanal(numeric_scaled_data, factors = 6,scores="regression")$scores
head(scores)

cm1 <- cor(scores, method="pearson")
corrplot::corrplot(cm1, method= "number", order = "hclust")

#factor analysis

ggpairs(numeric_df_transformed)
#According to the plots of pairs, we assume linearity assumption in satisfied.

df_transformed_factor <- df_transformed
df_transformed_factor$Category <- NULL
df_transformed_factor$CategoricalBinary <- NULL
df_transformed_factor$Obese <- NULL

GGally::ggpairs(df_transformed_factor,  aes(color = ObeseSituation,  # Color by group (cat. variable)
                                          alpha = 0.5))
colnames(df_transformed_factor)
model <- lda(ObeseSituation~.,data = df_transformed_factor)
model

par(mar = c(4, 4, 2, 1))
plot(model)

model.values <- predict(model)
names(model.values)
par(mar = c(4, 4, 2, 1))
partimat(as.factor(ObeseSituation)~.,data=df_transformed_factor,method="lda") 
colnames(df_transformed_factor)


# Cluster Analysis
# We are interested in anemia in women, so we proceed with the necessary variables.
head(df_transformed_factor)

# Sadece sayD1sal deDiEkenleri seC'
numeric_data_clus <- select_if(df_transformed_factor, is.numeric)
# Veriyi standardize et (normalizasyon)
scaled_data <- scale(numeric_data_clus)
#Visualization Level Plot
# Compute the distance matrix
data_dist <- dist(scaled_data)
par(mfrow = c(2, 2), mar = c(1, 2, 1, 2))  # 2x2 grid ve margin ayarD1

# Single Linkage
plot(hclust(data_dist, method = "single"), main = "Single Linkage")

# Complete Linkage
plot(hclust(data_dist, method = "complete"), main = "Complete Linkage")

# Average Linkage
plot(hclust(data_dist, method = "average"), main = "Average Linkage")

# Ward's Method
plot(hclust(data_dist, method = "ward.D2"), main = "Ward Method")

par(mfrow=c(1,1))
cc <- hclust(data_dist)
plot(cc, main = "Clustering Dendrogram", xlab = "Observations", ylab = "Height")
library(lattice)

# Convert the distance matrix to a matrix
dist_matrix <- as.matrix(data_dist)

# Plot the distance matrix
levelplot(dist_matrix, xlab = "Observations", ylab = "Observations", 
          main = "Distance Matrix Heatmap")


# Calculate the range for each variable in the original data
rge <- apply(numeric_data_clus, 2, function(x) max(x) - min(x))
print(rge)
scaled_data_s<-sweep(numeric_data_clus,2,rge,FUN="/")
sapply(scaled_data_s, var)

#For determine cluster number
# Number of observations in your dataset
n <- nrow(scaled_data)
# Initialize the WSS (Within-Cluster Sum of Squares) vector
wss <- rep(0, 6)

# Calculate WSS for 1 cluster
wss[1] <- (n - 1) * sum(sapply(scaled_data, var))

# Calculate WSS for 2 to 6 clusters
for (i in 2:6) {
  wss[i] <- sum(kmeans(scaled_data, centers = i)$withinss)
}

# Plot WSS for each number of clusters
plot(1:6, wss, type = "b", xlab = "Number of clusters", ylab = "Within groups sum of squares", 
     main = "Elbow Method for Optimal Clusters")
#K-MEANS with 3 cluster
kmeans(scaled_data_s, centers = 3)$centers*rge
#cluster number for each row
kmeans(scaled_data_s, centers = 3)$cluster

# Multidimensional Scaling (MDS)
# Perform classical multidimensional scaling
mds_result <- cmdscale(data_dist, k = 13, eig = TRUE)

# View eigenvalues
print(mds_result$eig)
cumsum(abs(mds_result$eig)) / sum(abs(mds_result$eig))
cumsum((mds_result$eig)^2) / sum((mds_result$eig)^2)

#visiualization of boyut
x <- mds_result$points[, 1]
y <- mds_result$points[, 2]

plot(x, y, xlab = "Coordinate 1", ylab = "Coordinate 2", 
     main = "MDS Plot", type = "n")
text(x, y, labels = rownames(scaled_data), cex = 0.7)

group_labels <- factor(kmeans_result$cluster)  # Crnek grup etiketleri
plot(x, y, col = group_labels, pch = 19, xlab = "Coordinate 1", ylab = "Coordinate 2")
legend("topright", legend = levels(group_labels), col = 1:length(levels(group_labels)), pch = 19)



#cannonical, x=nonpregnant, y=obese
X <- data.frame(numeric_df_transformed$HighBloodSugarLevel,
                           numeric_df_transformed$VeryHighBloodSugarLevel, numeric_df_transformed$MildlyBloodPressure)
                          
Y <- data.frame(numeric_df_transformed$WaisToHipRatio, numeric_df_transformed$BMI, numeric_df_transformed$Tobacco)
cca_result <- cancor(X, Y)

print(cca_result)

# Canonical Correlations
cca_result$cor

cca_result$xcoef  # NonPregnant coefficients
cca_result$ycoef  # Obese coefficients

X_scores <- as.matrix(X) %*% cca_result$xcoef[, 1]
Y_scores <- as.matrix(Y) %*% cca_result$ycoef[, 1]

plot(X_scores, Y_scores, main = "Canonical Correlation",
     xlab = "Canonical Variable 1 (X)", ylab = "Canonical Variable 1 (Y)")
abline(lm(Y_scores ~ X_scores), col = "red")

# Wilks' Lambda testi
cca_test <- matcor(X, Y)  # X and Y  correlation matrisleri
summary(cca_result)       # Wilks' Lambda
