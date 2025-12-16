# コード1
# データの読み込み
data1 = read.table("./data.tsv", header=T)

# Wald法
res0_1 <- lm(X ~ G5, data = data1)
res0_2 <- lm(Y ~ G5, data = data1)
summary(res0_2)$coefficients[2] / summary(res0_1)$coefficients[2]

# PRSの計算
GX1 <- lm(X ~ G1, data=data1)
GX2 <- lm(X ~ G2, data=data1)
GX3 <- lm(X ~ G3, data=data1)
GX4 <- lm(X ~ G4, data=data1)
GX5 <- lm(X ~ G5, data=data1)

b_GX <- c(
    summary(GX1)$coefficients[2],
    summary(GX2)$coefficients[2],
    summary(GX3)$coefficients[2],
    summary(GX4)$coefficients[2],
    summary(GX5)$coefficients[2]
)

data1$PRS <- rep(0, nrow(data1))
for (i in 1:nrow(data1)) {
    data1$PRS[i] <- b_GX[1]*data1$G1[i] +
        b_GX[2]*data1$G2[i] + b_GX[3]*data1$G3[i] +
        b_GX[4]*data1$G4[i] + b_GX[5]*data1$G5[i]
}
head(data1)

# コード2
# 2SLS
res1 <- lm(X ~ PRS, data=data1)
data1$pX <- predict(res1, newdata=data1)
head(data1)

res2 <- lm(Y ~ pX, data=data1)
summary(res2)

# コード3
# MR-Egger
GY1 <- lm(Y ~ G1, data=data1)
GY2 <- lm(Y ~ G2, data=data1)
GY3 <- lm(Y ~ G3, data=data1)
GY4 <- lm(Y ~ G4, data=data1)
GY5 <- lm(Y ~ G5, data=data1)
b_GY <- c(
    summary(GY1)$coefficients[2],
    summary(GY2)$coefficients[2],
    summary(GY3)$coefficients[2],
    summary(GY4)$coefficients[2],
    summary(GY5)$coefficients[2]
)

data2 <- data.frame(t(rbind(b_GX, b_GY)))
head(data2)

res3 <- lm(b_GY ~ b_GX, data = data2)
summary(res3)

plot(data2$b_GX, data2$b_GY)

