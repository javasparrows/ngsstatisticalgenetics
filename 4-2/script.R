#コード6
install.packages("bigsnpr")

#コード7
install.packages(c("RcppArmadillo", "data.table", "Matrix"), dependencies=TRUE)
install.packages("devtools")
library("devtools")
devtools::install_github("tshmak/lassosum")

#コード8
install.packages("stmgp")

#コード12
dat <- read.table(gzfile("Height.QC.gz"), header=T)
dat$BETA <- log(dat$OR)
write.table(dat, "Height.QC.Transformed", quote=F, row.names=F)

#コード16
dat <- read.table(gzfile("Height.QC.gz"), header=T)
dat$BETA <- log(dat$OR)
dat <- dat[, c("SNP", "A1", "A2", "BETA", "P")]
write.table(dat, "Height.QC.Transformed2", sep="\t", quote=F, row.names=F)

#コード20
library(bigsnpr)
options(bigstatsr.check.parallel.blas=FALSE)
options(default.nproc.blas=NULL)

library(data.table)
library(magrittr)
phenotype <- fread("EUR.height")
covariate <- fread("EUR.cov")
pheno <- merge(phenotype, covariate)

# HapMap3 SNPsの読み込み
info <- readRDS(runonce::download_file("https://ndownloader.figshare.com/files/25503788", fname="map_hm3_ldpred2.rds"))

# GWAS要約統計量の読み込み
sumstats <- bigreadr::fread2("Height.QC.gz")
names(sumstats) <- c("chr", "pos", "rsid", "a1", "a0", "n_eff", "beta_se", "p", "OR", "INFO", "MAF")
sumstats$beta <- log(sumstats$OR)
# GWAS要約統計量のうちHapMap3に存在するSNPsに限定
sumstats <- sumstats[sumstats$rsid %in% info$rsid,]
### LD行列の計算（ここから）###
NCORES <- nb_cores()
tmp <- tempfile(tmpdir="tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
corr <- NULL
ld <- NULL
fam.order <- NULL
snp_readBed("EUR.QC.bed")
obj.bigSNP <- snp_attach("EUR.QC.rds")
map <- obj.bigSNP$map[-3]
names(map) <- c("chr", "rsid", "pos", "a1", "a0")

info_snp <- snp_match(sumstats, map)
genotype <- obj.bigSNP$genotypes
CHR <- map$chr
POS <- map$pos

# 1000ゲノムプロジェクトの遺伝子地図からSNPの位置情報を取得
POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")

for (chr in 1:22) {
    ind.chr <- which(info_snp$chr == chr)
    ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]

    corr0 <- snp_cor(
        genotype,
        ind.col = ind.chr2,
        ncores = NCORES,
        infos.pos = POS2[ind.chr2],
        size = 3 / 1000
    )
    if (chr == 1) {
        ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, tmp)
    } else {
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
    }
}
### LD行列の計算（ここまで）###
fam.order <- as.data.table(obj.bigSNP$fam)
setnames(fam.order, c("family.ID", "sample.ID"), c("FID", "IID"))

# LDスコア回帰の実施
df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ldsc <- snp_ldsc(ld, length(ld), chi2 = (df_beta$beta / df_beta$beta_se)^2,
sample_size = df_beta$n_eff, blocks = NULL)
h2_est <- ldsc[["h2"]]

# LDpred2-infによる効果サイズβの推定
beta_inf <- snp_ldpred2_inf(corr, df_beta, h2=h2_est)

if(is.null(obj.bigSNP)){
    obj.bigSNP <- snp_attach("EUR.QC.rds")
}
genotype <- obj.bigSNP$genotypes

# 全サンプルに対してPRSを計算
ind.test <- 1:nrow(genotype)
pred_inf <- big_prodVec(genotype, beta_inf, ind.row = ind.test, ind.col = info_snp$`_NUM_ID_`)
pred_inf <- cbind(fam.order, pred_inf)
write.table(pred_inf, paste("inf-answer.txt", sep=""), quote=F, row.names=F, append=F)

#コード23
library(lassosum)
library(data.table)
library(methods)
library(magrittr)
library(parallel)
cl <- makeCluster(2)

sum.stat <- "Height.QC.gz"
bfile <- "EUR.QC"
# Read in and process the covariates
covariate <- fread("EUR.cov")
pcs <- fread("EUR.eigenvec") %>%
    setnames(., colnames(.), c("FID", "IID", paste0("PC", 1:6)))
# Need as.data.frame here as lassosum doesn't handle data.table
# covariates very well
cov <- merge(covariate, pcs)

ld.file <- "EUR.hg19"
prefix <- "EUR"
target.pheno <- fread("EUR.height")[, c("FID", "IID", "Height")]
ss <- fread(sum.stat)
# p値がゼロのSNPは後続の解析でエラーが発生するため除外
ss <- ss[!P == 0]

# p値からSNP間相関を計算
cor <- p2cor(p=ss$P, n=ss$N, sign=log(ss$OR))
fam <- fread(paste0(bfile, ".fam"))
fam[, ID:=do.call(paste, c(.SD, sep=":")), .SDcols=c(1:2)]

# lassosumパイプラインの実行
out <- lassosum.pipeline(
    cor=cor,
    chr=ss$CHR,
    pos=ss$BP,
    A1=ss$A1,
    A2=ss$A2,
    ref.bfile=bfile,
    test.bfile=bfile,
    LDblocks=ld.file,
    cluster=cl
)
target.res <- validate(out, pheno=as.data.frame(target.pheno), covar=as.data.frame(cov))
write.table(target.res$results.table, paste("lassosum.txt", sep=""), quote=F, row.names=F, append=F)

#コード26
library(stmgp)
library(tidyverse)

n.snp = length(readLines("EUR.QC.bim"))

pheno <- read.table("EUR.height", header=T)
fam <- read.table("EUR.QC.fam", header=F)
colnames(fam) <- c("FID", "IID", "V3", "V4", "V5", "V6")

fam <- left_join(fam, pheno, by=c("FID", "IID"))
fam <- fam[, c("FID", "IID", "V3", "V4", "V5", "Height")]
fam[is.na(fam)] <- -9
write.table(fam, paste("EUR.QC.b.fam", sep=""), quote=F, row.names=F, append=F, col.names=F)

stmgp = stmgplink(trainbed=paste0("EUR.QC", c(".bed", ".bim", ".b.fam")), testbed=paste0("EUR.QC", c(".bed", ".bim", ".b.fam")), maxal=5000/n.snp)

head(stmgp$PE)

write.table(stmgp$PEte, paste("stmgp.txt", sep=""), quote=F, row.names=F, append=F)



