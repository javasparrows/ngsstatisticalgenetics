# コード1
# ライブラリのインストール
install.packages("BeviMed") # バージョン5.4
install.packages("vcfR") # バージョン1.14.0
install.packages("tidyverse") # バージョン2.0.0
# ライブラリの読み込み
library(BeviMed)
library(vcfR)
library(tidyverse)

# コード2
genome_vcf <- read.vcfR("./FLG.vcf")
cadd <- read.table("./FLG_cadd.tsv", skip=1,
                   comment.char="", header=T) %>%
    rename(Chrom=X.Chrom) %>%
    mutate(ID=str_c(Chrom, Pos, Ref, Alt, sep="_")) %>%
    select(ID, RawScore, PHRED)
casectrl <- read.csv("./demodata_phenotype.csv")

# コード3
cadd_10list <- cadd %>%
    filter(PHRED <= 10) # CADDスコアが10以下のレアバリアントを選出．

# コード4
info <- vcfR2tidy(genome_vcf, info_only=TRUE)$fix %>% # VCFデータのINFO列を抽出．
    select(ID, CHROM, POS, ID, REF, ALT, ANN) # 指定した列を抽出．

# コード5
# HIGHバリアントグループを作成する場合
filtered_variant <- info %>%
    filter(str_detect(ANN, pattern="HIGH"))

# MODERATEバリアントグループを作成する場合
filtered_variant <- info %>%
    # IMPACT値がMODERATEかHIGHであるレアバリアントを抽出．
    filter(str_detect(ANN, pattern="MODERATE") |
           str_detect(ANN, pattern="HIGH")) %>%
    # IMPACT値がHIGHであるバリアントを1，
    # MODERATE であるバリアントを0とした列を追加．
    mutate(mode_high=if_else(str_detect(ANN, pattern="HIGH"), 1, 0))

# 5′UTRバリアントグループを作成する場合
filtered_variant <- info %>%
    filter(!(str_detect(ANN, pattern="MODERATE") |
             str_detect(ANN, pattern="HIGH")) & 
           str_detect(ANN, pattern="5_prime_UTR_variant"))

# 3′UTRバリアントグループを作成する場合
filtered_variant <- info %>%
    filter(!(str_detect(ANN, pattern="MODERATE") |
             str_detect(ANN, pattern="HIGH")) &
           str_detect(ANN, pattern="3_prime_UTR_variant"))

# コード6
genome_tidy <- vcfR2tidy(genome_vcf, format_fields="GT", single_frame=TRUE)
gt <- genome_tidy$dat %>%
    select(ID, Indiv, gt_GT) %>% # 解析に必要な列を抽出．
    # 遺伝型データの文字列を置換．
    mutate(allele_012=str_replace(gt_GT, pattern="0/0", replacement="0")) %>%
    mutate(allele_012=str_replace(allele_012, pattern="0/1", replacement="1")) %>%
    mutate(allele_012=str_replace(allele_012, pattern="1/1", replacement="2")) %>%
    mutate(allele_012=as.numeric(allele_012)) %>%
    rename(sampleID=Indiv)

# コード7
variant_table <- gt %>%
    filter(!(ID %in% cadd_10list$ID)) %>%
    filter(ID %in% filtered_variant$ID) %>%
    # 表形式に変換する．
    pivot_wider(id_cols="sampleID", names_from=ID, values_from=allele_012)

# コード8
# データの結合
bevimed_dataset <- casectrl %>%
    left_join(variant_table, by="sampleID") # 表現型データと遺伝型データを結合．
# データの作成
bevimed_y <- as.logical(bevimed_dataset$case_ctrl) # 表現型データのみを抽出．
bevimed_G <- bevimed_dataset %>%
    select(-sampleID, -case_ctrl) %>% # 遺伝型データのみを抽出．
    as.matrix() # データフレーム形式から行列形式に変換．
variant_prioritise <- filtered_variant$mode_high
# MODERATEバリアントグループの場合は，HIGH/MODERATEを示す列を抽出．

# コード9
# MODERATEバリアントグループの場合
bevimed(G=apply(bevimed_G, c(1,2), as.numeric), y=bevimed_y,
        omega_shape=c(2,8), variant_weights=variant_prioritise,
        standardise_weights=TRUE)
# HIGHバリアントグループの場合
bevimed(G=apply(bevimed_G, c(1,2), as.numeric), y=bevimed_y, omega_shape=c(2,1))
# 5′UTRと3′UTRバリアントグループの場合
bevimed(G=bevimed_G, y=bevimed_y, omega_shape=c(2,8))


# コード10
subject_count <- bevimed_dataset %>%
    # 集計に適した形式に変換．
    pivot_longer(cols=-c("sampleID", "case_ctrl"),
                 names_to="ID", values_to="gt") %>%
    # レアバリアントを保有していない人を除外．
    filter(gt != 0) %>%
    # 集計時に層別化する項目を設定．
    group_by(ID, case_ctrl, gt) %>%
    # 集計．
    summarise(count=n())
subject_count_withposition <- info %>%
    # アノテーションデータと集計したデータを結合．
    inner_join(subject_count, by="ID") %>% 
    # 集計結果にラベルを付与．
    mutate(case_ctrl_het_homo=case_when(
               case_ctrl == 0 & gt == 1 ~ "ctrl_hetero",
               case_ctrl == 0 & gt == 2 ~ "ctrl_homo",
               case_ctrl == 1 & gt == 1 ~ "case_hetero",
               case_ctrl == 1 & gt == 2 ~ "case_homo"))

# コード11
# 入力データと列を指定．
ggplot(data=subject_count_withposition,
       aes(x=reorder(x=ID, X=POS), y=count, fill=case_ctrl_het_homo)) +
    geom_bar(stat="identity") + # グラフの種類を指定．
    ylab("Subjects count") + # Y 軸のラベル名を設定．
    # 軸や軸ラベルを設定．
    theme(axis.title.x=element_blank(), axis.title=element_text(size=7), 
          axis.text=element_text(size=7), axis.text.x=element_text(angle=90),
          legend.title=element_blank(), legend.text=element_text(size=7)) +
    scale_fill_manual(values=c("firebrick1", "firebrick", "lightskyblue")) # 色

# コード12
result <- bevimed(G=bevimed_G, y=bevimed_y, omega_shape=c(2,8),
                  variant_weights=variant_prioritise, standardise_weights=TRUE)
extract_prob_pathogenic(result, by_model=FALSE)/ extract_prob_association(result)
# BeviMed の実行結果からレアバリアントごとの事後確率を取り出す．

# コード13
info_cadd <- info %>%
    left_join(cadd, by="ID") # snpEffとCADDスコアを結合．
info_cadd <- info_cadd %>%
    select(-RawScore) %>%
    replace_na(list(PHRED=30)) # CADDスコアのNAを30に置換．
filtered_variant <- info_cadd %>%
    filter(str_detect(ANN, pattern="MODERATE") | str_detect(ANN, pattern="HIGH"))
variant_prioritise_cadd <- filtered_variant$PHRED
bevimed(G=bevimed_G, y=bevimed_y, omega_shape=c(2,8),
        variant_weights=variant_prioritise)

# コード14
genename <- read.table("./genelist.txt", header=TRUE, sep="\t", fill=TRUE, quote="")
genelist <- genename$symbol # 遺伝子リストを作成．
# 結果をまとめる表の雛形を作成．
result_table <- data.frame(genename="test", group="test", prob="0")
for (i in 1:length(genelist)) { # 遺伝子リストで繰り返し処理を設定．
    gene_variant <- info_cadd %>%
        filter(str_detect(ANN, pattern=genelist[i]))
    variant_table <- gt %>%
        filter(ID %in% gene_variant$ID) %>%
        filter(!(ID %in% cadd_10list$ID)) %>%
        filter(ID %in% filtered_variant$ID) %>%
        pivot_wider(id_cols="sampleID", names_from=ID, values_from=allele_012)
    bevimed_dataset <- casectrl %>%
        inner_join(variant_table, by="sampleID")
    bevimed_y <- as.logical(bevimed_dataset$case_ctrl)
    bevimed_G <- as.matrix(bevimed_dataset[, -c(1,2)])
    colnames(bevimed_G) <- gsub("X", "", colnames(bevimed_G))
    variant_prioritise <- filtered_variant$mode_high
    result_temp <- data.frame(genename=genelist[i], group="moderate",
                              prior_prob_association=0.001,
                              prior_prob_dominant=0.5, omega_shape=c(2,8),
                              variant_weights=variant_prioritise,
                              standardise_weights=TRUE)
    # 格納する項目を選出．
    result_temp <- data.frame(genename=genelist[i], group="moderate",
                              prob=result)
    # 結果をresult_tableに格納．
    result_table <- bind_rows(result_table, result_temp)
}
result_table
