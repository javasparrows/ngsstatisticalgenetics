# コード10
### fit null model
obj_nullmodel <- fit_nullmodel(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+as.factor(study_race),
	data=phenotype,
	kins=sgrm,
	use_sparse=TRUE,kins_cutoff=0.022,
	id="sample.id",
	family=gaussian(link="identity"),
	verbose=TRUE)

# コード11
LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+as.factor(company)

# コード12
family=binomial(link="logit")

# コード14
jobs_num = get(load("/Path/To/STAAR/Association/s0_prep/jobs_num.Rdata"))
cumsum(jobs_num$individual_analysis_num)

# コード23
jobs_num <- get(load("/Path/To/STAAR/Association/s0_prep/jobs_num.Rdata"))
cumsum(jobs_num$sliding_window_num)

# コード26
jobs_num <- get(load("/Path/To/STAAR/Association/s0_prep/jobs_num.Rdata"))
cumsum(jobs_num$scang_num)
