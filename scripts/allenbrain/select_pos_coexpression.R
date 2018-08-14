# This script loads co-expression computed in brain regions DG, CA1,CA2,CA3,CA4,subiculun and SptN using Allen brain atlas data

# Load datasets 
# Co-expression in DG
load(file="~/absb/results/allenbrain/DG_coexp_int.RData")
DG_coexp_int_ps <- DG_coexp_int[DG_coexp_int$score>0,]

# Co-expression in CA1
load(file="~/absb/results/allenbrain/CA1_coexp_int.RData")
CA1_coexp_int_ps <- CA1_coexp_int[CA1_coexp_int$score>0,]

# Co-expression in CA2
load(file="~/absb/results/allenbrain/CA2_coexp_int.RData")
CA2_coexp_int_ps <- CA2_coexp_int[CA2_coexp_int$score>0,]

# Co-expression in CA3
load(file="~/absb/results/allenbrain/CA3_coexp_int.RData")
CA3_coexp_int_ps <- CA3_coexp_int[CA3_coexp_int$score>0,]

# Co-expression in CA4
#load(file="~/absb/results/allenbrain/CA4_coexp_int.RData")
#CA4_coexp_int_ps <- CA4_coexp_int[CA4_coexp_int$score>0,]
load(file="~/absb/results/allenbrain/CA4_coexp_int_ps.RData")

# Co-expression in subiculum
load(file="~/absb/results/allenbrain/subiculum_coexp_int.RData")
subiculum_coexp_int_ps <- subiculum_coexp_int[subiculum_coexp_int$score>0,]

# Co-expression in SptN
load(file="~/absb/results/allenbrain/SptN_coexp_int.RData")
SptN_coexp_int_ps <- SptN_coexp_int[SptN_coexp_int$score>0,]


coexp_aba_ps <- rbind(DG_coexp_int_ps, CA1_coexp_int_ps, CA2_coexp_int_ps, CA3_coexp_int_ps, CA4_coexp_int_ps, subiculum_coexp_int_ps, SptN_coexp_int_ps) 
save(coexp_aba_ps, file="/absb/results/allenbrain/coexp_aba_ps.RData") 
