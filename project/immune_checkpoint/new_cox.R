rt_esc <- cbind(t(rt_T_m[grep(gnam, row.names(rt_T_m)), ]), rt_S_m[, c("X_OS_IND", "X_OS")], 
                rt_C_m[, c("age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", "pathologic_N", "pathologic_M", "tumor_stage.diagnoses")])
colnames(rt_esc) <- c(gnam, "status", "time",  "age", "gender", "T", "N", "M", "stage")
coxph(formula = Surv(time, status) ~ ENSG00000188389 + age + gender + T +  N + M + stage, data = rt_esc)
