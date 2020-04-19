### Set OBJFUNC and validation function parameters
DF0 <- 4 #degrees of freedom for Scaled-Inverse Chi-square destribution of initial random effects
DF1 <- 1 #degrees of freedom for Scaled-Inverse Chi-square destribution of residuals
SCALE <- 1 #scale for Scaled-Inverse Chi-square destribution random effects and residuals
OBJFUNC.ARGS <- list(df0=DF0, df1=DF1, S0=SCALE, nIter=500, burnIn=300, thin=2, bs=50, saveAt='temp/MTM_')
VAL.ARGS <- list(df0=DF0, df1=DF1, S0=SCALE, nIter=1500, burnIn=1000, thin=5, bs=50, saveAt='temp/MTM_')
#
OBJFUNC_Parameters <- c(OBJFUNC.ARGS$df0, OBJFUNC.ARGS$df1, OBJFUNC.ARGS$S0, OBJFUNC.ARGS$nIter, OBJFUNC.ARGS$burnIn, OBJFUNC.ARGS$thin, OBJFUNC.ARGS$bs)
names(OBJFUNC_Parameters) <- c('OBJ_df0', 'OBJ_df1', 'OBJ_S0', 'OBJ_nIter', 'OBJ_burnIn', 'OBJ_thin', 'OBJ_bs')
VALID_Parameters <- c(VAL.ARGS$df0, VAL.ARGS$df1, VAL.ARGS$S0, VAL.ARGS$nIter, VAL.ARGS$burnIn, VAL.ARGS$thin, VAL.ARGS$bs)
names(VALID_Parameters) <- c('VAL_df0', 'VAL_df1', 'VAL_S0', 'VAL_nIter', 'VAL_burnIn', 'VAL_thin', 'VAL_bs')
