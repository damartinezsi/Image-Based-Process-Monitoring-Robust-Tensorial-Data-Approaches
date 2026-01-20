#==============================================================================#
# ---------------------------------------------------------------------------- #
#              Main - Simulaciones Trabajo de grado D.A.M.S.                   #
# ---------------------------------------------------------------------------- #
#==============================================================================#

# Generar los cluster para procesamiento paralelo

setwd("~/Tesis/Versión final - code")
library(parallel)

# Definir el número de clusters
num_cores <- detectCores() - 1 

# Crear el cluster
cl <- makeCluster(num_cores)

# Exportar los objetos necesarios al cluster
clusterExport(cl, varlist = c("run_simulation","generate_ellipse_image",
                              "generate_samples", "recon_trod_tensor",
                              "vecnorm", "makeDirections","unimcd",
                              "MacroPARAFAC2","Best_alpha","proj_qr"
                              ), envir = environment())

clusterEvalQ(cl, {
  library(grid)
  library(abind)
  library(rTensor)
  library(pracma)     
  library(stats)      
  library(cellWise)
  library(robustbase)
  library(rrcov)
  library(RobStatTM)
  library(StableMCD)
})

# ============================================================================ #
#                Hallar valores de alpha y varianza explicada                  #
# ============================================================================ #

clusterSetRNGStream(cl, 520)
N_iter <-1000
results_ba <- parLapply(cl, 1:N_iter, Best_alpha)
results_df_besta <- do.call(rbind, results_ba)

apply(results_df_besta,2,mean)

# ============================================================================ #
#        Medir el tiempo de descomposición - Hallar los UCL                    #
# ============================================================================ #

clusterSetRNGStream(cl, 520)
results <- parLapply(cl, 1:N_iter, run_simulation)
results_df <- do.call(rbind, results)
stopCluster(cl)

# Valores de alpha para n=80 Imágenes ---------------------------------------

alpha_TROD <-0.056        
alpha_TROD_MCD <-0.054  

alpha_MPCA <-0.054      
alpha_MPCA_MCD <-0.052    
alpha_MPCA_MCD_80 <-0.050   

alpha_MACRO <-0.05     
alpha_MACRO_MCD <-0.05    
alpha_MACRO_MCD_80 <-0.050


# UCLs basados en percentiles empíricos --------------------------------------
#Cartas T2 --------------------------------------------------------------------

UCL_T_TROD <- quantile(results_df$max_T2_TROD, probs = 1 - (alpha_TROD / 2))
UCL_T_TROD_MCD <- quantile(results_df$max_T2_TROD_MCD, probs = 1 - (alpha_TROD_MCD / 2))

UCL_T_MPCA <- quantile(results_df$max_T2_MPCA, probs = 1 - (alpha_MPCA / 2))
UCL_T_MPCA_MCD <- quantile(results_df$max_T2_MPCA_MCD, probs = 1 - (alpha_MPCA_MCD / 2))
UCL_T_MPCA_MCD_80 <- quantile(results_df$max_T2_MPCA_MCD_80, probs = 1 - (alpha_MPCA_MCD_80 / 2))

UCL_T_MACRO <- quantile(results_df$max_T2_MACRO, probs = 1 - (alpha_MACRO / 2))
UCL_T_MACRO_MCD <- quantile(results_df$max_T2_MACRO_MCD, probs = 1 - (alpha_MACRO_MCD / 2))
UCL_T_MACRO_MCD_80 <- quantile(results_df$max_T2_MACRO_MCD_80, probs = 1 - (alpha_MACRO_MCD_80 / 2))

#Cartas Q ----------------------------------------------------------------------
UCL_Q_TROD <- quantile(results_df$max_Q_TROD, probs = 1 - (alpha_TROD / 2))
UCL_Q_TROD_MCD <- quantile(results_df$max_Q_TROD, probs = 1 - (alpha_TROD_MCD / 2))

UCL_Q_MPCA <- quantile(results_df$max_Q_MPCA, probs = 1 - (alpha_MPCA / 2))
UCL_Q_MPCA_MCD <- quantile(results_df$max_Q_MPCA, probs = 1 - (alpha_MPCA_MCD / 2))
UCL_Q_MPCA_MCD_80 <- quantile(results_df$max_Q_MPCA, probs = 1 - (alpha_MPCA_MCD_80 / 2))

UCL_Q_MACRO <- quantile(results_df$max_Q_MACRO, probs = 1 - (alpha_MACRO / 2))
UCL_Q_MACRO_MCD <- quantile(results_df$max_Q_MACRO, probs = 1 - (alpha_MACRO_MCD / 2))
UCL_Q_MACRO_MCD_80 <- quantile(results_df$max_Q_MACRO, probs = 1 - (alpha_MACRO_MCD_80 / 2))


# ============================================================================ #
#           CARTAS CONJUNTAS T²–Q PARA TODOS LOS MÉTODOS Y VARIANTES           #
# ============================================================================ #

# ------------------ TROD ------------------
results_df$SIGNAL_TROD <- as.integer(
  (results_df$max_T2_TROD > UCL_T_TROD) |
    (results_df$max_Q_TROD  > UCL_Q_TROD)
)
results_df$SIGNAL_TROD_MCD <- as.integer(
  (results_df$max_T2_TROD_MCD > UCL_T_TROD_MCD) |
    (results_df$max_Q_TROD      > UCL_Q_TROD_MCD)
)

# ------------------ MPCA ------------------
results_df$SIGNAL_MPCA <- as.integer(
  (results_df$max_T2_MPCA > UCL_T_MPCA) |
    (results_df$max_Q_MPCA  > UCL_Q_MPCA)
)
results_df$SIGNAL_MPCA_MCD <- as.integer(
  (results_df$max_T2_MPCA_MCD > UCL_T_MPCA_MCD) |
    (results_df$max_Q_MPCA      > UCL_Q_MPCA_MCD)
)
results_df$SIGNAL_MPCA_MCD_80 <- as.integer(
  (results_df$max_T2_MPCA_MCD_80 > UCL_T_MPCA_MCD_80) |
    (results_df$max_Q_MPCA      > UCL_Q_MPCA_MCD_80)
)

# ------------------ MACRO-PARAFAC ------------------
results_df$SIGNAL_MACRO <- as.integer(
  (results_df$max_T2_MACRO > UCL_T_MACRO) |
    (results_df$max_Q_MACRO  > UCL_Q_MACRO)
)
results_df$SIGNAL_MACRO_MCD <- as.integer(
  (results_df$max_T2_MACRO_MCD > UCL_T_MACRO_MCD) |
    (results_df$max_Q_MACRO      > UCL_Q_MACRO_MCD)
)
results_df$SIGNAL_MACRO_MCD_80 <- as.integer(
  (results_df$max_T2_MACRO_MCD_80 > UCL_T_MACRO_MCD_80) |
    (results_df$max_Q_MACRO      > UCL_Q_MACRO_MCD_80)
)

# ============================================================================ #
#                           Tasas de señal promedio                            #
# ============================================================================ #

cat("Tasa de señal promedio (bajo control):\n")

mean_vals <- c(
  TROD        = mean(results_df$SIGNAL_TROD),
  TROD_MCD    = mean(results_df$SIGNAL_TROD_MCD),
  MPCA        = mean(results_df$SIGNAL_MPCA),
  MPCA_MCD    = mean(results_df$SIGNAL_MPCA_MCD),
  MPCA_MCD_80 = mean(results_df$SIGNAL_MPCA_MCD_80),
  MACRO       = mean(results_df$SIGNAL_MACRO),
  MACRO_MCD   = mean(results_df$SIGNAL_MACRO_MCD),
  MACRO_MCD_80= mean(results_df$SIGNAL_MACRO_MCD_80)
)

print(round(mean_vals, 3))




# ============================================================================ #
#           Crear los clusters para correr las diferentes simulaciones         #
# ============================================================================ #

num_cores <- detectCores() - 1  

# Crear el cluster
cl <- makeCluster(num_cores)

# Exportar los objetos necesarios al cluster
clusterExport(cl, varlist = c("run_simulation","generate_ellipse_image",
                              "generate_samples", "recon_trod_tensor",
                              "vecnorm", "makeDirections", "unimcd",
                              "MacroPARAFAC2"), envir = environment())

clusterEvalQ(cl, {
  library(grid)
  library(abind)
  library(rTensor)
  library(pracma)     
  library(stats)      
  library(cellWise)
  library(robustbase)
  library(rrcov)
  library(RobStatTM)
})

###################################################################
library(openxlsx)

# delta controles both delta5 and delta6
delta_values <- seq(1, 5, by = 1)  # Para delta 1
N_iter <- 1000

clusterSetRNGStream(cl, 520)

# DataFrame vacío para guardar todo
results_signal <- data.frame()

#########  Loop para diferenes valores de delta

for (delta_sim in delta_values) {
    resultados <- parLapply(
    cl, 1:N_iter,
    function(i, dval) {
      run_simulation(
        i,
        delta = 0,
        n_samples = 80,
        prop_in_control = 0.95,
        region_enable = FALSE,
        d1delta=dval
        noise_frac = 0.0)
    },
    dval = delta_sim
  )
  
  resultados_df2 <- do.call(rbind, resultados)
  
  señales_iter <- data.frame(
    delta = delta_sim,
    
    # ------------------ TROD ------------------
    señal_TROD      = as.integer(
      (resultados_df2$max_T2_TROD      > UCL_T_TROD) |
        (resultados_df2$max_Q_TROD       > UCL_Q_TROD)
    ),
    señal_TROD_MCD  = as.integer(
      (resultados_df2$max_T2_TROD_MCD  > UCL_T_TROD_MCD) |
        (resultados_df2$max_Q_TROD       > UCL_Q_TROD_MCD)
    ),
    
    # ------------------ MPCA ------------------
    señal_MPCA      = as.integer(
      (resultados_df2$max_T2_MPCA      > UCL_T_MPCA) |
        (resultados_df2$max_Q_MPCA       > UCL_Q_MPCA)
    ),
    señal_MPCA_MCD  = as.integer(
      (resultados_df2$max_T2_MPCA_MCD  > UCL_T_MPCA_MCD) |
        (resultados_df2$max_Q_MPCA       > UCL_Q_MPCA_MCD)
    ),
    señal_MPCA_MCD_80  = as.integer(
      (resultados_df2$max_T2_MPCA_MCD_80  > UCL_T_MPCA_MCD_80) |
        (resultados_df2$max_Q_MPCA       > UCL_Q_MPCA_MCD_80)
    ),
    
    # ------------------ MACRO-PARAFAC ------------------
    señal_MACRO     = as.integer(
      (resultados_df2$max_T2_MACRO     > UCL_T_MACRO) |
        (resultados_df2$max_Q_MACRO      > UCL_Q_MACRO)
    ),
    señal_MACRO_MCD = as.integer(
      (resultados_df2$max_T2_MACRO_MCD > UCL_T_MACRO_MCD) |
        (resultados_df2$max_Q_MACRO      > UCL_Q_MACRO_MCD)
    ),
    señal_MACRO_MCD_80 = as.integer(
      (resultados_df2$max_T2_MACRO_MCD_80 > UCL_T_MACRO_MCD_80) |
        (resultados_df2$max_Q_MACRO      > UCL_Q_MACRO_MCD_80)
    )
  )
  
  results_signal <- rbind(results_signal, señales_iter)
  cat("Fin de delta:", delta_sim, "\n")
}


# Calculamos tasa de detección por método y delta
tasa_senal <- aggregate(. ~ delta, data = results_signal, mean)
print(tasa_senal)



