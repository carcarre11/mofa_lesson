


#Cargamos los paquetes que necesitamos para la práctica

library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)
library(fitdistrplus)
library(pheatmap)
#cargamos los datos 
utils::data("CLL_data")       
#################################################
#################################################
#1Cargar y visualizar los datos

#Vamos a echar un vistazo a los datos
# esta formado por una lista de matrices

# respuesta a drogas exvivo
CLL_data$Drugs[1:10,1:10]

# Metilación 
CLL_data$Methylation[1:10,1:10]


# RNA 
CLL_data$mRNA[1:10,1:10] 

gene <- CLL_data$mRNA[1,]
gene <- gene[!is.na(gene)]

plotdist(gene, histo = TRUE, demp = TRUE)
plotdist(log2(gene), histo = TRUE, demp = TRUE)

# Mutations 
CLL_data$Mutations[1:10,1:10]

#Cargamos los metadatos
CLL_metadata <- fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/cll_vignette/sample_metadata.txt")

head(CLL_metadata)

# Creamos el objeto de MOFA
MOFAobject <- create_mofa(CLL_data)

#Vamos a ver el conjunto de nuestros datos. Vamos a ver la n de cada vista 
# y si hay algunos datos desconocidos
plot_data_overview(MOFAobject)

#################################################
#################################################
#2Generar el modelo
data_opts <- get_default_data_options(MOFAobject)
data_opts

# Dejamos las opciones de los datos como están.

# El escalado de las vistas permite ajustar
# el hecho de que haya una vista que aporte mucha mas información que otra
# Pero cuidado con la interpretación. Estas forzando a que las vistas
# sean comparables

model_opts <- get_default_model_options(MOFAobject)
model_opts

# En primer lugar cambiamos la distribución de
# las mutaciones que es de Bernouille (presencia ausencia)
model_opts$likelihoods["Mutations"] <- "bernoulli"
model_opts$num_factors <- 10


train_opts <- get_default_training_options(MOFAobject)

# los objetivos del training
# En esta práctica hemos seleccionado fast para
#que tarde menos pero recordad, que para vuestro modelo
#final necesita estar en slow para que sea mas preciso

#stochastic = TRUE puede ser util y necesario
# para dataset grandes como en single cell


train_opts$convergence_mode <- "fast"
train_opts$seed <- 123

# Entrenamos el modelo

MOFAobject <- prepare_mofa(MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

#Generamos el modelo con las opciones
#escogidas

MOFAobject <- run_mofa(MOFAobject, outfile="MOFA2_CLL.hdf5",use_basilisk = TRUE)

# Añadimos el metadata
samples_metadata(MOFAobject) <- CLL_metadata
#####################################################################################
#3 Investigar el modelo

# Buscamos factores que expliquen variabilidades únicas que no deben de correlacionar entre sí
plot_factor_cor(MOFAobject)

#Estudio de la varianza explicada por cada afctor en cada vista

plot_variance_explained(MOFAobject, max_r2=15)

# El Factor 1 captura una fuente de variabilidad presente en todas las modalidades de datos.  
# Por lo tanto, su etiología probablemente sea algo muy importante para la enfermedad.  

# El Factor 2 y 3 captura una fuente de variación presente en fármacos y mRNA.  

# El Factor 4 captura variaciones presentes en múltiples modalidades de datos, 
#mas escasa en la metilación del ADN.  
# Esto probablemente también sea importante.  

#Vamos a estudiar la variabildiad de que vistas estan mas explicadas por el modelo
plot_variance_explained(MOFAobject, plot_total = T)[[2]]
##############################################################################
#Asociación entre factores y variables clinicas con MOFA

correlate_factors_with_covariates(MOFAobject, 
  covariates = c("Gender","died","trisomy12","IGHV"), 
  plot="log_pval"
)

# Si no os convence los métodos usados por lo autores
# Sacad los pesos y usa tu propio código para hacer
# los análisis

Z <- get_expectations(MOFAobject, variable = "Z") 

meta <- MOFAobject@samples_metadata

meta <- cbind(meta,Z)

#Escribimos nuestras funciones para sacar la matriz de p valores usando  anova
source("src/functions.R")

#Sacamos la matriz de p-valores    
p.values <- matrix.pvalues(data = meta,
variables_of_interest = c("Gender","died","trisomy12", "IGHV"))

#Ajustamos los p valores por test multiple de cada variable
p_adjusted_matrix <- apply(p.values,  2, function(x) p.adjust(x, method = "fdr"))

#Transform to -log10(pval)
p_transformed_matrix <- -log10(p_adjusted_matrix)

# Convertir la matriz a un formato largo
p_matrix_long <- as.data.frame(as.table(p_transformed_matrix)) %>%
  rename(Row = Var1, Column = Var2, p_value = Freq) %>% 
  mutate(
    Row = gsub("group1\\.","",Row),
    Column = gsub("\\.p.value","",Column)
) 

p_matrix_long$p_value[p_matrix_long$p_value < -log10(0.05)] <- 0


# Crear el heatmap
ggplot(p_matrix_long, aes(x = Column, y = Row, fill = p_value)) +
  geom_tile(color = "black",lwd = 0.5) +
  scale_fill_gradient(low = "white", high = "red") +
  scale_fill_gradientn(colors = c("grey", "white", "red"), values = c(0, 0.000001, 1))+
  theme_minimal() +
  labs(title = "Heatmap de p-valores", x = "Columnas", y = "Filas") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


######################################################################################
# Una vez que sabemos que variables se relacionan con que factores podemos mostrarlo
#con diferentes gráficos

# En este caso el factor 1 y la mutaciñon de IGHV tiene una clara asociación

plot_factor(MOFAobject, 
  factors = 1, 
  color_by = "IGHV"
)


plot_factor(MOFAobject, 
  factors = 1, 
  color_by = "IGHV",
  add_violin = TRUE,
  dodge = TRUE
)
# Vamos a corroborar en la matriz de mutaciones que efectivamente IGHV esta explicando
# variabilidad en la vista de mutaciones

W.mutations <-  as.data.frame(get_expectations(MOFAobject, variable = "W")$Mutations)

head(W.mutations[order(W.mutations$Factor1, decreasing= T),])

#Podemos usar los plots del paquete para mostrarlo graficamente
plot_weights(MOFAobject,
 view = "Mutations",
 factor = 1,
 nfeatures = 10,     # Top number of features to highlight
 scale = T           # Scale weights from -1 to 1
)

plot_top_weights(MOFAobject,
 view = "Mutations",
 factor = 1,
 nfeatures = 10,     # Top number of features to highlight
 scale = T           # Scale weights from -1 to 1
)

#Ahora ya sabemos que el Factor1 basicamente diferencia las muestras en funciñon
# de su mutacion en IGHV.
#Ahora queremos saber como se relaciona esto con las demás capas

plot_weights(MOFAobject, 
  view = "mRNA", 
  factor = 1, 
  nfeatures = 10
)

#Vemos que hay un gen que claramente esta asociado a este factor al polo positivo
# y otro grupo hacia el otro sentido
# Hay un gen por lo tanto mas expresado en las muestras mutadas y otros genes
# menos expresados en las no mutadas.