# ==============================================================================
# PROYECTO: Análisis Automatizado de Cobertura Vegetal RGB
# AUTORA: Karla Cautivo
# DESCRIPCIÓN: Mitigación de luz mediante VARI y umbralización adaptativa (Otsu)
# ==============================================================================

# 1. Cargar librerías necesarias
if (!requireNamespace("terra", quietly = TRUE)) install.packages("terra")
library(terra)

# 2. Definir rutas del proyecto (Usa rutas relativas, ideales para GitHub)
carpeta_entrada <- "data/raw/"
carpeta_salida  <- "data/processed/"

# Listar todas las imágenes .tif recortadas de QGIS
imagenes_lista <- list.files(carpeta_entrada, pattern = "\\.tif$", full.names = TRUE)

# Crear un dataframe vacío para consolidar los resultados de la investigación
resultados_finales <- data.frame(Imagen = character(), Porcentaje_Vegetacion = numeric())

# 3. El Bucle Automatizado (Procesa cada imagen de forma independiente)
cat("--- Iniciando procesamiento de imágenes ---\n")

for (archivo in imagenes_lista) {
  
  nombre_archivo <- basename(archivo)
  cat("Procesando:", nombre_archivo, "...\n")
  
  # A. Cargar imagen y renombrar bandas a RGB estándar
  img <- rast(archivo)
  names(img) <- c("R", "G", "B")
  
  # B. Mitigar contaminación lumínica usando el índice normalizado VARI
  vari <- (img$G - img$R) / (img$G + img$R - img$B)
  
  # C. Umbralización de Otsu (Calcula el corte óptimo basado en la luz de ESTA foto)
  seg_binaria <- thresh(vari, method = "otsu")
  
  # D. Calcular porcentaje de cobertura de la máscara (Píxeles con valor 1)
  tabla_conteo <- freq(seg_binaria)
  total_pixeles      <- sum(tabla_conteo$count)
  pixeles_vegetacion <- tabla_conteo$count[tabla_conteo$value == 1]
  
  # Manejo de errores por si una imagen no tuviera vegetación
  if(length(pixeles_vegetacion) == 0) pixeles_vegetacion <- 0
  porcentaje_veg <- (pixeles_vegetacion / total_pixeles) * 100
  
  # E. Guardar los datos de esta imagen en nuestra tabla general
  resultados_finales <- rbind(resultados_finales, data.frame(
    Imagen = nombre_archivo, 
    Porcentaje_Vegetacion = round(porcentaje_veg, 2)
  ))
  
  # F. Guardar la máscara binaria final (1-0) en formato .tif procesado
  writeRaster(seg_binaria, paste0(carpeta_salida, "mascara_", nombre_archivo), overwrite = TRUE)
}

# 4. Exportar el reporte final de la investigación a un archivo CSV
write.csv(resultados_finales, paste0(carpeta_salida, "Reporte_Cobertura_Vegetal.csv"), row.names = FALSE)

cat("\n--- ¡PROCESO COMPLETADO EXITOSAMENTE! ---\n")
print(resultados_finales)
