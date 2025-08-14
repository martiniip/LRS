# LRS - Longest Run Subsequence Problem

Este repositorio contiene implementaciones de varios algoritmos para resolver el problema de la Subsecuencia de Carrera Más Larga (Longest Run Subsequence).

## Descripción del Problema

El problema LRS consiste en encontrar una subsecuencia de una cadena de entrada que maximice la suma de las longitudes de las carreras (runs) consecutivas de caracteres idénticos.

## Estructura del Proyecto

```
├── src/                    # Código fuente de los algoritmos
├── Instances/             # Instancias de prueba
├── CD/                    # Análisis estadístico y resultados
├── LLM_instances/         # Instancias generadas por LLM
└── LLM_probabilities/     # Probabilidades LLM
```

## Algoritmos Implementados

El directorio `src/` contiene las siguientes implementaciones:

1. **aco.cpp** - Algoritmo de Optimización por Colonia de Hormigas (ACO)
2. **barrakuda.cpp** - Algoritmo Barrakuda
3. **brkga.cpp** - Algoritmo Genético con Claves Aleatorias Sesgadas (BRKGA)
4. **brkgallm.cpp** - BRKGA con componentes LLM
5. **cmsa.cpp** - Algoritmo CMSA (Construct, Merge, Solve & Adapt)
6. **cmsalrs.cpp** - CMSA específico para LRS
7. **cplex.cpp** - Implementación usando CPLEX
8. **grasp.cpp** - Algoritmo GRASP
9. **greedy.cpp** - Algoritmo Greedy
10. **greedy_aco.cpp** - Híbrido Greedy-ACO
11. **lcmsa_b.cpp** - LCMSA versión B

## Compilación

Para compilar todos los algoritmos, use los siguientes comandos:

```bash
# Compilar algoritmo greedy
g++ -std=c++11 -O3 -o greedy src/greedy.cpp

# Compilar ACO
g++ -std=c++11 -O3 -o aco src/aco.cpp

# Compilar BRKGA
g++ -std=c++11 -O3 -o brkga src/brkga.cpp

# Compilar GRASP
g++ -std=c++11 -O3 -o grasp src/grasp.cpp

# Compilar CMSA
g++ -std=c++11 -O3 -o cmsa src/cmsa.cpp

# Compilar CMSA-LRS
g++ -std=c++11 -O3 -o cmsalrs src/cmsalrs.cpp

# Compilar Barrakuda
g++ -std=c++11 -O3 -o barrakuda src/barrakuda.cpp

# Compilar BRKGA-LLM
g++ -std=c++11 -O3 -o brkgallm src/brkgallm.cpp

# Compilar Greedy-ACO
g++ -std=c++11 -O3 -o greedy_aco src/greedy_aco.cpp

# Compilar LCMSA-B
g++ -std=c++11 -O3 -o lcmsa_b src/lcmsa_b.cpp

# Para CPLEX (requiere CPLEX instalado)
g++ -std=c++11 -O3 -o cplex src/cplex.cpp -lcplex -lilocplex -lconcert
```

## Ejecución de Algoritmos

### Formato de Instancias

Las instancias están en formato de texto plano, cada archivo contiene una cadena de caracteres. Ejemplo:
```
BAABAABABBAABAABBABBBABAABBBBAAAAAABBAAABBBAABABBABBABBBBABAAABBBAAABAABBAABBABBBABAABBBAAAABAABABAB
```

### Parámetros Comunes

Todos los algoritmos aceptan el parámetro `-i` para especificar el archivo de instancia:

```bash
./algoritmo -i Instances/len_100_sigma2_1.txt
```

### Parámetros Específicos por Algoritmo

#### 1. Greedy
```bash
./greedy -i <archivo_instancia> [-pond <valor>] [-d <nivel_determinismo>]
```
- `-pond`: Factor de ponderación (default: valor por defecto)
- `-d`: Nivel de determinismo (default: valor por defecto)

#### 2. ACO (Ant Colony Optimization)
```bash
./aco -i <archivo_instancia> [-n_ants <num>] [-l_rate <tasa>] [-d_rate <tasa>] [-t <tiempo>] [-trials <num>] [-tuning] [-newheuristic]
```
- `-n_ants`: Número de hormigas (default: 50)
- `-l_rate`: Tasa de aprendizaje (default: 0.1)
- `-d_rate`: Tasa de determinismo (default: 0.2)
- `-t`: Límite de tiempo en segundos (default: 100.0)
- `-trials`: Número de ejecuciones (default: 1)
- `-tuning`: Activar modo tuning
- `-newheuristic`: Usar nueva heurística

#### 3. BRKGA (Biased Random-Key Genetic Algorithm)
```bash
./brkga -i <archivo_instancia> [-p <tamaño_pob>] [-pe <prop_elite>] [-pm <prop_mutante>] [-rhoe <prob_herencia>] [-t <tiempo>] [-th <umbral>] [-trials <num>] [-id <id>] [-tuning] [-stnweb]
```
- `-p`: Tamaño de población (default: valor por defecto)
- `-pe`: Proporción de élite (default: valor por defecto)
- `-pm`: Proporción de mutantes (default: valor por defecto)
- `-rhoe`: Probabilidad de herencia élite (default: valor por defecto)
- `-t`: Límite de tiempo (default: valor por defecto)
- `-th`: Umbral (default: valor por defecto)
- `-trials`: Número de ejecuciones
- `-id`: Identificador
- `-tuning`: Modo tuning
- `-stnweb`: Modo stnweb

#### 4. CMSA
```bash
./cmsa -i <archivo_instancia> [-nruns <num>] [-t <tiempo>] [-cpl_t <tiempo_cplex>] [-drate <tasa>] [-nsols <num>] [-max_age <edad>] [-lsize <tamaño>] [-warm_start <0|1>] [-h_emph <0|1>] [-cpl_abort <0|1>] [-tuning]
```
- `-nruns`: Número de ejecuciones
- `-t`: Límite de tiempo total
- `-cpl_t`: Límite de tiempo para CPLEX
- `-drate`: Tasa de determinismo
- `-nsols`: Número de soluciones
- `-max_age`: Límite de edad
- `-lsize`: Tamaño de lista candidata
- `-warm_start`: Inicio en caliente (0 o 1)
- `-h_emph`: Énfasis heurístico (0 o 1)
- `-cpl_abort`: Abortar CPLEX (0 o 1)
- `-tuning`: Modo tuning

### Ejemplos de Ejecución

```bash
# Ejecutar greedy simple
./greedy -i Instances/len_100_sigma2_1.txt

# Ejecutar ACO con parámetros personalizados
./aco -i Instances/len_100_sigma2_1.txt -n_ants 100 -t 60 -trials 5

# Ejecutar BRKGA con población grande
./brkga -i Instances/len_100_sigma2_1.txt -p 200 -t 120

# Ejecutar CMSA con tiempo limitado
./cmsa -i Instances/len_100_sigma2_1.txt -t 300 -cpl_t 10
```

## Ejecución en Lote

Para ejecutar todos los algoritmos en todas las instancias, puede crear un script de PowerShell:

```powershell
# Script para ejecutar todos los algoritmos
$algoritmos = @("greedy", "aco", "brkga", "grasp", "cmsa")
$instancias = Get-ChildItem -Path "Instances\*.txt"

foreach ($algoritmo in $algoritmos) {
    foreach ($instancia in $instancias) {
        Write-Host "Ejecutando $algoritmo en $($instancia.Name)"
        & ".\$algoritmo" -i $instancia.FullName
    }
}
```

## Formato de Salida

Los algoritmos generalmente producen salida que incluye:
- Información de la instancia (cadena de entrada, alfabeto, etc.)
- Parámetros utilizados
- Mejor solución encontrada
- Tiempo de ejecución
- Valor objetivo

## Análisis Estadístico

El directorio `CD/` contiene herramientas para análisis estadístico:

1. Prepare un archivo `results.dat` con los resultados de todos los algoritmos
2. Use R con el script `statistical_evaluation.R` para generar análisis estadísticos

Ver `CD/Readme` para instrucciones detalladas sobre el análisis estadístico.

## Instancias

El directorio `Instances/` contiene instancias de prueba con diferentes características:
- `len_100_sigma2_*.txt`: Cadenas de longitud 100 con alfabeto de tamaño 2
- `len_100_sigma16_*.txt`: Cadenas de longitud 100 con alfabeto de tamaño 16
- `len_100_sigma32_*.txt`: Cadenas de longitud 100 con alfabeto de tamaño 32

## Requisitos

- Compilador C++ con soporte para C++11 o superior
- Para CPLEX: IBM ILOG CPLEX Optimization Studio
- Para análisis estadístico: R con paquete 'scmamp'

## Notas

- Algunos algoritmos pueden requerir librerías adicionales
- Los tiempos de ejecución pueden variar significativamente según la instancia y los parámetros
- Para experimentos de investigación, se recomienda ejecutar múltiples trials y reportar estadísticas

## Contribución

Para añadir nuevos algoritmos:
1. Implemente siguiendo la estructura de parámetros estándar (`-i` para archivo de instancia)
2. Añada documentación de parámetros específicos
3. Actualice este README con instrucciones de compilación y ejecución
