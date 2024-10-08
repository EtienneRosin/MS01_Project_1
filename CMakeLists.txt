cmake_minimum_required(VERSION 3.10)

# Nom du projet
project(MS01_PROJECT_1)

# Spécifier le standard C++ (par exemple, C++17)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Ajouter la recherche de MPI
find_package(MPI REQUIRED)

# Définir les répertoires de sources
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# Récupérer tous les fichiers .cpp dans les sous-dossiers utils et solvers
file(GLOB utils_srcs ${SOURCE_DIR}/utils/*.cpp)
file(GLOB solvers_srcs ${SOURCE_DIR}/solvers/*.cpp)

# Inclure les répertoires pour trouver les fichiers headers
include_directories(${SOURCE_DIR}/utils)
include_directories(${SOURCE_DIR}/solvers)
include_directories(${MPI_INCLUDE_PATH})  # Inclure les headers MPI

# Afficher les fichiers trouvés pour débogage
message(STATUS "Utils sources: ${utils_srcs}")
message(STATUS "Solvers sources: ${solvers_srcs}")

### 1. Main à la racine du projet ###
# Inclure le fichier main.cpp à la racine
file(GLOB main_src ${SOURCE_DIR}/main.cpp)

# Créer l'exécutable principal
add_executable(main ${main_src} ${utils_srcs} ${solvers_srcs})
target_link_libraries(main MPI::MPI_CXX)
target_compile_options(main PRIVATE -O2)

### 2. Cas d'étude 1 : jacobi_sequential_asymptotic_accuracy ###
# Ajouter un autre exécutable pour un cas d'étude spécifique
file(GLOB case1_srcs ${SOURCE_DIR}/study_cases/jacobi_sequential_asymptotic_accuracy/main.cpp)
add_executable(jacobi_sequential_asymptotic_accuracy ${case1_srcs} ${utils_srcs} ${solvers_srcs})
target_link_libraries(jacobi_sequential_asymptotic_accuracy MPI::MPI_CXX)
target_compile_options(jacobi_sequential_asymptotic_accuracy PRIVATE -O2)

# Cas d'étude 2 et 3 commentés pour l'instant
# file(GLOB case2_srcs ${SOURCE_DIR}/study_cases/jacobi_MPI_asymptotic_accuracy/main.cpp)
# add_executable(jacobi_MPI_asymptotic_accuracy ${case2_srcs} ${utils_srcs} ${solvers_srcs})
# target_link_libraries(jacobi_MPI_asymptotic_accuracy MPI::MPI_CXX)
# target_compile_options(jacobi_MPI_asymptotic_accuracy PRIVATE -O2)

# file(GLOB case3_srcs ${SOURCE_DIR}/study_cases/jacobi_sequential_MPI_comparison/main.cpp)
# add_executable(jacobi_sequential_MPI_comparison ${case3_srcs} ${utils_srcs} ${solvers_srcs})
# target_link_libraries(jacobi_sequential_MPI_comparison MPI::MPI_CXX)
# target_compile_options(jacobi_sequential_MPI_comparison PRIVATE -O2)

# Exporter les commandes de compilation (utile pour les outils comme clang-tidy)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")