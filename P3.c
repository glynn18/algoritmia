# Práctica 3: Programación Dinámica - Algoritmos Recursivos - Monte Carlo
# Algoritmia y Estructuras de Datos Avanzadas - Curso 2026-2027

import math
import random

# =============================================================================
# I. Alineamiento de Proteínas
# =============================================================================

# I.A Alineamiento Global: Needleman-Wunsch

def needleman_wunsch_matrix(seq1: str, seq2: str, match: int = 2, mismatch: int = -1, gap: int = -1) -> list[list[int]]:
    """
    Construye y devuelve la matriz completa de programación dinámica para el 
    alineamiento global de dos secuencias mediante el algoritmo Needleman-Wunsch.

    La matriz resultante tiene dimensiones (len(seq2) + 1) x (len(seq1) + 1), 
    donde seq1 se desarrolla a lo largo de las columnas y seq2 a lo largo de las filas.

    :param seq1: Primera secuencia a alinear (desarrollada en columnas).
    :param seq2: Segunda secuencia a alinear (desarrollada en filas).
    :param match: Puntuación acordada para una coincidencia de caracteres.
    :param mismatch: Penalización por diferencia de caracteres.
    :param gap: Penalización por introducción de hueco (gap).
    :return: Matriz de programación dinámica (lista de listas de enteros).
    """
    m = (len(seq2) + 1)
    n = (len(seq1) + 1)

    matriz = [[0 for j in range (n)] for i in range (m)]

    for i in range (0, m):
        matriz[i][0] = i*gap
    for j in range (0, n):
        matriz[0][j] = j*gap

    for i in range(1, m):
        for j in range(1, n):

            if seq1[j-1] == seq2[i-1]:
                puntuacion = match
            else:
                puntuacion = mismatch

            diagonal = matriz[i-1][j-1] + puntuacion
            arriba = matriz[i-1][j] + gap
            izquierda = matriz[i][j-1] + gap

            matriz[i][j] = max(diagonal, arriba, izquierda)

    return matriz     


def needleman_wunsch_backtracking(matrix: list[list[int]], seq1: str, seq2: str, match: int = 2, mismatch: int = -1, gap: int = -1) -> tuple[str, str, int]:
    """
    Reconstruye las secuencias alineadas globalmente y calcula el score total a partir 
    de la matriz de programación dinámica de Needleman-Wunsch mediante backtracking.

    El backtracking comienza en la esquina inferior derecha (len(seq2), len(seq1)) y 
    recorre la matriz hasta la esquina superior izquierda (0, 0). 
    En la estructura de la matriz, seq1 se desarrolla en columnas y seq2 en filas.

    :param matrix: Matriz de programación dinámica previamente calculada.
    :param seq1: Primera secuencia original (desarrollada en columnas).
    :param seq2: Segunda secuencia original (desarrollada en filas).
    :param match: Puntuación acordada para una coincidencia.
    :param mismatch: Penalización por diferencia.
    :param gap: Penalización por hueco.
    :return: Tupla con (aligned_seq1, aligned_seq2, total_score).
    """
    N = len(matrix)
    M = len(matrix[0])

    aligned_seq1 = ""
    aligned_seq2 = ""

    i = N - 1
    j = M - 1

    while i != 0 or j != 0:

        # OPCIÓN 1: Alineación / sustitución (diagonal)
        if i > 0 and j > 0:

            if seq1[j-1] == seq2[i-1]:
                puntuacion = match
            else:
                puntuacion = mismatch

            if matrix[i][j] == matrix[i-1][j-1] + puntuacion:
                aligned_seq1 += seq1[j-1]
                aligned_seq2 += seq2[i-1]

                i -= 1
                j -= 1
                continue

        # OPCIÓN 2: Hueco en Y (arriba)
        if i > 0:
            if matrix[i][j] == matrix[i-1][j] + gap:
                aligned_seq1 += "-"
                aligned_seq2 += seq2[i-1]

                i -= 1
                continue

        # OPCIÓN 3: Hueco en X (izquierda)
        if j > 0:
            if matrix[i][j] == matrix[i][j-1] + gap:
                aligned_seq1 += seq1[j-1]
                aligned_seq2 += "-"

                j -= 1
                continue

    aligned_seq1 = ''.join(reversed(aligned_seq1))
    aligned_seq2 = ''.join(reversed(aligned_seq2))

    total_score = matrix[N-1][M-1]

    return aligned_seq1, aligned_seq2, total_score


# I.B Alineamiento Local: Smith-Waterman

def smith_waterman_matrix(seq1: str, seq2: str, match: int = 2, mismatch: int = -1, gap: int = -1) -> list[list[int]]:
    """
    Construye y devuelve la matriz de programación dinámica para el alineamiento 
    local de dos secuencias mediante el algoritmo Smith-Waterman.

    Aplica inicialización a ceros en los bordes (primera fila y primera columna) 
    y truncamiento a cero en cada celda para evitar valores negativos. 
    La matriz tiene dimensiones (len(seq2) + 1) x (len(seq1) + 1), donde seq1 se 
    desarrolla a lo largo de las columnas y seq2 a lo largo de las filas.

    :param seq1: Primera secuencia a alinear (desarrollada en columnas).
    :param seq2: Segunda secuencia a alinear (desarrollada en filas).
    :param match: Puntuación por coincidencia.
    :param mismatch: Penalización por diferencia.
    :param gap: Penalización por hueco.
    :return: Matriz de programación dinámica local (lista de listas de enteros).
    """
    N = len(seq1)
    M = len(seq2)

    # Crear matriz inicializada a ceros
    dp = [[0 for j in range(M+1)] for i in range(N+1)]

    # Rellenar la matriz por recurrencia
    for i in range(1, N+1):
        for j in range(1, M+1):

            # OPCIÓN 1: Alineación / sustitución (diagonal)
            if seq1[i-1] == seq2[j-1]:
                puntuacion = match
            else:
                puntuacion = mismatch

            diagonal = dp[i-1][j-1] + puntuacion

            # OPCIÓN 2: Hueco en Y (arriba)
            hueco_y = dp[i-1][j] + gap

            # OPCIÓN 3: Hueco en X (izquierda)
            hueco_x = dp[i][j-1] + gap

            # Elegir el máximo, incluyendo el 0
            dp[i][j] = max(0, diagonal, hueco_y, hueco_x)

    return dp


def smith_waterman_backtracking(matrix: list[list[int]], seq1: str, seq2: str, match: int = 2, mismatch: int = -1, gap: int = -1) -> tuple[str, str, int]:
    """
    Reconstruye el alineamiento local óptimo y su puntuación a partir de la matriz 
    de Smith-Waterman mediante backtracking.

    Localiza el valor máximo global en la matriz y realiza el backtracking hasta 
    alcanzar una celda con valor 0. En la matriz, seq1 se desarrolla en columnas 
    y seq2 en filas.

    :param matrix: Matriz de programación dinámica local previamente calculada.
    :param seq1: Primera secuencia original (desarrollada en columnas).
    :param seq2: Segunda secuencia original (desarrollada en filas).
    :param match: Puntuación por coincidencia.
    :param mismatch: Penalización por diferencia.
    :param gap: Penalización por hueco.
    :return: Tupla con (sub_seq1, sub_seq2, max_score).
    """
     N = len(matrix)
    M = len(matrix[0])

    # Buscar el máximo global de la matriz
    max_score = 0
    i_max = 0
    j_max = 0

    for i in range(N):
        for j in range(M):
            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                i_max = i
                j_max = j

    # Inicializar el backtracking en el máximo
    i = i_max
    j = j_max

    aligned_seq1 = []
    aligned_seq2 = []

    # Retroceder hasta encontrar un cero
    while matrix[i][j] != 0:

        # OPCIÓN 1: Alineación / sustitución (diagonal)
        if i > 0 and j > 0:

            if seq1[i-1] == seq2[j-1]:
                puntuacion = match
            else:
                puntuacion = mismatch

            if matrix[i][j] == matrix[i-1][j-1] + puntuacion:

                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append(seq2[j-1])

                i -= 1
                j -= 1
                continue

        # OPCIÓN 2: Hueco en Y (arriba)
        if i > 0:
            if matrix[i][j] == matrix[i-1][j] + gap:

                aligned_seq1.append("-")
                aligned_seq2.append(seq1[i-1])

                i -= 1
                continue

        # OPCIÓN 3: Hueco en X (izquierda)
        if j > 0:
            if matrix[i][j] == matrix[i][j-1] + gap:

                aligned_seq1.append(seq2[j-1])
                aligned_seq2.append("-")

                j -= 1
                continue

    # Invertir las alineaciones
    aligned_seq1.reverse()
    aligned_seq2.reverse()

    # Convertir las listas en cadenas
    sub_seq1 = ''.join(aligned_seq1)
    sub_seq2 = ''.join(aligned_seq2)

    return sub_seq1, sub_seq2, max_score


# =============================================================================
# II. Optimización de Rutas en una Cuadrícula
# =============================================================================

# II.A Recursividad Ingenua (Enfoque de Fuerza Bruta)

def min_path_recursive(grid: list[list[int]], i: int, j: int) -> int:
    """
    Calcula el coste mínimo para alcanzar la celda (i, j) desde el origen (0, 0) 
    utilizando una estrategia recursiva pura (fuerza bruta sin memoización).

    :param grid: Matriz rectangular de costes numéricos enteros.
    :param i: Índice de la fila de la celda destino.
    :param j: Índice de la columna de la celda destino.
    :return: Coste acumulado mínimo para llegar a la celda (i, j).
    """

    #CASOS BASE
    #Para (0, 0), la función devolverá grid[0][0].
    #Si i == 0, solo podrá procederse desde la celda izquierda.
    #Si j == 0, solo podrá procederse desde la celda superior.    
    
    coste = 0
    if i==0 and j==0:
        return grid[0][0]
    elif i == 0:
        for k in range (0,j+1):
            coste += grid[0][k]
        return coste
    elif j == 0:
        for k in range (0,i+1):
            coste += grid[k][0]
        return coste
    else:  #FORMULA: F(i,j) = grid[i][j] + min(F(i−1,j), F(i,j−1))
        return grid[i][j] + min(min_path_recursive(grid, i-1, j), min_path_recursive(grid, i, j-1))   


# II.B Programación Dinámica (Enfoque Bottom-Up)

def min_path_dp(grid: list[list[int]]) -> tuple[int, list[list[int]]]:
    """
    Calcula el coste mínimo para ir desde la celda (0, 0) hasta (N-1, M-1) 
    mediante programación dinámica ascendente (bottom-up).

    :param grid: Matriz de costes de dimensiones N x M.
    :return: Tupla con el coste mínimo global (int) y la matriz dp de costes acumulados.
    """
    N = len(grid)
    M = len(grid[0])

    aux = [[0 for j in range(M)] for i in range(N)]

    for x in range (0, N):
        for y in range (0, M):
            if x==0 and y ==0:
                aux[x][y] = grid[x][y]
            elif x==0:
                aux[x][y] = grid[x][y] + aux[x][y-1]
            elif y==0:
                aux[x][y] = grid[x][y] + aux[x-1][y] 
            else:
                aux[x][y] = grid[x][y] + min(aux[x][y-1], aux[x-1][y])

    return aux[N-1][M-1], aux   


# II.C Reconstrucción de la Ruta Óptima

def get_optimal_path(grid: list[list[int]]) -> tuple[int, list[tuple[int, int]]]:
    """
    Reconstruye la secuencia exacta de coordenadas (0, 0) ... (N-1, M-1) que forman 
    el camino de coste mínimo en la cuadrícula.

    Invoca internamente a min_path_dp(grid) para generar la matriz de costes acumulados.

    :param grid: Matriz original de costes.
    :return: Tupla con el coste total mínimo y la lista de coordenadas (i, j) de la ruta.
    """
    (tot, aux) = min_path_dp(grid)

    N = len(grid)
    M = len(grid[0])

    pila = []

    i = N - 1
    j = M - 1

    pila.append((i, j))

    while i != 0 or j != 0:

        if i == 0:
            j -= 1
        elif j == 0:
            i -= 1
        elif aux[i][j-1] <= aux[i-1][j]:
            j -= 1
        else:
            i -= 1

        pila.append((i, j))

    pila.reverse()

    return tot, pila
  




# =============================================================================
# III. Simulación de Movimiento y Difusión mediante Caminata Aleatoria (Random Walk)
# =============================================================================

# III.A El comportamiento de una partícula individual

def simulate_single_walk(steps: int) -> tuple[float, float, float]:
    """
    Simula la trayectoria bidimensional de una sola partícula que parte desde el origen (0.0, 0.0) 
    y da un número determinado de pasos en direcciones aleatorias uniformes [0, 2*pi).

    :param steps: Número de pasos a simular.
    :return: Tupla con la posición final (x, y) y la distancia euclidiana final al origen.
    """
    pass

# III.B El enfoque de Monte Carlo y la Ley de Difusión

def monte_carlo_diffusion_simulation(num_particles: int, num_steps: int) -> float:
    """
    Aplica el método de Monte Carlo simulando un conjunto masivo de partículas independientes 
    para calcular la Distancia Cuadrática Media (R^2).

    :param num_particles: Número total de partículas independientes a simular.
    :param num_steps: Número de pasos por caminata.
    :return: Valor de la Distancia Cuadrática Media (R^2).
    """
    pass

# III.C Confinamiento y barreras absorbentes

def simulate_particle_with_barrier(steps: int, max_radius: float) -> bool:
    """
    Simula la caminata paso a paso de una partícula evaluando en cada paso si supera 
    un radio límite (barrera absorbente).

    :param steps: Número máximo de pasos.
    :param max_radius: Radio límite de la barrera desde el origen.
    :return: True si la partícula cruza la barrera en algún paso, False en caso contrario.
    """
    pass

def calculate_escape_probability(num_particles: int, num_steps: int, max_radius: float) -> float:
    """
    Estima la probabilidad de escape de un conjunto de partículas en presencia 
    de una barrera absorbente mediante el método de Monte Carlo.

    :param num_particles: Número de partículas a simular.
    :param num_steps: Número máximo de pasos por partícula.
    :param max_radius: Radio de la barrera absorbente.
    :return: Probabilidad de escape expresada en el intervalo [0.0, 1.0].
    """
    pass
