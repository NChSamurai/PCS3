#include <iostream>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <cstdlib>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <mpi.h>

double getTime() {
    auto now = std::chrono::high_resolution_clock::now();
    auto duration = now.time_since_epoch();
    double seconds = std::chrono::duration<double>(duration).count();
    seconds = std::round(seconds * 1e5) / 1e5;
    return seconds;
}

void appendTimeToBinaryFile(double timeValue) {
    std::ofstream outfile("GeneralTime.bin", std::ios::binary | std::ios::app);
    if (!outfile) {
        std::cerr << "Ошибка открытия файла GeneralTime.bin" << std::endl;
        return;
    }
    outfile.write(reinterpret_cast<const char*>(&timeValue), sizeof(double));
    outfile.close();
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const char* rows_env = std::getenv("ROWS");
    const char* cols_env = std::getenv("COLS");

    if (!rows_env || !cols_env) {
        if (rank == 0) {
            std::cerr << "Ошибка: переменные окружения ROWS и COLS должны быть установлены" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    const size_t ROWS = std::stoul(rows_env);
    const size_t COLS = std::stoul(cols_env);

    if (ROWS == 0 || COLS == 0) {
        if (rank == 0) {
            std::cerr << "Ошибка: размеры матрицы должны быть положительными числами" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    // Распределяем строки матрицы между процессами
    size_t rows_per_proc = ROWS / size;
    size_t remainder = ROWS % size;
    size_t local_rows = rows_per_proc + (rank < remainder ? 1 : 0);
    size_t offset = rank * rows_per_proc + (rank < remainder ? rank : remainder);

    // Выделяем память под локальные части матриц
    double** local_matrix1 = new double*[local_rows];
    double** local_matrix2 = new double*[local_rows];
    double** local_resultAdd = new double*[local_rows];
    double** local_resultSub = new double*[local_rows];
    double** local_resultMul = new double*[local_rows];
    double** local_resultDiv = new double*[local_rows];

    for (size_t i = 0; i < local_rows; ++i) {
        local_matrix1[i] = new double[COLS];
        local_matrix2[i] = new double[COLS];
        local_resultAdd[i] = new double[COLS];
        local_resultSub[i] = new double[COLS];
        local_resultMul[i] = new double[COLS];
        local_resultDiv[i] = new double[COLS];
    }

    // Процесс 0 читает полные матрицы и распределяет строки
    if (rank == 0) {
        double** matrix1 = new double*[ROWS];
        double** matrix2 = new double*[ROWS];
        for (size_t i = 0; i < ROWS; ++i) {
            matrix1[i] = new double[COLS];
            matrix2[i] = new double[COLS];
        }

        // Загрузка данных (в реальном коде нужно добавить обработку ошибок)
        std::ifstream file1("matrix.txt");
        std::ifstream file2("matrix2.txt");
        for (size_t i = 0; i < ROWS; ++i) {
            for (size_t j = 0; j < COLS; ++j) {
                file1 >> matrix1[i][j];
                file2 >> matrix2[i][j];
            }
        }

        // Отправляем данные другим процессам
        for (int dest = 1; dest < size; ++dest) {
            size_t dest_rows = rows_per_proc + (dest < remainder ? 1 : 0);
            size_t dest_offset = dest * rows_per_proc + (dest < remainder ? dest : remainder);

            for (size_t i = 0; i < dest_rows; ++i) {
                MPI_Send(matrix1[dest_offset + i], COLS, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
                MPI_Send(matrix2[dest_offset + i], COLS, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
            }
        }

        // Копируем свою часть
        for (size_t i = 0; i < local_rows; ++i) {
            std::copy(matrix1[i], matrix1[i] + COLS, local_matrix1[i]);
            std::copy(matrix2[i], matrix2[i] + COLS, local_matrix2[i]);
        }

        // Освобождаем память
        for (size_t i = 0; i < ROWS; ++i) {
            delete[] matrix1[i];
            delete[] matrix2[i];
        }
        delete[] matrix1;
        delete[] matrix2;
    } else {
        // Получаем данные от процесса 0
        for (size_t i = 0; i < local_rows; ++i) {
            MPI_Recv(local_matrix1[i], COLS, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(local_matrix2[i], COLS, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    double startTime = getTime();

    // Каждый процесс обрабатывает свои строки
    for (size_t i = 0; i < local_rows; ++i) {
        for (size_t j = 0; j < COLS; ++j) {
            local_resultAdd[i][j] = local_matrix1[i][j] + local_matrix2[i][j];
            local_resultSub[i][j] = local_matrix1[i][j] - local_matrix2[i][j];
            local_resultMul[i][j] = local_matrix1[i][j] * local_matrix2[i][j];
            local_resultDiv[i][j] = (local_matrix2[i][j] != 0) ? local_matrix1[i][j] / local_matrix2[i][j]
                                                              : std::numeric_limits<double>::quiet_NaN();
        }
    }

    double endTime = getTime();

    if (rank == 0) {
        std::cout << endTime - startTime << "\n";
        appendTimeToBinaryFile(endTime - startTime);
    }

    // Освобождаем память
    for (size_t i = 0; i < local_rows; ++i) {
        delete[] local_matrix1[i];
        delete[] local_matrix2[i];
        delete[] local_resultAdd[i];
        delete[] local_resultSub[i];
        delete[] local_resultMul[i];
        delete[] local_resultDiv[i];
    }
    delete[] local_matrix1;
    delete[] local_matrix2;
    delete[] local_resultAdd;
    delete[] local_resultSub;
    delete[] local_resultMul;
    delete[] local_resultDiv;

    MPI_Finalize();
    return 0;
}
