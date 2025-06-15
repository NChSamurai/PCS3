#include <iostream>
#include <fstream>
#include <stdexcept>
#include <limits>
#include <cstdlib>
#include <algorithm>
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
    if (!outfile.good()) {
        std::cerr << "Ошибка записи в файл" << std::endl;
    }
    outfile.close();
}

void loadArrayFromFile(const std::string& filename, double* arr, size_t size) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Не удалось открыть файл: " + filename);
    }

    for (size_t i = 0; i < size; ++i) {
        if (!(file >> arr[i])) {
            throw std::runtime_error("Ошибка чтения данных из файла: " + filename);
        }
    }
}

void printArray(const double* arr, size_t size) {
    const size_t print_count = std::min(size, static_cast<size_t>(10));
    for (size_t i = 0; i < print_count; ++i) {
        std::cout << arr[i] << " ";
    }
    if (size > 10) {
        std::cout << "...";
    }
    std::cout << std::endl;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Получаем размер из переменной окружения
    const char* size_env = std::getenv("SIZE");
    if (!size_env && rank == 0) {
        std::cerr << "Ошибка: переменная окружения SIZE не установлена" << std::endl;
        MPI_Finalize();
        return 1;
    }

    const size_t TOTAL_SIZE = std::stoul(size_env);
    if (TOTAL_SIZE <= 0 && rank == 0) {
        std::cerr << "Размер массива должен быть положительным числом" << std::endl;
        MPI_Finalize();
        return 1;
    }

    // Вычисляем размер части для каждого процесса
    size_t chunk_size = TOTAL_SIZE / size;
    size_t remainder = TOTAL_SIZE % size;

    // Распределяем остаток по процессам
    size_t local_size = chunk_size + (rank < remainder ? 1 : 0);
    size_t offset = rank * chunk_size + std::min(rank, (int)remainder);

    // Выделяем память для локальных данных
    double* local_array1 = new double[local_size];
    double* local_array2 = new double[local_size];
    double* local_resultAdd = new double[local_size];
    double* local_resultSub = new double[local_size];
    double* local_resultMul = new double[local_size];
    double* local_resultDiv = new double[local_size];

    // Главный процесс загружает данные и распределяет их
    if (rank == 0) {
        double* array1 = new double[TOTAL_SIZE];
        double* array2 = new double[TOTAL_SIZE];

        try {
            loadArrayFromFile("array.txt", array1, TOTAL_SIZE);
            loadArrayFromFile("array2.txt", array2, TOTAL_SIZE);
        } catch (const std::runtime_error& e) {
            std::cerr << e.what() << std::endl;
            delete[] array1;
            delete[] array2;
            MPI_Abort(MPI_COMM_WORLD, 1);
            MPI_Finalize();
            return 1;
        }

        // Отправляем данные всем процессам (включая себя)
        int* send_counts = new int[size];
        int* displacements = new int[size];

        for (int i = 0; i < size; ++i) {
            send_counts[i] = chunk_size + (i < remainder ? 1 : 0);
            displacements[i] = i * chunk_size + std::min(i, (int)remainder);
        }

        // Разослать данные всем процессам (используем MPI_IN_PLACE для root)
        MPI_Scatterv(array1, send_counts, displacements, MPI_DOUBLE,
                    local_array1, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Scatterv(array2, send_counts, displacements, MPI_DOUBLE,
                    local_array2, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        delete[] array1;
        delete[] array2;
        delete[] send_counts;
        delete[] displacements;
    } else {
        // Получаем данные от главного процесса
        MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE,
                    local_array1, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE,
                    local_array2, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
double startTime = getTime();
    // Синхронизируем перед началом вычислений
    MPI_Barrier(MPI_COMM_WORLD);


    // Локальные вычисления
    for (size_t i = 0; i < local_size; ++i) {
        local_resultAdd[i] = local_array1[i] + local_array2[i];
        local_resultSub[i] = local_array1[i] - local_array2[i];
        local_resultMul[i] = local_array1[i] * local_array2[i];
        if (local_array2[i] != 0) {
            local_resultDiv[i] = local_array1[i] / local_array2[i];
        } else {
            local_resultDiv[i] = std::numeric_limits<double>::quiet_NaN();
        }
    }

    // Синхронизируем перед окончанием вычислений
    MPI_Barrier(MPI_COMM_WORLD);
double endTime = getTime();
    // Только главный процесс выводит время и собирает результаты
    if (rank == 0) {
        std::cout << endTime - startTime << "\n";
        appendTimeToBinaryFile(endTime - startTime);

        // Выделяем память для полных результатов
        double* resultAdd = new double[TOTAL_SIZE];
        double* resultSub = new double[TOTAL_SIZE];
        double* resultMul = new double[TOTAL_SIZE];
        double* resultDiv = new double[TOTAL_SIZE];

        // Собираем результаты со всех процессов
        int* recv_counts = new int[size];
        int* displacements = new int[size];

        for (int i = 0; i < size; ++i) {
            recv_counts[i] = chunk_size + (i < remainder ? 1 : 0);
            displacements[i] = i * chunk_size + std::min(i, (int)remainder);
        }

        MPI_Gatherv(local_resultAdd, local_size, MPI_DOUBLE,
                   resultAdd, recv_counts, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Gatherv(local_resultSub, local_size, MPI_DOUBLE,
                   resultSub, recv_counts, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Gatherv(local_resultMul, local_size, MPI_DOUBLE,
                   resultMul, recv_counts, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Gatherv(local_resultDiv, local_size, MPI_DOUBLE,
                   resultDiv, recv_counts, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Здесь можно использовать полные результаты (resultAdd, resultSub и т.д.)
        // Например, вывести их или сохранить в файл

        delete[] resultAdd;
        delete[] resultSub;
        delete[] resultMul;
        delete[] resultDiv;
        delete[] recv_counts;
        delete[] displacements;
    } else {
        // Отправляем результаты главному процессу
        MPI_Gatherv(local_resultAdd, local_size, MPI_DOUBLE,
                   NULL, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Gatherv(local_resultSub, local_size, MPI_DOUBLE,
                   NULL, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Gatherv(local_resultMul, local_size, MPI_DOUBLE,
                   NULL, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Gatherv(local_resultDiv, local_size, MPI_DOUBLE,
                   NULL, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // Освобождаем локальную память
    delete[] local_array1;
    delete[] local_array2;
    delete[] local_resultAdd;
    delete[] local_resultSub;
    delete[] local_resultMul;
    delete[] local_resultDiv;

    MPI_Finalize();
    return 0;
}
