#include <iostream>
#include <fstream>
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

    if (argc < 2) {
        if (rank == 0) {
            std::cerr << "Использование: " << argv[0] << " <размер_массива>" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    const int SIZE = atoi(argv[1]);
    if (SIZE <= 0) {
        if (rank == 0) {
            std::cerr << "Размер массива должен быть положительным числом" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    int* array = nullptr;
    int index = 0;
    double startTime = 0, endTime = 0;

    // Только процесс 0 читает файл
    if (rank == 0) {
        std::ifstream file("array.txt");
        if (!file.is_open()) {
            std::cerr << "Ошибка открытия файла :(" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        array = new int[SIZE];
        while (index < SIZE && file >> array[index]) {
            index++;
        }
        file.close();

        // Проверяем, что считали достаточно элементов
        if (index < SIZE) {
            std::cerr << "В файле недостаточно элементов (" << index << " из " << SIZE << ")" << std::endl;
            delete[] array;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        startTime = getTime();
    }

    // Рассылаем фактическое количество элементов (index)
    MPI_Bcast(&index, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Вычисляем размер части для каждого процесса
    int chunk_size = index / size;
    int remainder = index % size;

    // Определяем границы для каждого процесса
    int start = rank * chunk_size + (rank < remainder ? rank : remainder);
    int end = start + chunk_size + (rank < remainder ? 1 : 0);
    int local_size = end - start;

    // Выделяем память под локальную часть массива
    int* local_array = new int[local_size];

    // Разделяем массив между процессами
    MPI_Scatterv(array, nullptr, nullptr, MPI_INT,  // Эти параметры будут установлены на процессе 0
                 local_array, local_size, MPI_INT,
                 0, MPI_COMM_WORLD);

    // Локальное суммирование
    int local_sum = 0;
    for (int i = 0; i < local_size; i++) {
        local_sum += local_array[i];
    }

    // Собираем результаты на процессе 0
    int global_sum = 0;
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // Замер времени и вывод результата на процессе 0
    if (rank == 0) {
        endTime = getTime();
        std::cout << "Sum: " << global_sum << std::endl;
        std::cout << "Time: " << endTime - startTime << " seconds\n";
        appendTimeToBinaryFile(endTime - startTime);
        delete[] array;
    }

    delete[] local_array;
    MPI_Finalize();
    return 0;
}
