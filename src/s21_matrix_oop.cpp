#include "s21_matrix_oop.h"

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::logic_error("Несоответствие размеров матриц");
  }

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::logic_error("Несоответствие размеров матриц");
  }

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols_ != other.rows_) {
    throw std::logic_error(
        "Столбцы первой матрицы не равны строкам второй матрицы!");
  }

  S21Matrix result(rows_, other.cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      for (int k = 0; k < cols_; k++) {
        result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
      }
    }
  }
  *this = result;
}

/*Аксесор*/
int S21Matrix::get_rows() { return rows_; }
int S21Matrix::get_cols() { return cols_; }

void S21Matrix::set_matrix(const std::vector<std::vector<double>>& values) {
  if ((int)values.size() != rows_) {
    throw std::invalid_argument("Несоответствие количества строк!");
  }
  if ((int)values[0].size() != cols_) {
    throw std::invalid_argument("Несоответствие количества столбцов!");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = values[i][j];
    }
  }
}

void S21Matrix::set_rows(int rows) {
  if (rows <= 0) {
    throw std::invalid_argument(
        "Количество строк должно больше нуля и положительным!\n");
  }
  if (rows != rows_) {
    S21Matrix tmp{rows, cols_};
    int min = std::min(rows_, rows);
    for (int i = 0; i < min; ++i) {
      for (int j = 0; j < cols_; ++j) {
        tmp(i, j) = (*this)(i, j);
      }
    }
    *this = std::move(tmp);
  }
}

void S21Matrix::set_cols(int cols) {
  if (cols <= 0) {
    throw std::invalid_argument(
        "Количество столбцов должно быть положительным и больше нуля!\n");
  }
  if (cols != cols_) {
    S21Matrix tmp{rows_, cols};
    int min = std::min(cols_, cols);
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < min; ++j) {
        tmp(i, j) = (*this)(i, j);
      }
    }

    *this = std::move(tmp);
  }
}

bool S21Matrix::EqMatrix(const S21Matrix& other) {
  bool flag = true;
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::logic_error("Несоответствие размеров");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      if (fabs(matrix_[i][j] - other.matrix_[i][j]) > EPSILON) {
        flag = false;
      }
    }
  }
  return flag;
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix res(cols_, rows_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      res.matrix_[j][i] = matrix_[i][j];
    }
  }
  return res;
}

S21Matrix S21Matrix::CalcComplements() {
  if (rows_ != cols_) {
    throw std::logic_error("Несоответствие размеров матриц");
  }

  S21Matrix result(rows_, cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      S21Matrix minor_matrix = get_minor(i, j);
      double minor_det = minor_matrix.Determinant();
      result.matrix_[i][j] = pow(-1, i + j) * minor_det;
    }
  }
  return result;
}

S21Matrix S21Matrix::get_minor(int row, int col) {
  S21Matrix minor(rows_ - 1, cols_ - 1);

  for (int i = 0, minor_i = 0; i < rows_; i++) {
    for (int j = 0, minor_j = 0; j < cols_; j++) {
      if (i != row && j != col) {
        minor.matrix_[minor_i][minor_j] = matrix_[i][j];
        minor_j++;
        if (minor_j == minor.cols_) {
          minor_j = 0;
          minor_i++;
        }
      }
    }
  }
  return minor;
}

double S21Matrix::Determinant() {
  if (rows_ != cols_) {
    throw std::logic_error("Несоответствие размеров матриц");
  }

  S21Matrix temp(*this);
  int sign = 1;
  double det = 1.0;

  for (int i = 0; i < rows_; i++) {
    // На шаге i если элемент на диагонали (temp.matrix_[i][i]) равен 0, то
    // меняем строки местами
    if (temp.matrix_[i][i] == 0) {
      int swap_row = i + 1;
      while (swap_row < temp.rows_ && temp.matrix_[swap_row][i] == 0) {
        swap_row++;  // Переходим к следующей строке
      }
      // если все элементы в строке равны 0, то матрица вырожденная
      if (swap_row == temp.rows_) {
        det = 0.0;
        break;
      }
      // меняем строки местами
      double* tmp = temp.matrix_[i];
      temp.matrix_[i] = temp.matrix_[swap_row];
      temp.matrix_[swap_row] = tmp;
      // Меняем знак определителя
      sign = -sign;
    }
    // приводим строку к верхнетреугольному виду
    for (int j = i + 1; j < temp.rows_; j++) {
      // вычисляем коэффициент
      double coeff = temp.matrix_[j][i] / temp.matrix_[i][i];
      // вычитаем коэффициент из строк
      for (int k = i; k < temp.cols_; k++) {
        temp.matrix_[j][k] -= coeff * temp.matrix_[i][k];
      }
    }

    det *= temp.matrix_[i][i];
  }
  return det * sign;
}

S21Matrix S21Matrix::InvertMatrix() {
  if (rows_ != cols_) {
    throw std::logic_error("Несоответствие размеров матриц");
  }
  double det = Determinant();
  if (fabs(det) < EPSILON) {
    throw std::logic_error("Матрица не имеет обратной");
  }
  S21Matrix invers(rows_, cols_);

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      S21Matrix minor = this->get_minor(i, j);
      double sign = ((i + j) % 2 == 0) ? 1.0 : -1.0;
      invers.matrix_[j][i] = sign * minor.Determinant() / det;
    }
  }
  return invers;
}

/*Вспомогательные функци*/
void S21Matrix::memory_allocation() {
  matrix_ = new double*[rows_];
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_];
  }
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (this == &other) {
    return *this;  // Защита от самоприсваивания
  }
  if (matrix_ != nullptr) {
    for (int i = 0; i < rows_; i++) {
      delete[] matrix_[i];
    }
    delete[] matrix_;
  }

  rows_ = other.rows_;
  cols_ = other.cols_;
  memory_allocation();
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
  return *this;
}

double& S21Matrix::operator()(int row, int col) {
  // Проверяем корректность индексов
  if (row < 0 || col < 0 || col >= cols_ || row >= rows_) {
    throw std::range_error("Segmentation fault!");
  }

  return matrix_[row][col];
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  S21Matrix res(rows_, cols_);
  res = *this;
  res.SumMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  S21Matrix res(rows_, cols_);
  res = *this;
  res.SubMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator*(const double num) {
  S21Matrix res(rows_, cols_);
  res = *this;
  res.MulNumber(num);
  return res;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  S21Matrix res(rows_, cols_);
  res = *this;
  res.MulMatrix(other);
  return res;
}

bool S21Matrix::operator==(const S21Matrix& other) { return EqMatrix(other); }
bool S21Matrix::operator!=(const S21Matrix& other) { return !EqMatrix(other); }

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  SumMatrix(other);
  return *this;
}
S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  SubMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  MulMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const double num) {
  MulNumber(num);
  return *this;
}