#ifndef __S21_MATRIX_OOP_H__
#define __S21_MATRIX_OOP_H__

#include <math.h>

#include <iostream>
#include <stdexcept>
#include <vector>

#define EPSILON 1e-7

class S21Matrix {
 private:
  // Attributes
  int rows_, cols_;  // Rows and columns
  double** matrix_;  // Pointer to the memory where the matrix is allocated

 public:
  // Default constructor
  S21Matrix() : S21Matrix(3, 3){};

  // Parameterized constructor
  S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
    if (rows < 1 || cols < 1) {
      throw std::range_error("Строки и столбцы дожны иметь значения больше 0");
    } else {
      memory_allocation();
      for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
          matrix_[i][j] = 0;
        }
      }
    }
  }
  // Copy constructor
  S21Matrix(const S21Matrix& other) : rows_(other.rows_), cols_(other.cols_) {
    memory_allocation();
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        matrix_[i][j] = other.matrix_[i][j];
      }
    }
  }
  // Move constructor
  S21Matrix(S21Matrix&& other) : rows_(other.rows_), cols_(other.cols_) {
    matrix_ = other.matrix_;
    other.matrix_ = nullptr;
    other.rows_ = 0;
    other.cols_ = 0;
  }
  // Destructor
  ~S21Matrix() {
    for (int i = 0; i < rows_; i++) {
      delete[] matrix_[i];
    }
    delete[] matrix_;
    matrix_ = nullptr;
    rows_ = 0;
    cols_ = 0;
  }

  // Operators
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);

  S21Matrix Transpose();
  S21Matrix CalcComplements();
  S21Matrix InvertMatrix();

  bool EqMatrix(const S21Matrix& other);
  double Determinant();

  /*Перегрузка операторов*/
  S21Matrix& operator=(const S21Matrix& other);
  double& operator()(int row, int col);
  S21Matrix operator+(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix operator*(const double num);
  S21Matrix operator*(const S21Matrix& other);
  bool operator==(const S21Matrix& other);
  bool operator!=(const S21Matrix& other);
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix& operator*=(const S21Matrix& other);
  S21Matrix& operator*=(const double num);

  /*Аксесор*/
  int get_rows();
  int get_cols();

  void set_rows(int rows);
  void set_cols(int cols);

  /*Вспомогательные функции*/
  void memory_allocation();
  void set_matrix(const std::vector<std::vector<double>>& values);
  S21Matrix get_minor(int rows, int cols);
};

#endif