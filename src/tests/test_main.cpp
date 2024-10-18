#include <gtest/gtest.h>

#include "../s21_matrix_oop.h"

TEST(CREATE_MATRIX, standart_constructor) {
  S21Matrix A;
  int rows = A.get_rows();
  int cols = A.get_cols();

  ASSERT_TRUE(rows == 3);
  ASSERT_TRUE(cols == 3);
}

TEST(CREATE_MATRIX, parameter_constructor) {
  S21Matrix A(5, 4);
  int rows = A.get_rows();
  int cols = A.get_cols();

  ASSERT_TRUE(rows == 5);
  ASSERT_TRUE(cols == 4);
}

TEST(CREATE_MATRIX, incorrect_parameters) {
  EXPECT_THROW(S21Matrix a(0, 0), std::range_error);
  EXPECT_THROW(S21Matrix a(-3, -3), std::range_error);
}

TEST(COPY_MATRIX, copy1) {
  S21Matrix A(2, 2);

  A(0, 0) = 1;
  A(0, 1) = 2;
  A(1, 0) = 3;
  A(1, 1) = 4;

  ASSERT_NO_THROW(S21Matrix err(A));
  S21Matrix B(A);
  ASSERT_TRUE(A.EqMatrix(B) == 1);
}

TEST(COPY_MATRIX, copy2) {
  S21Matrix A(2, 3);

  A(0, 0) = 1;
  A(0, 1) = 2;
  A(0, 2) = 3;
  A(1, 0) = 4;
  A(1, 1) = 5;
  A(1, 2) = 6;
  ASSERT_NO_THROW(S21Matrix err(A));
  S21Matrix B(A);
  ASSERT_TRUE(A.EqMatrix(B) == 1);
}

TEST(EQ_MATRIX, incorrect_parameters1) {
  S21Matrix A(2, 2);
  S21Matrix B(3, 4);

  EXPECT_THROW(A.EqMatrix(B), std::logic_error);
}

TEST(EQ_MATRIX, incorrect_parameters2) {
  S21Matrix A(2, 2);
  S21Matrix B(2, 2);

  A(0, 0) = 1;
  A(0, 1) = 2;
  A(1, 0) = 3;
  A(1, 1) = 4;

  B(0, 0) = 1;
  B(0, 1) = 2;
  B(1, 0) = 3;
  B(1, 1) = 5;

  ASSERT_FALSE(A.EqMatrix(B));
}

TEST(Transpose, test1) {
  S21Matrix A(2, 2);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(1, 0) = 3;
  A(1, 1) = 4;

  A.Transpose();

  S21Matrix B(2, 2);

  B(0, 0) = 1;
  B(0, 1) = 3;
  B(1, 0) = 2;
  B(1, 1) = 4;

  ASSERT_TRUE(A.Transpose() == B);
}

TEST(Determinant, test1) {
  S21Matrix A(2, 2);
  A.set_matrix({{1, 2}, {3, 4}});

  ASSERT_TRUE(A.Determinant() == -2);
}

TEST(Determinant, test2) {
  S21Matrix A(3, 3);
  A.set_matrix({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
  ASSERT_TRUE(A.Determinant() == 0);
}

TEST(Determinant, test3) {
  S21Matrix A(3, 3);
  A.set_matrix({{0, 0, 0}, {4, 5, 6}, {7, 8, 9}});
  ASSERT_TRUE(A.Determinant() == 0);
}
TEST(Determinant, test4) {
  S21Matrix A(3, 3);
  A.set_matrix({{0, 2, 3}, {0, 5, 6}, {7, 8, 9}});

  double det = A.Determinant();
  EXPECT_NEAR(det, -21.0, EPSILON);
}

TEST(Determinant, test6) {
  S21Matrix A(2, 2);
  A.set_matrix({{1, 2}, {3, 4}});
  ASSERT_TRUE(A.Determinant() == -2);
}

TEST(Determinant, test5) {
  S21Matrix A(2, 3);
  EXPECT_THROW(A.Determinant(), std::logic_error);
}

TEST(set_matrix, test1) {
  S21Matrix A(2, 2);
  EXPECT_THROW(A.set_matrix({{1, 2, 3}, {3, 4, 5}}), std::logic_error);
}

TEST(set_matrix, test2) {
  S21Matrix A(2, 2);
  EXPECT_THROW(A.set_matrix({{1, 2}, {3, 4}, {5, 6}}), std::logic_error);
}

TEST(set_cols, test1) {
  S21Matrix A(2, 2);
  A.set_matrix({{1, 2}, {3, 4}});
  A.set_cols(3);

  ASSERT_EQ(A.get_cols(), 3);
  ASSERT_EQ(A(0, 0), 1);
  ASSERT_EQ(A(0, 1), 2);
  ASSERT_EQ(A(1, 0), 3);
  ASSERT_EQ(A(1, 1), 4);
  ASSERT_EQ(A(0, 2), 0);  // Новый столбец должен быть инициализирован нулями
  ASSERT_EQ(A(1, 2), 0);
}

TEST(set_cols, test2) {
  S21Matrix A(2, 2);
  A.set_matrix({{1, 2}, {3, 4}});

  EXPECT_THROW(A.set_cols(0), std::invalid_argument);
}

TEST(set_rows, test1) {
  S21Matrix A(2, 2);
  A.set_matrix({{1, 2}, {3, 4}});
  A.set_rows(3);

  ASSERT_EQ(A.get_rows(), 3);
  ASSERT_EQ(A(0, 0), 1);
  ASSERT_EQ(A(0, 1), 2);
  ASSERT_EQ(A(1, 0), 3);
  ASSERT_EQ(A(1, 1), 4);
  ASSERT_EQ(A(2, 0), 0);
  ASSERT_EQ(A(2, 1), 0);
}

TEST(set_rows, test2) {
  S21Matrix A(2, 2);
  A.set_matrix({{1, 2}, {3, 4}});

  EXPECT_THROW(A.set_rows(0), std::invalid_argument);
}

TEST(Calc_matrix, incorrect_parameters) {
  S21Matrix A(4, 3);
  EXPECT_THROW(A.CalcComplements(), std::logic_error);
}

TEST(Calc_matrix, r2x2) {
  S21Matrix A(2, 2);
  A.set_matrix({{4, 7}, {2, 6}});
  S21Matrix B(2, 2);
  B.set_matrix({{6, -2}, {-7, 4}});
  ASSERT_TRUE(A.CalcComplements() == B);
}

TEST(Calc_matrix, test) {
  S21Matrix A(2, 2);
  A.set_matrix({{4, 3}, {2, 2}});
  S21Matrix B(2, 2);
  B.set_matrix({{2, -2}, {-3, 4}});
  ASSERT_TRUE(A.CalcComplements() == B);
}

TEST(INVERSE_MATRIX, ParametrWrong) {
  S21Matrix A(4, 3);
  EXPECT_THROW(A.InvertMatrix(), std::logic_error);
}

TEST(Inverse_matrix, r2x2) {
  S21Matrix A(2, 2);
  A.set_matrix({{4, 7}, {2, 6}});

  S21Matrix B(2, 2);
  B.set_matrix({{0.6, -0.7}, {-0.2, 0.4}});

  ASSERT_TRUE(A.InvertMatrix() == B);
}

TEST(Inverse_matrix, test) {
  S21Matrix A(3, 3);
  A.set_matrix({{5, 2, 1}, {2, 2, 2}, {4, 2, 4}});

  S21Matrix B(3, 3);
  B.set_matrix({{0.25, -0.375, 0.125}, {0, 1, -0.5}, {-0.25, -0.125, 0.375}});
  ASSERT_TRUE(A.InvertMatrix() == B);
}
TEST(REMOVE_MATRIX, remove1) {
  S21Matrix matrix_A(2, 2);
  ASSERT_NO_THROW(S21Matrix matrix_B = std::move(matrix_A));
}

TEST(REMOVE_MATRIX, remove2) {
  S21Matrix matrix_A(4, 3);
  ASSERT_NO_THROW(S21Matrix matrix_B = std::move(matrix_A));
}

TEST(OPERATORS, ij_1) {
  S21Matrix A(2, 2);

  A(0, 0) = 1;
  A(0, 1) = 2;
  A(1, 0) = 3;
  A(1, 1) = 4;

  ASSERT_TRUE(A(0, 0) == 1);
  ASSERT_TRUE(A(0, 1) == 2);
  ASSERT_TRUE(A(1, 0) == 3);
  ASSERT_TRUE(A(1, 1) == 4);
}

TEST(OPERATORS, ij_2) {
  S21Matrix A(2, 3);

  A(0, 0) = 1;
  A(0, 1) = 2;
  A(0, 2) = 3;
  A(1, 0) = 4;
  A(1, 1) = 5;
  A(1, 2) = 6;

  ASSERT_TRUE(A(0, 0) == 1);
  ASSERT_TRUE(A(0, 1) == 2);
  ASSERT_TRUE(A(0, 2) == 3);
  ASSERT_TRUE(A(1, 0) == 4);
  ASSERT_TRUE(A(1, 1) == 5);
  ASSERT_TRUE(A(1, 2) == 6);
}

TEST(OPERATORS, ij_3) {
  S21Matrix A(2, 2);

  A(0, 0) = 1;

  ASSERT_THROW(A(0, -1) = 1;, std::range_error);
}

TEST(OPERATORS, assignment1) {
  S21Matrix A(2, 2);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(1, 0) = 3;
  A(1, 1) = 4;

  S21Matrix B(2, 2);
  B = A;  // Присваивание

  ASSERT_EQ(B(0, 0), 1);
  ASSERT_EQ(B(0, 1), 2);
  ASSERT_EQ(B(1, 0), 3);
  ASSERT_EQ(B(1, 1), 4);
}

TEST(OPERATORS, assignment2) {
  S21Matrix A(2, 2);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(1, 0) = 3;
  A(1, 1) = 4;

  A = A;  // Присваивание

  ASSERT_EQ(A(0, 0), 1);
  ASSERT_EQ(A(0, 1), 2);
  ASSERT_EQ(A(1, 0), 3);
  ASSERT_EQ(A(1, 1), 4);
}

TEST(SUM_MATRIX, equal_matrix1) {
  S21Matrix A(2, 2);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(1, 0) = 3;
  A(1, 1) = 4;

  S21Matrix B(2, 2);
  B(0, 0) = 4;
  B(0, 1) = 3;
  B(1, 0) = 2;
  B(1, 1) = 1;

  A.SumMatrix(B);
  S21Matrix C(2, 2);
  C(0, 0) = 5;
  C(0, 1) = 5;
  C(1, 0) = 5;
  C(1, 1) = 5;
  ASSERT_TRUE(A.EqMatrix(C) == 1);
}
TEST(SUM_MATRIX, equal_matrix2) {
  S21Matrix A(2, 2);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(1, 0) = 3;
  A(1, 1) = 4;

  S21Matrix B(2, 2);
  B(0, 0) = 4;
  B(0, 1) = 3;
  B(1, 0) = 2;
  B(1, 1) = 1;

  S21Matrix C = A + B;
  ASSERT_EQ(C(0, 0), 5);
  ASSERT_EQ(C(0, 1), 5);
  ASSERT_EQ(C(1, 0), 5);
  ASSERT_EQ(C(1, 1), 5);
}

TEST(SUM_MATRIX, not_equal_matrix) {
  S21Matrix A(2, 2);
  S21Matrix B(3, 3);
  EXPECT_THROW(A.SumMatrix(B), std::logic_error);
}

TEST(SUB_MATRIX, equal_matrix1) {
  S21Matrix A(2, 2);
  A(0, 0) = 5;
  A(0, 1) = 5;
  A(1, 0) = 5;
  A(1, 1) = 5;

  S21Matrix B(2, 2);
  B(0, 0) = 4;
  B(0, 1) = 4;
  B(1, 0) = 4;
  B(1, 1) = 4;

  A.SubMatrix(B);
  S21Matrix C(2, 2);
  C(0, 0) = 1;
  C(0, 1) = 1;
  C(1, 0) = 1;
  C(1, 1) = 1;
  ASSERT_TRUE(A.EqMatrix(C) == 1);
}
TEST(SUB_MATRIX, equal_matrix2) {
  S21Matrix A(2, 2);
  A(0, 0) = 5;
  A(0, 1) = 5;
  A(1, 0) = 5;
  A(1, 1) = 5;

  S21Matrix B(2, 2);
  B(0, 0) = 4;
  B(0, 1) = 4;
  B(1, 0) = 4;
  B(1, 1) = 4;

  S21Matrix C = A - B;
  ASSERT_EQ(C(0, 0), 1);
  ASSERT_EQ(C(0, 1), 1);
  ASSERT_EQ(C(1, 0), 1);
  ASSERT_EQ(C(1, 1), 1);
}

TEST(SUB_MATRIX, not_equal_matrix) {
  S21Matrix A(2, 2);
  S21Matrix B(3, 3);

  EXPECT_THROW(A.SubMatrix(B), std::logic_error);
}

TEST(MULT_NUMBER, test1) {
  S21Matrix A(2, 2);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(1, 0) = 3;
  A(1, 1) = 4;

  A.MulNumber(2);

  S21Matrix B(2, 2);
  B(0, 0) = 2;
  B(0, 1) = 4;
  B(1, 0) = 6;
  B(1, 1) = 8;

  ASSERT_TRUE(A.EqMatrix(B) == 1);
}

TEST(MULT_NUMBER, test2) {
  S21Matrix A(2, 2);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(1, 0) = 3;
  A(1, 1) = 4;

  S21Matrix B(2, 2);
  B(0, 0) = 2;
  B(0, 1) = 4;
  B(1, 0) = 6;
  B(1, 1) = 8;

  S21Matrix C(2, 2);
  C = A * 2;
  ASSERT_TRUE(B.EqMatrix(C) == 1);
}

TEST(MULT_MATRIX, test1) {
  S21Matrix A(2, 2);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(1, 0) = 3;
  A(1, 1) = 4;

  S21Matrix B(2, 2);
  B(0, 0) = 4;
  B(0, 1) = 3;
  B(1, 0) = 2;
  B(1, 1) = 1;

  A.MulMatrix(B);

  S21Matrix C(2, 2);

  C(0, 0) = 8;
  C(0, 1) = 5;
  C(1, 0) = 20;
  C(1, 1) = 13;

  ASSERT_TRUE(C.EqMatrix(A) == true);
}
TEST(MULT_MATRIX, test2) {
  S21Matrix A(2, 2);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(1, 0) = 3;
  A(1, 1) = 4;

  S21Matrix B(2, 2);
  B(0, 0) = 4;
  B(0, 1) = 3;
  B(1, 0) = 2;
  B(1, 1) = 1;

  S21Matrix C(2, 2);
  C = A * B;
  S21Matrix res(2, 2);
  res(0, 0) = 8;
  res(0, 1) = 5;
  res(1, 0) = 20;
  res(1, 1) = 13;

  ASSERT_TRUE(res.EqMatrix(C) == 1);
}

TEST(MULT_MATRIX, incorrect_parameters) {
  S21Matrix A(2, 2);
  S21Matrix B(3, 3);
  EXPECT_THROW(A.MulMatrix(B), std::logic_error);
}

TEST(OPERATORS, Equal) {
  S21Matrix A(2, 2);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(1, 0) = 3;
  A(1, 1) = 4;

  S21Matrix B(2, 2);
  B(0, 0) = 1;
  B(0, 1) = 2;
  B(1, 0) = 3;
  B(1, 1) = 4;
  ASSERT_TRUE(A == B);
}

TEST(OPERATORS, NotEqual) {
  S21Matrix A(2, 2);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(1, 0) = 3;
  A(1, 1) = 4;

  S21Matrix B(2, 2);
  B(0, 0) = 4;
  B(0, 1) = 3;
  B(1, 0) = 2;
  B(1, 1) = 1;
  ASSERT_TRUE(A != B);
}

TEST(OPERATORS, PlusEqual) {
  S21Matrix A(2, 2);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(1, 0) = 3;
  A(1, 1) = 4;

  S21Matrix B(2, 2);
  B(0, 0) = 4;
  B(0, 1) = 3;
  B(1, 0) = 2;
  B(1, 1) = 1;
  S21Matrix C(2, 2);
  C(0, 0) = 5;
  C(0, 1) = 5;
  C(1, 0) = 5;
  C(1, 1) = 5;
  ASSERT_TRUE((A += B) == C);
}

TEST(OPERATORS, MinusEqual) {
  S21Matrix A(2, 2);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(1, 0) = 3;
  A(1, 1) = 4;

  S21Matrix B(2, 2);
  B(0, 0) = 4;
  B(0, 1) = 3;
  B(1, 0) = 2;
  B(1, 1) = 1;
  S21Matrix C(2, 2);
  C(0, 0) = -3;
  C(0, 1) = -1;
  C(1, 0) = 1;
  C(1, 1) = 3;
  ASSERT_TRUE((A -= B) == C);
}

TEST(OPERATORS, multiply_matrix) {
  S21Matrix A(2, 2);
  A(0, 0) = 1;
  A(0, 1) = 1;
  A(1, 0) = 1;
  A(1, 1) = 1;

  S21Matrix B(2, 2);
  B(0, 0) = 1;
  B(0, 1) = 1;
  B(1, 0) = 1;
  B(1, 1) = 1;
  S21Matrix C(2, 2);
  C(0, 0) = 2;
  C(0, 1) = 2;
  C(1, 0) = 2;
  C(1, 1) = 2;
  ASSERT_TRUE((A *= B) == C);
}

TEST(OPERATORS, multiply_num) {
  S21Matrix A(2, 2);
  A(0, 0) = 1;
  A(0, 1) = 1;
  A(1, 0) = 1;
  A(1, 1) = 1;

  int B = 6;
  S21Matrix C(2, 2);
  C(0, 0) = 6;
  C(0, 1) = 6;
  C(1, 0) = 6;
  C(1, 1) = 6;
  ASSERT_TRUE((A *= B) == C);
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}