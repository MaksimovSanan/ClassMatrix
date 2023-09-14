#ifndef CLASSMATRIX_CLASS_MATRIX_H
#define CLASSMATRIX_CLASS_MATRIX_H

#include <cmath>
#include <iostream>

namespace songrash {
class Matrix {
 public:
  explicit Matrix() noexcept;
  Matrix(int n, int m);
  Matrix(const Matrix& other);
  Matrix(Matrix&& other) noexcept;

  ~Matrix() noexcept;
  Matrix& operator=(const Matrix& other);
  Matrix& operator=(Matrix&& other) noexcept;

  const double& operator()(int row, int col) const;
  double& operator()(int row, int col);

  bool EqMatrix(const Matrix& other) const;
  bool operator==(const Matrix& other) const;
  void SumMatrix(const Matrix& other);
  Matrix operator+(const Matrix& other) const;
  Matrix& operator+=(const Matrix& other);
  void SubMatrix(const Matrix& other);
  Matrix operator-(const Matrix& other) const;
  Matrix& operator-=(const Matrix& other);
  void MulNumber(const double num) noexcept;
  Matrix operator*(const double num) const;
  friend Matrix operator*(const double num, const Matrix& other);
  Matrix& operator*=(const double num) noexcept;
  void MulMatrix(const Matrix& other);
  Matrix operator*(const Matrix& other) const;
  Matrix& operator*=(const Matrix& other);
  Matrix Transpose() const;
  Matrix CalcComplements() const;
  double Determinant() const;
  Matrix InverseMatrix() const;
  int GetRows() const noexcept;
  int GetCols() const noexcept;
  void SetRows(const int rows);
  void SetCols(const int cols);
  Matrix CreateMinor(const int i, const int j) const;
  friend std::ostream& operator<<(std::ostream& out, const Matrix& other);
  friend std::istream& operator>>(std::istream& in, Matrix& other);

 private:
  int rows_, cols_;
  double** matrix_;
};
}  // namespace songrash

#endif  // CLASSMATRIX_CLASS_MATRIX_H