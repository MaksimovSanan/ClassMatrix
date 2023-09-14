#include "class_matrix.h"

namespace songrash {
Matrix::Matrix() noexcept {
  rows_ = 3;
  cols_ = 3;
  matrix_ = new double*[rows_]();
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_]();
  }
}

Matrix::Matrix(int n, int m) : rows_(n), cols_(m) {
  if (n <= 0 || m <= 0) {
    throw std::length_error("Incorrect input");
  }
  matrix_ = new double*[rows_]();
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_]();
  }
}

Matrix::Matrix(const Matrix& other)
    : rows_(other.rows_), cols_(other.cols_) {
  matrix_ = new double*[rows_]();
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_]();
  }

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

Matrix::Matrix(Matrix&& other) noexcept
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.matrix_ = nullptr;
  other.rows_ = 0;
  other.cols_ = 0;
}

Matrix::~Matrix() noexcept {
  for (int i = 0; i < rows_; i++) {
    delete matrix_[i];
  }
  delete matrix_;
  rows_ = 0;
  cols_ = 0;
}

Matrix& Matrix::operator=(const Matrix& other) {
  if (this == &other) return *this;

  for (int i = 0; i < rows_; i++) {
    delete matrix_[i];
  }
  delete matrix_;
  rows_ = other.rows_;
  cols_ = other.cols_;

  matrix_ = new double*[rows_]();
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_]();
  }

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
  return *this;
}

Matrix& Matrix::operator=(Matrix&& other) noexcept {
  if (this == &other) return *this;

  if (rows_ != other.rows_ || cols_ != other.cols_) {
    for (int i = 0; i < rows_; i++) {
      delete matrix_[i];
    }
    delete matrix_;
    matrix_ = nullptr;
    rows_ = other.rows_;
    cols_ = other.cols_;
  }

  matrix_ = other.matrix_;
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;

  return *this;
}

const double& Matrix::operator()(int row, int col) const {
  if (row < rows_ && col < cols_ && row < 0 && col < 0) {
    throw std::out_of_range("Incorrect input, index is out of range");
  }
  return matrix_[row][col];
}

double& Matrix::operator()(int row, int col) {
  if (row < rows_ && col < cols_ && row < 0 && col < 0) {
    throw std::out_of_range("Incorrect input, index is out of range");
  }
  return matrix_[row][col];
}

bool Matrix::EqMatrix(const Matrix& other) const {
  bool status = true;
  if (rows_ != other.rows_ || cols_ != other.cols_)
    status = false;
  else {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        if (std::fabs(matrix_[i][j] - other.matrix_[i][j]) > 1e-7)
          status = false;
      }
    }
  }
  return status;
}

bool Matrix::operator==(const Matrix& other) const {
  bool status = this->EqMatrix(other);
  return status;
}

void Matrix::SumMatrix(const Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::logic_error("Different matrix dimensions");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matrix_[i][j] + other.matrix_[i][j];
    }
  }
}

Matrix Matrix::operator+(const Matrix& other) const {
  Matrix res(*this);
  res.SumMatrix(other);
  return res;
}

Matrix& Matrix::operator+=(const Matrix& other) {
  this->SumMatrix(other);
  return *this;
}

void Matrix::SubMatrix(const Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::logic_error("Different matrix dimensions");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matrix_[i][j] - other.matrix_[i][j];
    }
  }
}

Matrix Matrix::operator-(const Matrix& other) const {
  Matrix res(*this);
  res.SubMatrix(other);
  return res;
}

Matrix& Matrix::operator-=(const Matrix& other) {
  this->SubMatrix(other);
  return *this;
}

void Matrix::MulNumber(const double num) noexcept {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matrix_[i][j] * num;
    }
  }
}

Matrix Matrix::operator*(const double num) const {
  Matrix res(*this);
  res.MulNumber(num);
  return res;
}

Matrix operator*(const double num, const Matrix& other) {
  Matrix res(other);
  res.MulNumber(num);
  return res;
}

Matrix& Matrix::operator*=(const double num) noexcept {
  this->MulNumber(num);
  return *this;
}

void Matrix::MulMatrix(const Matrix& other) {
  if (cols_ != other.rows_) {
    throw std::logic_error(
        "The number of columns of the first matrix does not equal the number "
        "of rows of the second matrix");
  }
  Matrix tmp(rows_, other.cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      double res = 0;
      for (int k = 0; k < cols_; k++) {
        res += (matrix_[i][k] * other.matrix_[k][j]);
      }
      tmp.matrix_[i][j] = res;
    }
  }
  *this = std::move(tmp);
}

Matrix Matrix::operator*(const Matrix& other) const {
  Matrix res(*this);
  res.MulMatrix(other);
  return res;
}

Matrix& Matrix::operator*=(const Matrix& other) {
  this->MulMatrix(other);
  return *this;
}

Matrix Matrix::Transpose() const {
  Matrix res(cols_, rows_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      res.matrix_[j][i] = matrix_[i][j];
    }
  }
  return res;
}

Matrix Matrix::CalcComplements() const {
  if (rows_ != cols_) {
    throw std::logic_error("the matrix is not square");
  }
  Matrix res(rows_, cols_);
  int correction_r = (rows_ > 1) ? 1 : 0;
  int correction_c = (cols_ > 1) ? 1 : 0;
  Matrix minor(rows_ - correction_r, cols_ - correction_c);

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      if (rows_ > 1)
        minor = this->CreateMinor(i, j);
      else
        minor.matrix_[i][j] = 1;
      double det = minor.Determinant();
      res.matrix_[i][j] = pow(-1, (i + j)) * det;
    }
  }
  return res;
}

double Matrix::Determinant() const {
  if (rows_ != cols_) {
    throw std::logic_error("the matrix is not square");
  }
  Matrix tmp(rows_, cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      tmp.matrix_[i][j] = matrix_[i][j];
    }
  }
  double result = 1;

  for (int i = 0; i < tmp.rows_; i++) {
    if (tmp.matrix_[i][i] == 0) {
      for (int j = i + 1; j < tmp.rows_; j++) {
        if (tmp.matrix_[j][i] != 0) {
          std::swap(tmp.matrix_[i], tmp.matrix_[j]);
          result *= -1;
        }
      }
    }
    if (tmp.matrix_[i][i] == 0) {
      result = 0;
      break;
    }
    result *= tmp.matrix_[i][i];

    double divider = tmp.matrix_[i][i];
    for (int j = i; j < tmp.cols_; j++) {
      tmp.matrix_[i][j] = tmp.matrix_[i][j] / divider;
    }

    for (int j = i + 1; j < tmp.rows_; j++) {
      double multipliyer = tmp.matrix_[j][i] / tmp.matrix_[i][i];
      for (int k = i; k < tmp.cols_; k++) {
        tmp.matrix_[j][k] -= tmp.matrix_[i][k] * multipliyer;
      }
    }
  }
  return result;
}

Matrix Matrix::InverseMatrix() const {
  double det = this->Determinant();

  if (det == 0) {
    throw std::logic_error("Matrix determinant is 0");
  }

  Matrix res(rows_, cols_);
  Matrix calc_comp = std::move(this->CalcComplements());
  Matrix calc_comp_trans = std::move(calc_comp.Transpose());

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      res.matrix_[i][j] = (1 / det) * calc_comp_trans.matrix_[i][j];
    }
  }

  return res;
}

int Matrix::GetRows() const noexcept { return rows_; }

int Matrix::GetCols() const noexcept { return cols_; }

void Matrix::SetRows(const int rows) {
  if (rows < 1) {
    throw std::out_of_range("Incorrect input, index is out of range");
  }
  Matrix tmp(rows, cols_);
  for (int i = 0; i < tmp.rows_; i++) {
    for (int j = 0; j < tmp.cols_; j++) {
      if (rows_ > i && cols_ > j) tmp.matrix_[i][j] = matrix_[i][j];
    }
  }
  *this = std::move(tmp);
}

void Matrix::SetCols(const int cols) {
  if (cols < 1) {
    throw std::out_of_range("Incorrect input, index is out of range");
  }
  Matrix tmp(rows_, cols);
  for (int i = 0; i < tmp.rows_; i++) {
    for (int j = 0; j < tmp.cols_; j++) {
      if (rows_ > i && cols_ > j) tmp.matrix_[i][j] = matrix_[i][j];
    }
  }
  *this = std::move(tmp);
}

Matrix Matrix::CreateMinor(const int i, const int j) const {
  Matrix res(cols_ - 1, rows_ - 1);
  for (int r = 0; r < rows_; r++) {
    for (int c = 0; c < cols_; c++) {
      if (r != i && c != j) {
        int correction_r = (r > i) ? -1 : 0;
        int correction_c = (c > j) ? -1 : 0;

        res.matrix_[r + correction_r][c + correction_c] = matrix_[r][c];
      }
    }
  }
  return res;
}

std::ostream& operator<<(std::ostream& out, const Matrix& other) {
  for (int i = 0; i < other.GetRows(); i++) {
    for (int j = 0; j < other.GetCols(); j++) {
      out << other.matrix_[i][j] << "\t";
    }
    out << std::endl;
  }
  return out;
}

std::istream& operator>>(std::istream& in, Matrix& other) {
  for (int i = 0; i < other.GetRows(); i++) {
    for (int j = 0; j < other.GetCols(); j++) {
      in >> other(i, j);
    }
  }
  return in;
}
}  // namespace songrash