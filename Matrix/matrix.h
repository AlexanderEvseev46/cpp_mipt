#include <array>
#include <cmath>
#include <type_traits>

#define CPP23 1

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

class BigInteger {
  friend std::istream& operator>>(std::istream& in, BigInteger& bi);

  friend BigInteger abs(const BigInteger&);

 private:
  enum class Sign { Negative = -1, Zero = 0, Positive = 1 };

  const static int BASE_ = 1000000000;

  const static size_t RANK_SIZE = 9;

  Sign sign_ = Sign::Positive;

  std::vector<int> data_;

  void delete_nulls();

  void decreaseBiggest(const BigInteger& other);

  int Binsearch(const BigInteger& other);

  void turnSign();

  void setSign(int value);

  void takeBase(int index);

  void getBase(int index);

 public:
  BigInteger();

  BigInteger(const std::string& s);

  BigInteger(int value);

  const std::vector<int>& data() const;

  int sign() const;

  BigInteger operator-() const;

  BigInteger& operator+=(const BigInteger& other);

  BigInteger& operator-=(const BigInteger& other);

  BigInteger& operator*=(const int& other);

  BigInteger& operator*=(const BigInteger& other);

  BigInteger& operator<<(const int& rank);

  BigInteger& operator/=(const BigInteger& other);

  BigInteger& operator%=(const BigInteger& other);

  BigInteger& operator++();

  BigInteger operator++(int);

  BigInteger& operator--();

  BigInteger operator--(int);

  explicit operator bool();

  std::string toString() const;
};

void BigInteger::turnSign() {
  if (sign_ == Sign::Positive) {
    sign_ = Sign::Negative;
    return;
  }

  if (sign_ == Sign::Negative) {
    sign_ = Sign::Positive;
  }
}

void BigInteger::setSign(int value) {
  if (value == 0) {
    sign_ = BigInteger::Sign::Zero;
    return;
  }
  if (value < 0) {
    sign_ = BigInteger::Sign::Negative;
    return;
  }
  if (value > 0) {
    sign_ = BigInteger::Sign::Positive;
  }
}

BigInteger::BigInteger() {
  data_.resize(1, 0);
  sign_ = Sign::Positive;
}

BigInteger::BigInteger(const std::string& s) {
  if (s[0] == '0') {
    sign_ = Sign::Zero;
  }
  if (s[0] != '0' && s[0] != '-') {
    sign_ = Sign::Positive;
  }
  if (s[0] == '-') {
    sign_ = Sign::Negative;
    data_.resize((s.length() - 2) / RANK_SIZE + 1);
  } else {
    data_.resize((s.length() - 1) / RANK_SIZE + 1);
  }
  int i = s.length() - 9;
  size_t k = 0;
  while (i >= 0) {
    data_[k] = stoi(s.substr(i, RANK_SIZE));
    ++k;
    i -= 9;
  }
  if (k < data_.size()) {
    data_[k] = abs(stoi(s.substr(0, i + 9)));
  }
}

BigInteger::BigInteger(int value) {
  data_.resize(1, 0);
  data_[0] = (std::abs(value) % BASE_);
  if (value >= BASE_) {
    data_.resize(2, 0);
    data_[1] = (std::abs(value) / BASE_);
  }
  setSign(value);
}

void BigInteger::delete_nulls() {
  int i = data_.size() - 1;
  size_t cnt_nonnulls = data_.size();
  while (i >= 0 && data_[i] == 0) {
    --cnt_nonnulls;
    --i;
  }
  if (cnt_nonnulls == 0) {
    *this = 0;
    return;
  }
  if (cnt_nonnulls != data_.size()) {
    data_.resize(cnt_nonnulls);
  }
}

const std::vector<int>& BigInteger::data() const { return data_; }

int BigInteger::sign() const { return static_cast<int>(sign_); }

BigInteger abs(const BigInteger& bi) {
  BigInteger result = bi;
  if (result.sign_ != BigInteger::Sign::Zero) {
    result.sign_ = BigInteger::Sign::Positive;
  }

  return result;
}

BigInteger BigInteger::operator-() const {
  BigInteger result = *this;
  result.turnSign();
  return result;
}

bool module_compare(const BigInteger& first, const BigInteger& second) {
  if (first.data().size() > second.data().size()) {
    return false;
  }

  if (first.data().size() < second.data().size()) {
    return true;
  }

  for (int i = first.data().size() - 1; i >= 0; --i) {
    if (first.data()[i] < second.data()[i]) {
      return true;
    }
    if (first.data()[i] > second.data()[i]) {
      return false;
    }
  }
  return false;
}

bool operator<(const BigInteger& first, const BigInteger& second) {
  if (first.sign() < second.sign()) {
    return true;
  }

  if (first.sign() > second.sign()) {
    return false;
  }

  if (first.sign() == 1) {
    return module_compare(first, second);
  }

  if (first.sign() == -1) {
    return !module_compare(first, second);
  }

  // if (sign == 0)
  return false;
}

bool operator>(const BigInteger& first, const BigInteger& second) {
  return second < first;
}

bool operator==(const BigInteger& first, const BigInteger& second) {
  return (!(first < second)) && (!(first > second));
}

bool operator!=(const BigInteger& first, const BigInteger& second) {
  return !(first == second);
}

bool operator<=(const BigInteger& first, const BigInteger& second) {
  return !(first > second);
}

bool operator>=(const BigInteger& first, const BigInteger& second) {
  return !(first < second);
}

void BigInteger::takeBase(int index) {
  data_[index] += BASE_;
  --data_[index + 1];
}

void BigInteger::getBase(int index) {
  data_[index + 1] += data_[index] / BASE_;
  data_[index] %= BASE_;
}

void BigInteger::decreaseBiggest(const BigInteger& other) {
  if (other.sign_ == Sign::Zero) {
    return;
  }

  if (sign_ == other.sign_) {
    for (size_t i = 0; i < other.data_.size() - 1; ++i) {
      data_[i] += other.data_[i];
      getBase(i);
    }
    data_[other.data_.size() - 1] += other.data_[other.data_.size() - 1];
    if (data_[other.data_.size() - 1] >= BASE_) {
      if (data_.size() == other.data_.size()) {
        data_.resize(data_.size() + 1, 0);
      }
      getBase(other.data_.size() - 1);
    }
    return;
  }

  for (size_t i = 0; i < other.data_.size(); ++i) {
    if (other.data_[i] > data_[i]) {
      takeBase(i);
    }
    data_[i] -= other.data_[i];
  }

  for (size_t i = other.data_.size(); i < data_.size(); ++i) {
    if (data_[i] < 0) {
      takeBase(i);
    } else {
      break;
    }
  }

  delete_nulls();
  return;
}

BigInteger& BigInteger::operator+=(const BigInteger& other) {
  if (module_compare(*this, other)) {
    BigInteger bi = other;
    std::swap(*this, bi);
    decreaseBiggest(bi);
    return *this;
  }

  decreaseBiggest(other);
  return *this;
}

BigInteger operator+(const BigInteger& first, const BigInteger& second) {
  BigInteger answer = first;
  answer += second;
  return answer;
}

BigInteger& BigInteger::operator-=(const BigInteger& other) {
  *this += (-other);
  return *this;
}

BigInteger operator-(const BigInteger& first, const BigInteger& second) {
  BigInteger answer = first;
  answer -= second;
  return answer;
}

BigInteger& BigInteger::operator*=(const int& other) {
  unsigned long long other_module = std::abs(other);
  int next = 0;
  size_t sz = data_.size();
  for (size_t i = 0; i < sz; ++i) {
    long long composition = other_module * (data_[i] - next) + next;
    data_[i] = composition % BASE_;
    next = composition / BASE_;
    if (i + 1 == data_.size() && next > 0) {
      data_.resize(data_.size() + 1, 0);
    }
    if (next > 0) {
      data_[i + 1] += next;
    }
  }
  if (other == 0) {
    sign_ = Sign::Zero;
  }
  if (other < 0) {
    turnSign();
  }
  delete_nulls();
  return *this;
}

BigInteger& BigInteger::operator*=(const BigInteger& other) {
  BigInteger answer = 0;
  for (size_t i = 0; i < other.data_.size(); ++i) {
    BigInteger part = *this;
    part *= other.data_[i];
    part << static_cast<int>(i);
    answer += part;
  }

  *this = answer;
  if (other.sign_ == Sign::Negative) {
    turnSign();
  }
  if (other.sign_ == Sign::Zero) {
    sign_ = Sign::Zero;
  }

  return *this;
}

BigInteger operator*(const BigInteger& first, const BigInteger& second) {
  BigInteger result = first;
  result *= second;
  return result;
}

BigInteger& BigInteger::operator<<(const int& rank) {
  if (sign_ == Sign::Zero) {
    return *this;
  }
  data_.resize(data_.size() + rank);
  for (int i = data_.size() - rank - 1; i >= 0; --i) {
    data_[rank + i] = data_[i];
  }
  for (int i = 0; i < rank; ++i) {
    data_[i] = 0;
  }
  return *this;
};

int BigInteger::Binsearch(const BigInteger& other) {
  int left = 0;
  int right = BASE_;
  while (right - left > 1) {
    int mid = (left + right) / 2;
    if (!module_compare(*this, other * mid)) {
      left = mid;
    } else {
      right = mid;
    }
  }
  return left;
}

BigInteger& BigInteger::operator/=(const BigInteger& other) {
  if (module_compare(*this, other)) {
    *this = 0;
    return *this;
  }

  BigInteger result = 0;
  BigInteger current = 0;
  for (int i = data_.size() - 1; i >= 0; --i) {
    if (module_compare(current, other)) {
      current << 1;
      current += data_[i];
    }
    if (module_compare(current, other)) {
      result << 1;
      continue;
    }
    result << 1;
    BigInteger part = current.Binsearch(other);
    result += part;
    BigInteger x = (part * other * other.sign());
    current -= x;
  }

  if ((sign_ == Sign::Negative && other.sign_ == Sign::Positive) ||
      (sign_ == Sign::Positive && other.sign_ == Sign::Negative)) {
    result.turnSign();
  }

  *this = result;
  delete_nulls();
  return *this;
}

BigInteger operator/(const BigInteger& first, const BigInteger& second) {
  BigInteger result = first;
  result /= second;
  return result;
}

BigInteger& BigInteger::operator%=(const BigInteger& other) {
  *this -= (*this / other) * other;
  delete_nulls();
  return *this;
}

BigInteger operator%(const BigInteger& first, const BigInteger& second) {
  BigInteger result = first;
  result %= second;
  return result;
}

BigInteger& BigInteger::operator++() {
  *this += 1;
  return *this;
}

BigInteger BigInteger::operator++(int) {
  BigInteger result = *this;
  *this += 1;
  return result;
}

BigInteger& BigInteger::operator--() {
  *this -= 1;
  return *this;
}

BigInteger BigInteger::operator--(int) {
  BigInteger result = *this;
  *this -= 1;
  return result;
}

BigInteger::operator bool() {
  if (data_.size() == 1 && data_[0] == 0) {
    return false;
  }
  return true;
}

std::string BigInteger::toString() const {
  if (sign_ == Sign::Zero) {
    return "0";
  }

  std::string s;
  if (sign_ == Sign::Negative) {
    s = '-';
  }
  s += std::to_string(data_[data_.size() - 1]);
  for (int i = data_.size() - 2; i >= 0; --i) {
    std::string part = "000000000" + std::to_string(data_[i]);
    s += part.substr(part.length() - 9, 9);
  }
  return s;
}

BigInteger operator"" _bi(unsigned long long x) {
  BigInteger result = std::to_string(x);
  return result;
}

std::ostream& operator<<(std::ostream& out, const BigInteger& bi) {
  out << bi.toString();
  return out;
}

std::istream& operator>>(std::istream& in, BigInteger& bi) {
  std::string s;
  in >> s;
  BigInteger x(s);
  bi = x;
  return in;
}

class Rational {
  friend bool operator<(const Rational&, const Rational&);

 private:
  BigInteger numerator_;
  BigInteger denominator_;

 public:
  Rational(BigInteger x) : numerator_(x), denominator_(1) {}

  Rational(int x) : numerator_(x), denominator_(1) {}

  Rational() : numerator_(0), denominator_(1) {}

  Rational operator-() const {
    Rational bi = *this;
    bi.numerator_ = -numerator_;
    return bi;
  }

  void divide_on_gsd();

  Rational& operator+=(const Rational& other);

  Rational& operator-=(const Rational& other);

  Rational& operator*=(const Rational& other);

  Rational& operator/=(const Rational& other);

  std::string toString() const;

  std::string asDecimal(size_t precision) const;

  explicit operator double() const;
};

BigInteger GCD(BigInteger first, BigInteger second) {
  if (first == 0) {
    return second;
  }
  if (second == 1) {
    return 1;
  }
  first *= first.sign();
  second *= second.sign();
  while (first > 0 && second > 0) {
    if (first > second) {
      first %= second;
    } else {
      second %= first;
    }
  }
  if (first == 0) {
    return second;
  }
  return first;
}

void Rational::divide_on_gsd() {
  BigInteger gsd = GCD(numerator_, denominator_);
  numerator_ /= gsd;
  denominator_ /= gsd;
}

Rational& Rational::operator+=(const Rational& other) {
  numerator_ =
      numerator_ * other.denominator_ + denominator_ * other.numerator_;
  denominator_ *= other.denominator_;
  divide_on_gsd();
  return *this;
}

Rational& Rational::operator-=(const Rational& other) {
  *this += (-other);
  return *this;
}

Rational& Rational::operator*=(const Rational& other) {
  numerator_ *= other.numerator_;
  denominator_ *= other.denominator_;
  divide_on_gsd();
  return *this;
}

Rational& Rational::operator/=(const Rational& other) {
  numerator_ *= other.denominator_;
  numerator_ *= other.numerator_.sign();
  denominator_ *= other.numerator_;
  denominator_ *= other.numerator_.sign();
  divide_on_gsd();
  return *this;
}

Rational operator+(const Rational& first, const Rational& second) {
  Rational result = first;
  result += second;
  return result;
}

Rational operator-(const Rational& first, const Rational& second) {
  Rational result = first;
  result -= second;
  return result;
}

Rational operator*(const Rational& first, const Rational& second) {
  Rational result = first;
  result *= second;
  return result;
}

Rational operator/(const Rational& first, const Rational& second) {
  Rational result = first;
  result /= second;
  return result;
}

bool operator<(const Rational& first, const Rational& second) {
  return (first - second).numerator_ < 0;
}

bool operator>(const Rational& first, const Rational& second) {
  return second < first;
}

bool operator==(const Rational& first, const Rational& second) {
  return !(first < second) && !(second < first);
}

bool operator!=(const Rational& first, const Rational& second) {
  return !(first == second);
}

bool operator<=(const Rational& first, const Rational& second) {
  return !(first > second);
}

bool operator>=(const Rational& first, const Rational& second) {
  return !(first < second);
}

std::string Rational::toString() const {
  if (denominator_ == 1) {
    return numerator_.toString();
  }
  return (numerator_.toString() + '/' + denominator_.toString());
}

std::string Rational::asDecimal(size_t precision) const {
  std::string factor = "1";
  for (size_t i = 0; i < precision; ++i) {
    factor += "0";
  }
  BigInteger multiplier = factor;
  std::string str = abs((numerator_ * multiplier) / denominator_).toString();

  int cnt_pref_nulls = std::max(
      0, static_cast<int>(precision) - static_cast<int>(str.length()) + 1);
  std::string str1;
  for (int i = 0; i < cnt_pref_nulls; ++i) {
    str1 += "0";
  }
  str1 += str;

  std::string result;
  if (numerator_.sign() == -1) {
    result += '-';
  }
  result += str1.substr(0, str1.length() - precision);
  result += '.';
  result += str1.substr(str1.length() - precision, precision);
  return result;
}

Rational::operator double() const { return std::stod(asDecimal(40)); }

// a

bool constexpr is_prime(size_t N) {
  for (size_t i = 2; i <= std::sqrt(N); ++i) {
    if (N % i == 0) {
      return false;
    }
  }
  return true;
}

template <size_t N>
class Residue {
 private:
  size_t value_ = 0;

  static constexpr bool is_prime_ = is_prime(N);

  Residue<N> operator-() const;

  Residue<N> inverted() const;

 public:
  Residue() : value_(0){};

  Residue(int x);

  explicit operator int() const;

  Residue<N>& operator+=(Residue<N> other);

  Residue<N>& operator-=(Residue<N> other);

  Residue<N>& operator*=(Residue<N> other);

  Residue<N>& operator/=(Residue<N> other);
};

template <size_t N>
Residue<N>::Residue(int x) {
  long long value = static_cast<long long>(x) % static_cast<long long>(N);
  value_ = value > 0 ? value : value + N;
};

template <size_t N>
Residue<N>::operator int() const {
  return static_cast<int>(value_);
};

template <size_t N>
Residue<N> Residue<N>::operator-() const {
  Residue<N> result = *this;
  result.value_ = N - value_;
  return result;
}

template <size_t N>
Residue<N>& Residue<N>::operator+=(Residue<N> other) {
  value_ += other.value_;
  value_ %= N;
  return *this;
};

template <size_t N>
Residue<N>& Residue<N>::operator-=(Residue<N> other) {
  *this += -other;
  return *this;
}

template <size_t N>
Residue<N>& Residue<N>::operator*=(Residue<N> other) {
  value_ *= other.value_;
  value_ %= N;
  return *this;
};

template <size_t N>
Residue<N> Residue<N>::inverted() const {
  static_assert(Residue<N>::is_prime_);
  Residue<N> result(1);
  for (int i = 0; i < N - 2; ++i) {
    result *= *this;
  }

  return result;
}

template <size_t N>
Residue<N>& Residue<N>::operator/=(Residue<N> other) {
  *this *= other.inverted();
  return *this;
};

template <size_t N>
Residue<N> operator+(const Residue<N>& first, const Residue<N>& second) {
  Residue<N> result = first;
  result += second;
  return result;
}

template <size_t N>
Residue<N> operator-(const Residue<N>& first, const Residue<N>& second) {
  Residue<N> result = first;
  result -= second;
  return result;
}

template <size_t N>
Residue<N> operator*(const Residue<N>& first, const Residue<N>& second) {
  Residue<N> result = first;
  result *= second;
  return result;
}

template <size_t N>
Residue<N> operator/(const Residue<N>& first, const Residue<N>& second) {
  Residue<N> result = first;
  result /= second;
  return result;
}

template <size_t N, typename Field>
std::array<Field, N> multiplication(const std::array<Field, N>& first,
                                    const std::array<Field, N>& second) {
  std::array<Field, N> result;
  for (size_t i = 0; i < N; ++i) {
    result[i] = first[i] * second[i];
  }
  return result;
}

template <size_t M, size_t N, typename Field = Rational>
class Matrix {
 private:
  std::array<Field, M * N> arr;

  template <typename T>
  using pass_t = std::conditional_t<std::is_class_v<T>, const T&, T>;

  void swapRows(size_t first, size_t second);  // OK

  void increaseRow(size_t index, pass_t<Field> coefficient);

  void substractRow(size_t diminutive, size_t deductible,
                    pass_t<Field> coefficient);  // OK

  size_t toStep();  // OK

  void gaussMethod(Matrix<M, M, Field>& additional);  // OK

  void transpose();  // OK

  Field& data(size_t i, size_t j) { return arr[i * M + j]; }  // OK

  void reflect();

 public:
  Matrix();  // OK

  Matrix(std::initializer_list<std::array<Field, N>> list) : Matrix() {
    size_t i = 0;
    for (std::array<Field, N> row : list) {
      for (size_t j = 0; j < N; ++j) {
        arr[i + j] = row[j];
      }
      i += N;
    }
  }

  static Matrix<M, M, Field> unityMatrix();  // OK

  Matrix<M, N, Field>& operator+=(const Matrix<M, N, Field>& other);  // OK

  Matrix<M, N, Field>& operator-=(const Matrix<M, N, Field>& other);  // OK

  Matrix<M, N, Field>& operator*=(pass_t<Field> other);  // OK

  // template <size_t K>
  // Matrix<M, K, Field>& operator*=(const Matrix<N, K, Field>& other);

  Field det() const;  // OK

  Matrix<N, M, Field> transposed() const;  // OK

  size_t rank() const;  // OK

  Matrix<M, M, Field> inverted() const;  // OK

  void invert();  // OK

  Field trace() const;  // OK

  std::array<Field, N> getRow(size_t index) const;  // OK

  std::array<Field, M> getColumn(size_t index) const;  // OK

  Field& operator[](size_t i, size_t j);  // OK

  const Field& operator[](size_t i, size_t j) const;  // OK
};

template <size_t N, typename Field>
using SquareMatrix = Matrix<N, N, Field>;

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field>::Matrix() {
  arr.fill(0);
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field>& Matrix<M, N, Field>::operator+=(
    const Matrix<M, N, Field>& other) {
  for (size_t i = 0; i < arr.size(); ++i) {
    arr[i] += other.arr[i];
  }

  return *this;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field>& Matrix<M, N, Field>::operator-=(
    const Matrix<M, N, Field>& other) {
  for (size_t i = 0; i < arr.size(); ++i) {
    arr[i] -= other.arr[i];
  }

  return *this;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field>& Matrix<M, N, Field>::operator*=(pass_t<Field> other) {
  for (size_t i = 0; i < arr.size(); ++i) {
    arr[i] *= other;
  }

  return *this;
}

template <size_t M, size_t N, typename Field>
Field& Matrix<M, N, Field>::operator[](size_t i, size_t j) {
  return arr[i * M + j];
};

template <size_t M, size_t N, typename Field>
const Field& Matrix<M, N, Field>::operator[](size_t i, size_t j) const {
  return arr[i * M + j];
};

template <size_t M, size_t N, typename Field>
std::array<Field, N> Matrix<M, N, Field>::getRow(size_t index) const {
  std::array<Field, N> row;
  std::copy(arr[index * M], arr[index * (M + 1)], row);
  return row;
};

template <size_t M, size_t N, typename Field>
std::array<Field, M> Matrix<M, N, Field>::getColumn(size_t index) const {
  std::array<Field, M> column;
  for (size_t i = 0; i < M; ++i) {
    column[i] = arr[i, index];
  }
  return column;
};

template <size_t M, size_t N, size_t K, typename Field>
Matrix<M, K, Field> operator*(const Matrix<M, N, Field>& first,
                              const Matrix<N, K, Field>& second) {
  Matrix<M, K, Field> result;
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < K; ++j) {
      result[i, j] = multiplication(first.getRow(i), second.getColumn(j));
    }
  }
  return result;
}

template <size_t M, size_t N, typename Field>
void Matrix<M, N, Field>::swapRows(size_t first, size_t second) {
  size_t first_begin = first * M;
  size_t second_begin = second * M;
  for (size_t i = 0; i < N; ++i) {
    std::swap(arr[first_begin + i], arr[second_begin + i]);
  }
}

template <size_t M, size_t N, typename Field>
void Matrix<M, N, Field>::substractRow(size_t diminutive, size_t deductible,
                                       pass_t<Field> coefficient) {
  for (size_t i = 0; i < N; ++i) {
    data(diminutive, i) -= data(deductible, i) * coefficient;
  }
}

template <size_t M, size_t N, typename Field>
size_t Matrix<M, N, Field>::toStep() {
  size_t count_swaps = 0;
  size_t current_row = 0;
  for (size_t i = 0; i < N - 1 && current_row < M - 1; ++i) {
    bool is_null = true;
    if (data(current_row, i) == 0) {
      for (size_t j = current_row + 1; j < M; ++j) {
        if (data(j, i) != 0) {
          swapRows(current_row, j);
          ++count_swaps;
          is_null = false;
          break;
        }
      }
    } else {
      is_null = false;
    }

    if (is_null) {
      continue;
    }

    for (size_t j = current_row + 1; j < M; ++j) {
      if (data(j, i) != 0) {
        substractRow(j, current_row, data(j, i) / data(current_row, i));
      }
    }
    ++current_row;
  }

  return count_swaps;
}

template <size_t M, size_t N, typename Field>
Field Matrix<M, N, Field>::det() const {
  static_assert(M == N);
  Matrix<M, N, Field> copy = *this;
  size_t count_swaps = copy.toStep();
  Field result = 1;
  for (size_t i = 0; i < M; ++i) {
    result *= copy[i, i];
  }
  return (count_swaps % 2 == 0 ? result : -result);
}

template <size_t M, size_t N, typename Field>
size_t Matrix<M, N, Field>::rank() const {
  static_assert(M == N);
  Matrix<M, N, Field> copy = *this;
  copy.toStep();
  size_t rank = 0;
  size_t current_row = 0;
  size_t current_column = 0;

  while (current_column < N && current_row < M) {
    if (copy[current_row, current_column] == 0) {
      ++current_column;
      continue;
    } else {
      ++current_column;
      ++current_row;
      ++rank;
    }
  }
  return rank;
}

template <size_t M, size_t N, typename Field>
Matrix<N, M, Field> Matrix<M, N, Field>::transposed() const {
  Matrix<N, M, Field> result;
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      result[j, i] = data(i, j);
    }
  }
}

template <size_t M, size_t N, typename Field>
Field Matrix<M, N, Field>::trace() const {
  static_assert(M == N);
  Field trace = 0;
  for (size_t i = 0; i < M; ++i) {
    trace += data(i, i);
  }
  return trace;
}

template <size_t M, size_t N, typename Field>
Matrix<M, M, Field> Matrix<M, N, Field>::unityMatrix() {
  static_assert(M == N);
  Matrix<M, M, Field> unityMatrix;
  for (size_t i = 0; i < M; ++i) {
    unityMatrix[i, i] = 1;
  }
  return unityMatrix;
}

template <size_t M, size_t N, typename Field>
void Matrix<M, N, Field>::gaussMethod(Matrix<M, M, Field>& additional) {
  static_assert(M == N);
  for (size_t i = 0; i < N; ++i) {
    if (data(i, i) == 0) {
      for (size_t j = i + 1; j < M; ++j) {
        if (data(j, i) != 0) {
          swapRows(i, j);
          additional.swapRows(i, j);
        }
      }
    }

    for (size_t j = i + 1; j < M; ++j) {
      if (data(j, i) != 0) {
        Field coefficient = data(j, i) / data(i, i);
        substractRow(j, i, coefficient);
        additional.substractRow(j, i, coefficient);
      }
    }
  }
}

template <size_t M, size_t N, typename Field>
void Matrix<M, N, Field>::transpose() {
  static_assert(M == N);
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < i; ++j) {
      std::swap(data(i, j), data(j, i));
    }
  }
}

template <size_t M, size_t N, typename Field>
void Matrix<M, N, Field>::reflect() {
  static_assert(M == N);
  for (size_t i = 0; i < M / 2; ++i) {
    for (size_t j = 0; j < M; ++j) {
      std::swap(data(i, j), data(M - i - 1, M - j - 1));
    }
  }

  if (M % 2 == 1) {
    for (size_t i = 0; i < M / 2; ++i) {
      std::swap(data(M / 2, i), data(M / 2, M - i - 1));
    }
  }
}

template <size_t M, size_t N, typename Field>
void Matrix<M, N, Field>::increaseRow(size_t index, pass_t<Field> coefficient) {
  for (size_t i = 0; i < N; ++i) {
    data(index, i) *= coefficient;
  }
}

template <size_t M, size_t N, typename Field>
void Matrix<M, N, Field>::invert() {
  static_assert(M == N);
  Matrix<M, M, Field> result = Matrix<M, M, Field>::unityMatrix();

  gaussMethod(result);
  reflect();
  result.reflect();
  gaussMethod(result);
  for (size_t i = 0; i < M; ++i) {
    result.increaseRow(i, 1 / data(i, i));
  }

  result.reflect();
  std::swap(arr, result.arr);
}

template <size_t M, size_t N, typename Field>
Matrix<M, M, Field> Matrix<M, N, Field>::inverted() const {
  static_assert(M == N);
  Matrix<M, M, Field> result = *this;
  result.invert();
  return result;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator+(const Matrix<M, N, Field>& first,
                              const Matrix<M, N, Field>& second) {
  Matrix<M, N, Field> result = first;
  result += second;
  return result;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator-(const Matrix<M, N, Field>& first,
                              const Matrix<M, N, Field>& second) {
  Matrix<M, N, Field> result = first;
  result -= second;
  return result;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator*(const Matrix<M, N, Field>& first,
                              const Field& other) {
  Matrix<M, N, Field> result = first;
  result *= other;
  return result;
}
