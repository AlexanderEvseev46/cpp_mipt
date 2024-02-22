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
