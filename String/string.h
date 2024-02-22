#include <algorithm>
#include <cstring>
#include <iostream>

class String {
 private:
  size_t sz = 0;

  size_t cap = 1;

  char* arr = nullptr;

  String(size_t count);

  bool is_substr_equal(const String& other, size_t index) const;

 public:
  String(char c);

  String();

  String(const char* str);

  String(size_t count, char c);

  String(const String& other);

  void swap(String& other);

  String& operator=(String other);

  const char& operator[](size_t index) const;

  char& operator[](size_t index);

  size_t length() const;

  size_t size() const;

  size_t capacity() const;

  const char* data() const;

  char* data();

  void push_back(char c);

  void pop_back();

  const char& front() const;

  char& front();

  const char& back() const;

  char& back();

  String& operator+=(char c);

  String& operator+=(const String& str);

  String substr(size_t start, size_t count) const;

  size_t find(const String& substr) const;

  size_t rfind(const String& substr) const;

  bool empty() const;

  void clear();

  void shrink_to_fit();

  ~String();
};

String::String(size_t count)
    : sz(count), cap(count + 1), arr(new char[count + 1]) {
  arr[sz] = '\0';
}

String::String(char c) : String(static_cast<size_t>(1)) { arr[0] = c; }

String::String() : String(static_cast<size_t>(0)) {}

String::String(const char* str) : String(strlen(str)) { memcpy(arr, str, sz); }

String::String(size_t count, char c) : String(count) { memset(arr, c, count); }

String::String(const String& other) : String(other.sz) {
  memcpy(arr, other.arr, sz);
}

void String::swap(String& other) {
  std::swap(sz, other.sz);
  std::swap(cap, other.cap);
  std::swap(arr, other.arr);
}

String& String::operator=(String other) {
  swap(other);
  return *this;
}

const char& String::operator[](size_t index) const { return arr[index]; }

char& String::operator[](size_t index) { return arr[index]; }

size_t String::length() const { return sz; }

size_t String::size() const { return sz; }

size_t String::capacity() const { return cap - 1; }

const char* String::data() const { return arr; }

char* String::data() { return arr; }

void String::push_back(char c) {
  if (sz + 1 == cap) {
    char* str = new char[cap * 2];
    memcpy(str, arr, sz);
    delete[] arr;
    arr = str;
    cap *= 2;
  }
  arr[sz++] = c;
  arr[sz] = '\0';
}

void String::pop_back() { arr[--sz] = '\0'; }

const char& String::front() const { return arr[0]; }

char& String::front() { return arr[0]; }

const char& String::back() const { return arr[sz - 1]; }

char& String::back() { return arr[sz - 1]; }

String& String::operator+=(char c) {
  push_back(c);
  return *this;
}

String& String::operator+=(const String& str) {
  if (sz + str.length() + 1 > cap) {
    char* new_arr = new char[(sz + str.length() + 1) * 2];
    if (sz > 0) {
      std::copy(&arr[0], &arr[sz], new_arr);
    }
    delete[] arr;
    arr = new_arr;
    cap = (sz + str.length() + 1) * 2;
  }
  std::copy(&str.data()[0], &str.data()[str.length() + 1], &arr[sz]);
  sz += str.length();
  return *this;
}

String String::substr(size_t start, size_t count) const {
  if (start + count > sz) {
    count = sz - start;
  }

  String ans = String(count);
  std::copy(&arr[start], &arr[start + count], ans.data());
  return ans;
}

bool String::empty() const { return sz == 0; }

void String::clear() {
  sz = 0;
  arr[0] = '\0';
}

void String::shrink_to_fit() {
  String str = *this;
  swap(str);
}
String::~String() { delete[] arr; }

String operator+(const String& str1, const String& str2) {
  String ans = str1;
  ans += str2;
  return ans;
}

bool operator<(const String& str1, const String& str2) {
  if (strcmp(str1.data(), str2.data()) < 0) {
    return true;
  }
  return false;
}

bool operator>(const String& str1, const String& str2) { return str2 < str1; }

bool operator!=(const String& str1, const String& str2) {
  return str2 < str1 || str1 < str2;
}

bool operator==(const String& str1, const String& str2) {
  return !(str1 != str2);
}

bool operator<=(const String& str1, const String& str2) {
  return !(str2 < str1);
}

bool operator>=(const String& str1, const String& str2) {
  return !(str2 > str1);
}

bool String::is_substr_equal(const String& other, size_t index) const {
  for (size_t j = 0; j < other.length(); ++j) {
    if (other[j] != this->data()[index + j]) {
      return false;
    }
  }
  return true;
}

size_t String::find(const String& substr) const {
  if (substr.length() > sz) {
    return sz;
  }

  for (size_t i = 0; i <= sz - substr.length(); ++i) {
    if (is_substr_equal(substr, i)) {
      return i;
    }
  }
  return sz;
}

size_t String::rfind(const String& substr) const {
  if (substr.length() > sz) {
    return sz;
  }

  size_t result = sz;
  for (size_t i = 0; i <= sz - substr.length(); ++i) {
    if (is_substr_equal(substr, i)) {
      result = i;
    }
  }
  return result;
}

std::ostream& operator<<(std::ostream& out, const String& str) {
  out << str.data();
  return out;
}

std::istream& operator>>(std::istream& in, String& str) {
  char c = '1';
  while (c != '\0' && c != '\n' && c != ' ' && c != EOF) {
    c = in.peek();
    if (c != '\0' && c != '\n' && c != ' ' && c != EOF) {
      str.push_back(in.get());
    }
  }
  c = in.get();
  return in;
}