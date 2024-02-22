#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

const double pi = 3.141592653589793;
const double epsilon = 0.0000000000001;

struct Vector;

class Line;

struct Point {
  double x;
  double y;

  Point() = default;

  Point(double x, double y) : x(x), y(y) {}

  Point shift(const Vector& vector) const;

  Point projection(const Line& line) const;

  void rotate(const Point& center, double angle);

  Point reflect(const Point& center) const;

  Point reflect(const Line& axis) const;

  Point scale(const Point& center, double coefficient) const;
};

std::ostream& operator<<(std::ostream& out, Point p) {
  out << p.x << " " << p.y;
  return out;
}

bool operator==(const Point& first, const Point& second) {
  return (std::abs(first.x - second.x) < epsilon &&
          std::abs(first.y - second.y) < epsilon);
}

bool operator!=(const Point& first, const Point& second) {
  return !(first == second);
}

double distanse(const Point& first, const Point& second) {
  return std::sqrt((first.x - second.x) * (first.x - second.x) +
                   (first.y - second.y) * (first.y - second.y));
}

Point middle(const Point& first, const Point& second) {
  Point middle((first.x + second.x) / 2, (first.y + second.y) / 2);
  return middle;
}

class Line {
 public:
  double angular_coefficient;
  double y_shift;

  Line(const Point& a, const Point& b) {
    if (std::abs(a.x - b.x) < epsilon) {
      angular_coefficient = INFINITY;
      y_shift = a.x;
      return;
    }
    angular_coefficient = (b.y - a.y) / (b.x - a.x);
    y_shift = a.y - angular_coefficient * a.x;
  }

  Line(const double& angular_coefficient, const double& y_shift)
      : angular_coefficient(angular_coefficient), y_shift(y_shift) {}

  Line(const Point& point, const double& angular_coefficient)
      : angular_coefficient(angular_coefficient) {
    if (angular_coefficient == INFINITY) {
      y_shift = point.x;
      return;
    }
    y_shift = point.y - angular_coefficient * point.x;
  }

  Line rotate(const Point& point, double angle) const {
    if (std::abs(std::tan(std::atan(angular_coefficient) + angle)) >
        100000000000) {
      return Line(point, INFINITY);
    }

    double coefficient = std::tan(std::atan(angular_coefficient) + angle);
    Line result(point, coefficient);
    return result;
  }
};

std::ostream& operator<<(std::ostream& out, Line l) {
  out << l.angular_coefficient << " " << l.y_shift;
  return out;
}

Point splitInRatio(const Point& first, const Point& second,
                   double coefficient_1, double coefficient_2) {
  double x = (first.x * coefficient_2 + second.x * coefficient_1) /
             (coefficient_1 + coefficient_2);

  double y = (first.y * coefficient_2 + second.y * coefficient_1) /
             (coefficient_1 + coefficient_2);

  return Point(x, y);
}

struct Vector {
  double x;
  double y;

  double length() const { return std::sqrt(x * x + y * y); }

  double angular_coefficient() const { return y / x; }

  Vector() : x(0), y(0) {}

  Vector(const Point& first, const Point& second)
      : x(second.x - first.x), y(second.y - first.y) {}

  Vector(double angular_coefficient, double length) {
    if (angular_coefficient == INFINITY) {
      x = 0;
      y = length;
      return;
    }
    x = length / sqrt(1 + angular_coefficient * angular_coefficient);
    y = sqrt(length * length - x * x);
  }

  Vector& operator*=(double a) {
    x *= a;
    y *= a;
    return *this;
  }

  Vector& operator/=(double a) {
    *this *= (1 / a);
    return *this;
  }

  Vector& operator+=(const Vector& other) {
    x += other.x;
    y += other.y;
    return *this;
  }

  void rotate(double angle) {
    double x1 = x * cos(angle) - y * sin(angle);
    y = x * sin(angle) + y * cos(angle);
    x = x1;
  }
};

Vector operator+(const Vector& first, const Vector& second) {
  Vector result = first;
  result += second;
  return result;
}

Vector operator*(Vector vector, double a) {
  vector *= a;
  return vector;
}

Vector operator/(Vector vector, double a) {
  vector /= a;
  return vector;
}

Point Point::shift(const Vector& vector) const {
  Point result(x + vector.x, y + vector.y);
  return result;
}

bool operator==(const Line& first, const Line& second) {
  return (std::abs(first.angular_coefficient - second.angular_coefficient) <
              epsilon &&
          std::abs(first.y_shift - second.y_shift) < epsilon);
}

bool operator!=(const Line& first, const Line& second) {
  return !(first == second);
}

Point cross(const Line& first, const Line& second) {
  if (first.angular_coefficient == INFINITY) {
    double y = second.y_shift + second.angular_coefficient * first.y_shift;
    return Point(first.y_shift, y);
  }

  if (second.angular_coefficient == INFINITY) {
    double y = first.y_shift + first.angular_coefficient * second.y_shift;
    return Point(second.y_shift, y);
  }

  double x = (first.y_shift - second.y_shift) /
             (second.angular_coefficient - first.angular_coefficient);
  double y = first.y_shift + first.angular_coefficient * x;
  return Point(x, y);
}

// double orientedAngle(const Line& first, const Line& second) {
//   Point crossing = cross(first, second);
//   double angle =
//       atan(second.angular_coefficient) - atan(first.angular_coefficient);
//   if (angle < -pi / 2) {
//     angle += pi;
//   }
//   if (angle > pi / 2) {
//     angle -= pi;
//   }
//   if (std::abs(std::abs(angle) - pi / 2) < epsilon) {
//     return pi / 2;
//   }
//   Line line = first.rotate(crossing, angle);
//   return line == second ? angle : -angle;
// }

Point Point::projection(const Line& line) const {
  if (line.angular_coefficient == INFINITY) {
    return Point(line.y_shift, y);
  }
  Line parallel(*this, line.angular_coefficient);
  Line perpendicular = parallel.rotate(*this, pi / 4);
  return cross(line, perpendicular);
}

double scalarProduct(const Vector& first, const Vector& second) {
  return (first.x * second.x + first.y * second.y);
}

double vectorProduct(const Vector& first, const Vector& second) {
  return first.x * second.y - first.y * second.x;
}

double angle(const Vector& first, const Vector& second) {
  double angle =
      acos(scalarProduct(first, second) / (first.length() * second.length()));

  return vectorProduct(first, second) >= 0 ? angle : -angle;
}

void Point::rotate(const Point& center, double angle) {
  Vector vector(center, *this);
  vector.rotate(angle);
  *this = center.shift(vector);
}

Point Point::reflect(const Point& center) const {
  Vector vector(*this, center);
  vector *= 2;
  return shift(vector);
}

Point Point::reflect(const Line& axis) const {
  if (axis.angular_coefficient == INFINITY) {
    return Point(2 * axis.y_shift - x, y);
  }

  if (std::abs(axis.angular_coefficient) < epsilon) {
    return Point(x, 2 * axis.y_shift - y);
  }
  double coefficient = -1 / axis.angular_coefficient;

  Point crossing = cross(Line(*this, coefficient), axis);
  return reflect(crossing);
}

Point Point::scale(const Point& center, double coefficient) const {
  Vector vector(center, *this);
  vector *= coefficient;
  return center.shift(vector);
}

class Shape {
 public:
  virtual double perimeter() const = 0;

  virtual double area() const = 0;

  virtual bool operator==(const Shape& other) const = 0;

  virtual bool operator!=(const Shape& other) const = 0;

  virtual bool isCongruentTo(const Shape& another) const = 0;

  virtual bool isSimilarTo(const Shape& another) const = 0;

  virtual bool containsPoint(const Point& point) const = 0;

  virtual void rotate(const Point& center, double angle) = 0;

  virtual void reflect(const Point& center) = 0;

  virtual void reflect(const Line& axis) = 0;

  virtual void scale(const Point& center, double coefficient) = 0;

  virtual ~Shape() = default;
};

class Polygon : public Shape {
 protected:
  std::vector<Point> points_;

 public:
  ~Polygon() override = default;

  Polygon(std::vector<Point> points) { points_ = points; }

  void fillVector() {}

  template <class... Tail>
  void fillVector(Point head, Tail... tail) {
    points_.push_back(head);
    fillVector(tail...);
  }

  template <class... Tail>
  Polygon(Point head, Tail... tail) {
    fillVector(head, tail...);
  }

  int verticesCount() const { return points_.size(); }

  std::vector<Point> getVertices() const { return points_; }

  bool isConvex() const {
    Vector side_2(points_[0], points_[1]);
    Vector side_1(points_[1], points_[2]);
    double area = vectorProduct(side_2, side_1);
    for (size_t i = 1; i < points_.size(); ++i) {
      side_2 = Vector(points_[i + 1], points_[(i + 2) % points_.size()]);
      double new_area = vectorProduct(side_1, side_2);
      side_1 = side_2;
      if ((area >= 0 && new_area < 0) || (area <= 0 && new_area > 0)) {
        return false;
      }
    }

    return true;
  }

  double perimeter() const final {
    double perimeter = 0;
    for (size_t i = 0; i < points_.size() - 1; ++i) {
      perimeter += distanse(points_[i], points_[i + 1]);
    }
    perimeter += distanse(points_[0], points_[points_.size() - 1]);
    return perimeter;
  }

  double area() const final {
    double area = 0;
    Vector segment(points_[0], points_[1]);
    for (size_t i = 0; i < points_.size() - 2; ++i) {
      Vector segment_2(points_[0], points_[i + 2]);
      area += vectorProduct(segment, segment_2);
      segment = segment_2;
    }
    return std::abs(area) / 2;
  }

  bool isShift(const Polygon* polygon) const {
    for (size_t i = 0; i < points_.size(); ++i) {
      bool is_same = true;
      for (size_t j = 0; j < points_.size(); ++j) {
        if (polygon->points_[j] != points_[(i + j) % points_.size()]) {
          is_same = false;
          break;
        }
      }
      if (is_same) {
        return true;
      }
    }

    return false;
  }

  bool operator==(const Shape& other) const final {
    const Polygon* polygon = dynamic_cast<const Polygon*>(&other);
    if (polygon == nullptr) {
      return false;
    }

    if (polygon->points_.size() != points_.size()) {
      return false;
    }

    std::vector<Point> reversed_points = polygon->points_;
    std::reverse(reversed_points.begin(), reversed_points.end());
    Polygon reversed_polygon(reversed_points);

    return isShift(polygon) || isShift(&reversed_polygon);
  };

  bool operator!=(const Shape& other) const final { return !(*this == other); }

  std::vector<double> sides() const {
    std::vector<double> sides(points_.size());
    for (size_t i = 0; i < points_.size(); ++i) {
      sides[i] = distanse(points_[i], points_[(i + 1) % points_.size()]);
    }
    return sides;
  }

  std::vector<double> angles() const {
    std::vector<double> angles(points_.size());
    for (size_t i = 0; i < points_.size(); ++i) {
      Vector first_side(points_[i], points_[(i + 1) % points_.size()]);
      Vector second_side(points_[(i + 1) % points_.size()],
                         points_[(i + 2) % points_.size()]);
      angles[i] = angle(first_side, second_side);
    }
    return angles;
  }

  bool isSimilarShift(const std::vector<double>& sides_1,
                      const std::vector<double>& sides_2,
                      const std::vector<double>& angles_1,
                      const std::vector<double>& angles_2,
                      double coefficient) const {
    for (size_t i = 0; i < points_.size(); ++i) {
      bool is_same = true;
      int angle_direct = (angles_2[0] * angles_1[i] < 0 ? -1 : 1);
      for (size_t j = 0; j < points_.size(); ++j) {
        if (std::abs(sides_1[(i + j) % points_.size()] * coefficient -
                     sides_2[j]) > epsilon ||
            std::abs(angles_1[(i + j) % points_.size()] * angle_direct -
                     angles_2[j]) > epsilon) {
          is_same = false;
          break;
        }
      }
      if (is_same) {
        return true;
      }
    }
    return false;
  }

  bool isSimilarTo(const Shape& another) const final {
    const Polygon* polygon = dynamic_cast<const Polygon*>(&another);
    if (polygon == nullptr) {
      return false;
    }

    if (polygon->points_.size() != points_.size()) {
      return false;
    }

    std::vector<double> sides_1 = sides();
    std::vector<double> angles_1 = angles();
    std::vector<double> sides_2 = polygon->sides();
    std::vector<double> angles_2 = polygon->angles();

    double coefficient = polygon->perimeter() / perimeter();

    if (isSimilarShift(sides_1, sides_2, angles_1, angles_2, coefficient)) {
      return true;
    }

    std::vector<Point> reversed_points = polygon->points_;
    std::reverse(reversed_points.begin(), reversed_points.end());
    Polygon reversed_polygon(reversed_points);
    std::vector<double> sides_3 = reversed_polygon.sides();
    std::vector<double> angles_3 = reversed_polygon.angles();

    return isSimilarShift(sides_1, sides_3, angles_1, angles_3, coefficient);
  };

  bool isCongruentTo(const Shape& another) const final {
    return isSimilarTo(another) &&
           std::abs(another.perimeter() - perimeter()) < epsilon;
  }

  bool containsPoint(const Point& point) const final {
    // for (size_t i = 0; i < points_.size(); ++i) {
    //   std::cout << points_[i] << std::endl;
    // }
    // std::cout << point << "aaa" << std::endl;

    std::vector<std::pair<Point, Point>> sides(points_.size());
    for (size_t i = 0; i < points_.size() - 1; ++i) {
      sides[i] = {points_[i], points_[i + 1]};
    }
    sides[points_.size() - 1] = {points_[points_.size() - 1], points_[0]};

    Point start = splitInRatio(points_[0], points_[1], 119, 244);

    Line line(point, start);
    Vector vector(point, start);

    int count_crosses = 0;
    for (size_t i = 0; i < points_.size(); ++i) {
      Line side(sides[i].first, sides[i].second);
      if (line.angular_coefficient == side.angular_coefficient) {
        continue;
      }
      Point crossing = cross(line, side);
      if (crossing.x >= std::min(sides[i].first.x, sides[i].second.x) &&
          crossing.x <= std::max(sides[i].first.x, sides[i].second.x) &&
          crossing.y >= std::min(sides[i].first.y, sides[i].second.y) &&
          crossing.y <= std::max(sides[i].first.y, sides[i].second.y)) {
        Vector cross_vector(point, crossing);
        if (scalarProduct(vector, cross_vector) >= 0) {
          ++count_crosses;
        }
      }
    }

    return ((count_crosses % 2 == 0) ? false : true);
  };

  void rotate(const Point& center, double angle) final {
    for (size_t i = 0; i < points_.size(); ++i) {
      points_[i].rotate(center, angle);
    }
  };

  void reflect(const Point& center) final {
    for (size_t i = 0; i < points_.size(); ++i) {
      points_[i] = points_[i].reflect(center);
    }
  };

  void reflect(const Line& axis) final {
    for (size_t i = 0; i < points_.size(); ++i) {
      points_[i] = points_[i].reflect(axis);
    }
  };

  void scale(const Point& center, double coefficient) final {
    for (size_t i = 0; i < points_.size(); ++i) {
      points_[i] = points_[i].scale(center, coefficient);
    }
  };
};
class Ellipse : public Shape {
 protected:
  std::pair<Point, Point> focuses_;
  double sum_of_length_;

 public:
  ~Ellipse() override = default;

  Ellipse(const Point& first, const Point& second, double sum_of_length)
      : focuses_({first, second}), sum_of_length_(sum_of_length) {}

  std::pair<Point, Point> focuses() const { return focuses_; }

  Point center() const { return middle(focuses_.first, focuses_.second); }

  double eccentricity() const {
    return distanse(focuses_.first, focuses_.second) / sum_of_length_;
  }

  std::pair<Line, Line> directrices() {
    Vector vector(focuses_.first, focuses_.second);
    vector /= vector.length();
    vector *= sum_of_length_ / eccentricity() / 2;
    Point point_on_first_directrice = center().shift(vector);
    Point point_on_second_directrice = center().shift(vector * (-1));
    Line axis(focuses_.first, focuses_.second);

    std::pair<Line, Line> directrices = {
        axis.rotate(point_on_first_directrice, pi / 2),
        axis.rotate(point_on_second_directrice, pi / 2)};
    return directrices;
  }

  double small_axis() const {
    return sqrt(sum_of_length_ * sum_of_length_ / 4 -
                distanse(focuses_.first, focuses_.second) *
                    distanse(focuses_.first, focuses_.second) / 4);
  }

  double big_axis() const { return sum_of_length_ / 2; }

  double perimeter() const final {
    return 2 * sum_of_length_ * std::comp_ellint_2(eccentricity());
  }

  double area() const final { return pi * big_axis() * small_axis(); }

  bool operator==(const Shape& other) const final {
    if (const Ellipse* another = dynamic_cast<const Ellipse*>(&other);
        another) {
      return (focuses_ == another->focuses_ ||
              (focuses_.second == another->focuses_.first &&
               focuses_.first == another->focuses_.second)) &&
             std::abs(sum_of_length_ - another->sum_of_length_) < epsilon;
    }

    return false;
  }

  bool operator!=(const Shape& other) const final { return !(*this == other); }

  bool isCongruentTo(const Shape& another) const final {
    if (const Ellipse* other = dynamic_cast<const Ellipse*>(&another); other) {
      return std::abs(small_axis() - other->small_axis()) < epsilon &&
             std::abs(big_axis() - other->big_axis()) < epsilon;
    }

    return false;
  }

  bool isSimilarTo(const Shape& another) const final {
    if (const Ellipse* other = dynamic_cast<const Ellipse*>(&another); other) {
      return std::abs(small_axis() / other->small_axis() -
                      big_axis() / other->big_axis()) < epsilon;
    }

    return false;
  }

  bool containsPoint(const Point& point) const final {
    return distanse(point, focuses_.first) + distanse(point, focuses_.second) <=
           sum_of_length_ + epsilon;
  }

  void rotate(const Point& center, double angle) final {
    focuses_.first.rotate(center, angle);
    focuses_.second.rotate(center, angle);
  }

  void reflect(const Point& center) final {
    focuses_.first = focuses_.first.reflect(center);
    focuses_.second = focuses_.second.reflect(center);
  }

  void reflect(const Line& axis) final {
    focuses_.first = focuses_.first.reflect(axis);
    focuses_.second = focuses_.second.reflect(axis);
  }

  void scale(const Point& center, double coefficient) final {
    focuses_.first = focuses_.first.scale(center, coefficient);
    focuses_.second = focuses_.second.scale(center, coefficient);
    sum_of_length_ *= std::abs(coefficient);
  }
};

class Circle : public Ellipse {
 public:
  ~Circle() override = default;

  double radius() const { return sum_of_length_ / 2; }

  Circle(const Point& center, double radius)
      : Ellipse(center, center, radius * 2) {}
};

class Rectangle : public Polygon {
 public:
  ~Rectangle() override = default;

  Point center() const { return middle(points_[0], points_[2]); }

  std::pair<Line, Line> diagonals() const {
    return {Line(points_[0], points_[2]), Line(points_[1], points_[3])};
  }

  std::vector<Point> RectanglePoints(const Point& first, const Point& second,
                                     double ratio) {
    double min_side = distanse(first, second) / std::sqrt(1 + ratio * ratio);
    double max_side = ratio * min_side;

    if (ratio < 1) {
      std::swap(min_side, max_side);
    }

    Vector min(first, second);
    Vector max(first, second);
    min.rotate(atan(max_side / min_side));
    max.rotate(atan(max_side / min_side) - pi / 2);
    min *= (min_side / distanse(first, second));
    max *= (max_side / distanse(first, second));

    Point third = first.shift(min);
    Point fourth = first.shift(max);

    return {first, third, second, fourth};
  }

  Rectangle(const Point& first, const Point& second, double ratio)
      : Polygon(RectanglePoints(first, second, ratio)) {}
};

class Square : public Rectangle {
 public:
  ~Square() override = default;

  Square(const Point& first, const Point& second)
      : Rectangle(first, second, 1) {}

  Circle circumscribedCircle() const {
    double radius = distanse(points_[0], points_[1]) / 2;
    return Circle(center(), radius);
  }

  Circle inscribedCircle() const {
    double radius = distanse(points_[0], points_[2]) / 2;
    return Circle(center(), radius);
  }
};

class Triangle : public Polygon {
 public:
  Triangle(const Point& a, const Point& b, const Point& c) : Polygon(a, b, c) {}

  Circle circumscribedCircle() const {
    Point middle_1 = middle(points_[0], points_[1]);
    Point middle_2 = middle(points_[0], points_[2]);

    Line side_1(points_[0], points_[1]);
    Line side_2(points_[0], points_[2]);

    Line middle_perpendicular_1 = side_1.rotate(middle_1, pi / 2);
    Line middle_perpendicular_2 = side_2.rotate(middle_2, pi / 2);

    Point center = cross(middle_perpendicular_1, middle_perpendicular_2);
    double radius = distanse(center, points_[0]);
    return Circle(center, radius);
  }

  Circle inscribedCircle() const {
    double a = distanse(points_[0], points_[1]);
    double b = distanse(points_[1], points_[2]);
    double c = distanse(points_[2], points_[0]);
    double p = perimeter() / 2;
    Point first = splitInRatio(points_[0], points_[1], p - b, p - c);
    Point second = splitInRatio(points_[1], points_[2], p - c, p - a);
    Point third = splitInRatio(points_[2], points_[0], p - a, p - b);
    Triangle triangle(first, second, third);
    return triangle.circumscribedCircle();
  }

  Point centroid() const {
    double x = (points_[0].x + points_[1].x + points_[2].x) / 3;
    double y = (points_[0].y + points_[1].y + points_[2].y) / 3;
    return Point(x, y);
  }

  Point orthocenter() const {
    Point o = circumscribedCircle().center();
    return centroid().scale(o, 3);
  }

  Line EulerLine() const { return Line(centroid(), orthocenter()); }

  Circle ninePointsCircle() const {
    Point first = middle(points_[0], points_[1]);
    Point second = middle(points_[1], points_[2]);
    Point third = middle(points_[2], points_[0]);
    Triangle triangle(first, second, third);
    // std::cout << points_[0] << std::endl
    //           << points_[1] << std::endl
    //           << points_[2] << std::endl
    //           << centroid() << std::endl
    //           << circumscribedCircle().center() << std::endl
    //           << orthocenter() << std::endl
    //           << triangle.circumscribedCircle().center() << std::endl;
    return triangle.circumscribedCircle();
  }
};