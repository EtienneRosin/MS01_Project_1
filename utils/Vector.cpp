/*!
  \file Vector.cpp
  \authors E. Rosin
  \date 30 september 2024

  \brief Implementation of the Vector class
*/

#include "Vector.hpp"

// constructors ---------------------------------------------------------------
Vector::Vector() : std::vector<double>() {}

Vector::Vector(int length) : std::vector<double>(length) {}

Vector::Vector(int length, const double &value) : std::vector<double>(length, value) {}

Vector::Vector(const Vector& u) : std::vector<double>(u) {}


template<typename Iterator>
Vector::Vector(Iterator begin, Iterator end) : std::vector<double>(begin, end) {
    if (begin > end) {
        throw std::invalid_argument("Invalid iterator range");
    }
}

// Vector::Vector(const vector<double>& v) : std::vector<double>(v) {}
Vector Vector::subVector(size_t start, size_t end) const {
    if (start >= end || end > this->size()) {
        throw std::invalid_argument("Invalid range for subVector");
    }
    return Vector(this->begin() + start, this->begin() + end); // Utilise le constructeur par itérateur
}

// I/O ------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const Vector& vector) {
    os << "[";
    for (size_t i = 0; i < vector.size(); ++i) {
        os << vector[i];
        if (i < vector.size() - 1) {
            os << ", ";
        }
    }
    os << "]";
    return os;
}

// member operators -----------------------------------------------------------
Vector& Vector::operator=(const Vector &u) { // Assignment operator
    if (this != &u) {
        std::copy(u.begin(), u.end(), this->begin());
    }
    return *this;
}

// internal Algebra -----------------------------------------------------------
Vector& Vector::operator+=(const Vector &u) {
    if (this->size() != u.size()) {
        throw std::invalid_argument("Vector sizes do not match for addition");
    }
    for (size_t i = 0; i < this->size(); ++i) {
        (*this)[i] += u[i];
    }
    return *this;
}

Vector& Vector::operator-=(const Vector &u) {
    if (this->size() != u.size()) {
        throw std::invalid_argument("Vector sizes do not match for subtraction");
    }
    for (size_t i = 0; i < this->size(); ++i) {
        (*this)[i] -= u[i];
    }
    return *this;
}

Vector& Vector::operator*=(const double &a) {
    for (size_t i = 0; i < this->size(); ++i) {
        (*this)[i] *= a;
    }
    return *this;
}

Vector& Vector::operator/=(const double &a) {
    if (a == 0.0) {
        throw std::runtime_error("Division by zero is not allowed.");
    }
    for (size_t i = 0; i < this->size(); ++i) {
        (*this)[i] /= a;
    }
    return *this;
}

Vector Vector::operator+(const Vector &u) {
    Vector result = *this;
    result += u; // Use the += operator we already defined
    return result;
}

Vector Vector::operator-(const Vector &u) {
    Vector result = *this;
    result -= u; // Use the -= operator we already defined
    return result;
}

Vector Vector::operator*(const double &a) {
    Vector result = *this;
    result *= a; // Use the *= operator we already defined
    return result;
}

Vector Vector::operator/(const double &a) {
    if (a == 0.0) {
        throw std::runtime_error("Division by zero is not allowed.");
    }
    Vector result = *this;
    result /= a; // Use the /= operator we already defined
    return result;
}

double Vector::dot(const Vector &other) const {
    if (this->size() != other.size()) {
        throw std::invalid_argument("Vector sizes do not match for dot product");
    }
    
    double result = 0.0;
    for (size_t i = 0; i < this->size(); ++i) {
        result += (*this)[i] * other[i];
    }
    return result;
}

double Vector::NormSquared() const{
    return this->dot(*this);
}

double Vector::norm() const {
    // double result = this->dot(*this);
    return std::sqrt(this->NormSquared()); 
}