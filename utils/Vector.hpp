/*!
  \file Vector.hpp
  \authors E. Rosin
  \date 30 september 2024

  \brief Definition of the Vector class
*/
#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <vector>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <iterator>


class Vector : public std::vector<double>
{
public:

    // constructors -------------------------------------------------
    Vector();                                    //! default constructor
    Vector(int length);                          //! constructor by length
    Vector(int length, const double &value);     //! constructor by length for constant vector
    Vector(const Vector &u);                     //! constructor by copy

    template<typename Iterator>
    Vector(Iterator begin, Iterator end);        //! constructor by iterator

    // member operators ---------------------------------------------
    Vector &operator=(const Vector &u);
    Vector subVector(size_t start, size_t end) const;

    // internal Algebra ---------------------------------------------
    Vector &operator+=(const Vector &u);     //! Element-wise addition 'v = v + u'
    Vector &operator-=(const Vector &u);     //!  Element-wise subtraction 'v = v - u'
    Vector &operator*=(const double &a);     //! Element-wise multiplication 'v = a * v'
    Vector &operator/=(const double &a);     //! Element-wise division 'v = (1/a) * v'
    Vector operator+(const Vector &u);       //! Element-wise addition 'v + u'
    Vector operator-(const Vector &u);       //! Element-wise subtraction 'v - u'
    Vector operator*(const double &a);       //! Scalar multiplication 'v * a'
    Vector operator/(const double &a);       //! Scalar division 'v / a'
    double dot(const Vector &other) const;   //! dot product
    double NormSquared() const;              //! square euclidian norm
    double norm() const;                     //! Euclidean norm
};

// I/O ------------------------------------------------------------------------
std::ostream &operator<<(std::ostream &os, const Vector &vector); //! stream output operator

// External Algebra -----------------------------------------------------------

#endif  // VECTOR_HPP