/***********************************************
 *              DESCRIPTION                    *
 * File contains a definition of a class
 * TriangularMatrix in the namespace sd.
 *
 * The class represents a two-dimensional
 * triangular matrix (upper or lower), which is
 * a special case of a square matrix. It inherits
 * from SquareMatrix and enforces the triangular
 * structure by restricting write access to the
 * elements outside the allowed triangular part.
 *
 * Two types of triangular matrices are supported:
 *  - Lower triangular matrix
 *  - Upper triangular matrix
 *
 * The triangular type is specified by the
 * TriangleType enum.
 *
 * Attempts to modify elements outside the valid
 * triangular region result in a runtime error.
 *
 * The class preserves all read operations defined
 * for square matrices while providing additional
 * safety for write operations, but using 
 * operator () for read an element and send it 
 * to output stream  can lead to unwanted execptions
 * and it is better to use getValue() function.
 ***********************************************/


#pragma once

#include "square.hpp"

namespace sd{

  enum class TriangleType { Lower, Upper }; // Matrix Types  

  template <typename T>
  requires std::is_arithmetic_v<T>
  class TriangularMatrix : public SquareMatrix<T> {
  private:

    using TVector = std::vector<T>;         
    using TMatrix = std::vector<TVector>;
    using SquareMatrix<T>::m_dim;   

    TriangleType m_type;

    bool writePermission (size_t row, size_t column) const {
      return (( m_type == TriangleType::Lower ) && ( row >= column ))
          || (( m_type == TriangleType::Upper ) && ( row <= column ));
    }

  public:
      /***********************************************
    *              CONSTRUCTORS                   *
    * 
    * No 1: Default constructor 
    * Create a 1x1 matrix with T-type 0 value.
    * 
    * Ex.:
    *    TriangularMatrix<double> matA;
    * 
    * Result:
    *     
    *   matA = [ 0.0 ]
    *          
    * No 2: Intializing a matrix with a given dimension, 
    * type and value.
    * .
    * Ex.:
    *    TriangularMatrix matA(3, sd::TriangleType::Upper, 1);
    * or better (with explicit type specification)
         TriangularMatrix<int> matA(3, sd::TriangleType::Upper, 1);

    * Result:
    *          [ 1 1 1 ]
    *   matA = | 0 1 1 |
    *          [ 0 0 1 ]
    * 
    * Note: the dimensions argument are of type int.
    * Providing a negative value would result in an
    * implicit conversion to a very large size_t
    * value, potentially causing excessive memory
    * allocation and program termination.
    * 
    * No 3: Copy constructor is the default one - it 
    * need no implementation.
    ***********************************************/  
    TriangularMatrix( int dim = 1, TriangleType type = TriangleType::Lower, T value = T{0}) : SquareMatrix<T>( dim, T{0} ), m_type( type ) {
      
      for ( size_t row = 0; row < m_dim; row++ ){
        for ( size_t column = 0; column < m_dim; column ++ ){
          if ( writePermission( row, column ) ){
            SquareMatrix<T>::operator()( row, column ) = value;
          }
        }
      }
    }

    TriangularMatrix( const SquareMatrix<T> & sqMat) : SquareMatrix<T>( sqMat ) {
      
      bool is_upper{true}, is_lower{true};

      //check if matrix is lower triangular matrix? -> above diagonal
      for( size_t row = 0; (row < m_dim && is_lower) ; row++ ) {
        for ( size_t col = row + 1; col < m_dim; col++ ) {
            if ( sqMat( row, col ) != T{} ) {
                is_lower = false;   //it will exit the outer for-loop as well 
                break;
            }
        }
      }
      //check if matrix is upper triangular matrix? -> below diagonal
      for ( size_t row = 1; (row < m_dim && is_upper); row++ ){
        for ( size_t col = 0; col < row; col++ ) {
            if ( sqMat( row, col ) != T{} ) {
                is_upper = false;      //it will exit the outer for-loop as well 
                break;
            }
        }
      }
      

      if ( is_upper ){
        m_type = TriangleType::Upper;
      } else if ( is_lower ){
        m_type = TriangleType::Lower;
      } else{
        throw std::invalid_argument( "Given matrix is not triangular!" );
      }
    }
    
    bool is_upper() const {
      return ( m_type == TriangleType::Upper );
    }

    bool is_lower() const {
      return ( m_type == TriangleType::Lower );
    }

    /***********************************************
    *          SETTERS and operator()              *
    * Write access to matrix elements is restricted
    * to the valid triangular region. Attempts to
    * modify elements outside the allowed region
    * result in a runtime exception.
    *
    * The non-const version of operator() is intended
    * for write access and therefore enforces the
    * triangular structure. Read-only access can be
    * safely performed via the const operator() or
    * getter functions.
    ***********************************************/

    T& operator () (size_t r, size_t c) { 
      
       if ( !writePermission( r, c ) ){
         throw std::domain_error( "Access denied to that index of triangular matrix!");
       }
      
      return SquareMatrix<T>::operator()( r, c );
      
    }

    T operator()(size_t r, size_t c) const {
    
      return SquareMatrix<T>::operator()(r, c);

    }

    void setValue( size_t row, size_t column, T value ) {

      //Validation
      if ( (row >= m_dim) || (column >= m_dim) ){
        throw std::out_of_range( "Out of range!" );
      }

      (*this)( row, column ) = value;  //operator() checks write permission

    }

    void setRow( size_t row, const TVector & vector ){

      if ( row >= m_dim ){
        throw std::out_of_range( "Out of range!" );  
      }
      if ( vector.size() != m_dim ){
        throw std::invalid_argument( "Rows dimensions mismatch!" );
      }

      for ( size_t column = 0; column < m_dim; column++ ){
        
        if ( !writePermission(row, column) ){
          
          if ( vector[ column ] != T{0} ){   //trying to write in forbiden place
            throw std::domain_error( "Access denied to that index of triangular matrix!");
          }
          continue;   //0 in proper place
        }
        (*this)( row, column ) = vector[ column ]; //writting in proper place
      }   
    }


    void setColumn( size_t column, const TVector & vector ){

      if ( column >= m_dim ){
        throw std::out_of_range( "Out of range!" );  
      }
      if ( vector.size() != m_dim ){
        throw std::invalid_argument( "Columns dimensions mismatch!" );
      }

      for ( size_t row = 0; row < m_dim; row++ ){
        if ( !writePermission(row, column) ){
          
          if ( vector.at( row ) != 0){   //trying to write in forbiden place
            throw std::domain_error( "Access denied to that index of triangular matrix!");
          }
          continue;   //0 in proper place
        }
        (*this)( row, column ) = vector.at( row ); //writting in proper place
      }
    }

    /***********************************************
    *              Transposition                   *
    * Transposition of a triangular matrix converts
    * an upper triangular matrix into a lower one
    * and vice versa.
    ***********************************************/

    TriangularMatrix transpose() const {
      
     return TriangularMatrix<T> (  SquareMatrix<T>::transpose() );

    }

    void transpose_in_place() {

      for ( size_t row = 0; row < m_dim; row++ ) {
        for ( size_t column = row + 1; column < m_dim; column++) {
            std::swap( SquareMatrix<T>::operator()( row, column ),  // matrix is converted, so there is need of using SquareMatrix operator() 
                       SquareMatrix<T>::operator()( column, row )); // it allows to change normany unwritable elements.
        }
      }

      if ( m_type == TriangleType::Lower ){
        m_type = TriangleType::Upper;
      } else {
        m_type = TriangleType::Lower;
      }
    }

    /***********************************************
    *              Determinant                   *
    * For triangular matrix determinant can be 
    * calculated in a simplified way as the
    * product of diagonal elements.
    * The determinant is always of type T, as the
    * computation involves only multiplication of
    * diagonal elements and does not introduce any
    * sign changes.
    ***********************************************/
    T det() const {
      T result = T{1};
      for (size_t diagonal = 0; diagonal < m_dim; diagonal++){
        result *= (*this)( diagonal, diagonal );
      }
      return result;
    }

  };  //TriangularMatrix


}  //namespace