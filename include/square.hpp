/***********************************************
 * DESCRIPTION
 *
 * This file contains the definition of the
 * SquareMatrix class in the sd namespace.
 * SquareMatrix is a class derived from Matrix
 * and represents a two-dimensional square matrix,
 * i.e. a matrix with an equal number of rows
 * and columns.
 ***********************************************/

#pragma once

#include "matrix.hpp"


namespace sd{

  static constexpr double eps = 1e-12;   // near 0.0 definition

  /***********************************************
  * Determinant result type policy
  *
  * The determinant of a matrix may take negative
  * values even if the matrix element type is
  * unsigned. Therefore, for unsigned integral
  * matrix element types, the return type of the
  * determinant is promoted to the corresponding
  * signed type.
  *
  * For all other arithmetic types (signed integral
  * and floating-point), the return type remains
  * unchanged.
  *
  * This trait defines the return type of the
  * determinant in a safe and mathematically
  * consistent way.
  ***********************************************/
  template <typename T, typename = void>  //if type is proper for determinant leave it
  struct determinant_result {
    using type = T;
  };

  template <typename T>
  struct determinant_result<
    T,
    std::enable_if_t<std::is_integral_v<T> && std::is_unsigned_v<T>>  //for unsigned int change type to signed
  > {
    using type = std::make_signed_t<T>;
  };

  template <typename T>
  using determinant_result_t = typename determinant_result<T>::type;

  /***********************************************
    *                 CLASS                    *
  ***********************************************/
  template <typename T>
  requires std::is_arithmetic_v<T>
  class SquareMatrix : public Matrix<T> {
  private:
  
  protected:
    size_t m_dim;
  public:
    /***********************************************
    *              CONSTRUCTORS                   *
    * 
    * No 1: Default constructor 
    * Creatse a 1×1 matrix with a default-initialized 
    * T value.
    * 
    * Ex.:
    *    SquareMatrix<double> matA;
    * Result:
    *     
    *   matA = [ 0.0 ]
    *          
    * No 2: Intializing a matrix with a given value
    * and dimension.
    * 
    * Ex.:
    *    SquareMatrix matA(3, 1);
    * or better (with explicit type specification)
    *     SquareMatrix<int> matA(3, 1);
    * 
    * Result:
    *          [ 1 1 1 ]
    *   matA = | 1 1 1 |
    *          [ 1 1 1 ]
    * 
    * Note: the dimensions argument are int, because
    * by giving a negative value it will be convert to
    * very large size_t value and programm will try to
    * alocate a big amount of memory - what propably
    * will cause a termination of programm.  
    * No 3: Copy constructor is the default one - it 
    * need no implementation.
    ***********************************************/

    SquareMatrix( int dim = 1, T value = T{} ) 
      : Matrix<T>( dim, dim, value), m_dim( static_cast<size_t>( dim ) ) {}

    /***********************************************
    * No 4: Costructor converting Matrix 
    * into SquareMatrix. 
    ***********************************************/
    
    SquareMatrix( const Matrix<T> & matrix) 
            : Matrix<T>( matrix ), m_dim( matrix.getRowDim() ) {

      if ( matrix.getRowDim() != matrix.getColumnDim() ){
        throw std::invalid_argument( "Dimension error!" );
      }
    
    }
    
    /***********************************************
    *            Matrix operations                *
    * 
    * The following comments provide short examples
    * of each operation along with the mathematical
    * assumptions that must be satisfied.
    ***********************************************/
    /***********************************************
    * Transposition
    *
    * For a square matrix, transposition does not
    * change the matrix dimensions. Therefore,
    * two versions of the transpose operation are
    * provided:
    *
    * 1) transpose()        – returns a new matrix
    * 2) transpose_in_place() – transposes the matrix
    *                           in place
    ***********************************************/

    SquareMatrix transpose() const {
      
      return SquareMatrix<T> (  Matrix<T>::transpose() );
      
    }

    void transpose_in_place() {

      for ( size_t row = 0; row < m_dim; row++ ) {
        for ( size_t column = row + 1; column < m_dim; column++) {  //iteration only in upper triangle is sufficent
            std::swap((*this)(row, column), (*this)(column, row));
        }
      }

    }

    /***********************************************
    *                   Minor                      *
    *                matA.minor( r, c )
    * If A is an n × n matrix, then its minor is an
    * (n−1) × (n−1) matrix obtained from A by removing
    * row r and column c (0-based indexing).
    * 
    * Note: The minor of a 1×1 matrix is undefined and
    * calling this function for such a matrix will
    * result in an exception.
    * 
    * Ex.:
    *          [ 1 2 3 ]
    *   matA = | 4 5 6 |
    *          [ 7 8 9 ]
    * 
    * 
    *                    [ 1 2 ]
    *   matA.min(1, 2) = [ 7 8 ]
    *                    
    ***********************************************/

    SquareMatrix minor( size_t row, size_t column ) const {

      if ( (row >= m_dim) || (column >= m_dim) ){
        throw std::out_of_range( "Out of range!" );
      }

      SquareMatrix result( m_dim - 1 ); 

      size_t rr = 0;    //result matrix row index
      for( size_t r = 0; r < m_dim; r++ ){
        
        if (r == row){
            continue;     //skip given row
        }
        
        size_t cc = 0;  //result matrix column index
        for( size_t c =0; c < m_dim; c++ ){

          if (c == column){
            continue;     //skip given column
          }

          result( rr, cc ) = (*this)( r, c );
          cc++;     //increment column index
        }

        rr++;     //increment row index
      } 

      return result;
    }

    /***********************************************
    *                 Determinant                  *
    *                  matA.det()
    * 
    * If A is n x n matix, then determinant can be 
    * computed as follow:
    * 
    * det(A) = sum_{j=0}^{n-1} (-1)^j*a_{0j}*det(M_{0j}) 
    * 
    * where M_{0j} is the minor obtained by removing
    * row 0 and column j.
    * 
    * Note: The determinant can be expanded along any row
    * or column. In this implementation, expansion
    * is always performed along the first row.
    * 
    * The given algorithm is not efficient for large
    * matrices. 
    ***********************************************/
    
    auto det() const {
      
      using ResultType = determinant_result_t<T>; //determinant type can be different from matrix type

      //special cases for dim = 1 and dim = 2:
      if ( m_dim == 1 ) {
        return static_cast<ResultType>( (*this)(0, 0) );
      }

      if (m_dim == 2) {
        return static_cast<ResultType>(
               (*this)( 0, 0 ) * (*this)( 1, 1 )
               -(*this)( 0, 1 ) * (*this)( 1, 0 )
              );
      }

      ResultType result {};

      for ( size_t column = 0; column < m_dim; column++ ){
        
        ResultType sign =  ( column % 2 == 0 ) ? ResultType{1} : ResultType{-1}; //even with +, odd with -
        SquareMatrix<T> minorMat = minor( 0, column );
        result +=  sign *  static_cast<ResultType>( (*this)( 0, column ) )
                        *  static_cast<ResultType>( minorMat.det() );
      }

      return result;
    }

    /***********************************************
    *             Cofactor Matrix               *
    *                matA.cofactor()
    * The cofactor matrix of an n x n matrix is an n x n
    * matrix defined as follows:
    * 
    * C_{ij} = (-1)^(i+j) * det( Minor_{ij} )
    * 
    ***********************************************/

    auto cofactor() const {

      using ResultType = determinant_result_t<T>; //determinant type can be different from matrix type
      SquareMatrix<ResultType> result( m_dim, ResultType{0} );

      for (size_t row = 0; row < m_dim; row++ ){
        for (size_t column = 0; column < m_dim; column++ ){
          ResultType sign =  ( (row + column ) % 2 == 0 ) ? ResultType{1} : ResultType{-1}; //even with +, odd with  
          result( row, column ) = sign * minor(row, column).det(); 
        }
      }

      return result;
    }

    /***********************************************
    *             Matrix inversion               *
    *                matA.inv()
    * A matrix B is inverse to matrix A, if 
    *     matA * matB = matB * matA = I_n 
    * whre I_n is n-dimensional identitty matrix.
    * Matrix B is often designated as (matA)^(-1).
    * 
    * It can be calculated as:
    *   matB = 1/det(matA) * C^T
    * where C^T is trasposed cofactor matrix.
    *   
    * Assumption:
    *   1) The matrix must be invertible, i.e.
    *      its determinat must be non-zero.
    * 
    * Note: The inverse matrix is returned as a matrix
    * of type double, even if the original matrix
    * has an integral element type; due to 1/det 
    * operation inversion generally produces 
    * non-integer values.
    ***********************************************/

    SquareMatrix<double> inverse() const {
      
      auto detA = det(); //type will be specified thanks to det() function, but it will be cast to double 
      
      if ( std::abs( (static_cast<double>( detA )) ) <= eps ){
        throw std::invalid_argument( "Math error. Matrix is not invertible!" );
      }

      return ( 1.0 / (static_cast<double>( detA )) ) * ( ( cofactor() ).transpose() );
    }



  };    //class 



} //namespace