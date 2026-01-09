/***********************************************
 *              DESCRIPTION                    *
 * This file contains a definition of the Matrix
 * class in the namespace sd. The class 
 * implements the most general case of matrix:
 * two-dimensional rectangular matrix along with
 * appropriate basic operations.
 * 
 * It is a template class and because of it, it
 * must contain all definitions in .hpp file
 * (according to prof. Cyganek, sec. 3.15.5).
 ***********************************************/

#pragma once

#include <type_traits> //is_arithmetic_v
#include <vector>
#include <iomanip>     //for better-looking output stream 
#include <stdexcept>

namespace sd {


  template <typename T>
  requires std::is_arithmetic_v<T>
  class Matrix {

  private:
    /***********************************************
    * Matrix is implemeted as a vector of T-type
    * vectors. This representation supports all
    * arithmetic types handled by the class and 
    * provides user-friendly access to matrix
    * elements using matA(row, column).
    * 
    * Memory management is handled automatically
    * by std::vector, which takes care of allocation
    * and deallocation of resources.
    * 
    * Note that indexing starts from zero, i.e.
    * the first element of the first row has index
    * matA(0, 0).
    ***********************************************/
  protected:
    using TVector = std::vector<T>;         
    using TMatrix = std::vector<TVector>;

    size_t m_rows, m_columns;  // number of rows and columns in matrix
    TMatrix m_matrix;       // object containg matrix's data

    unsigned int elementDisplayWidth = 6;  // default value for << operator

  public:

    /***********************************************
    *              CONSTRUCTORS                   *
    * No. 1: Default constructor 
    * Creates a scalar (1x1) matrix with  a 
    * default-initialized T value.
    * 
    * Ex.:
    *    Matrix<double> matA;
    * 
    * Result:
    *     
    *   matA = [ 0.0 ]
    *          
    * No 2: Parametric constructor
    * Constructor initializing  a matrix with a given 
    * numer of rows, colums and initial value.
    * 
    * Ex.:
    *    Matrix matA(3, 2, 1);
    * or better (with explicit type specification)
    *      Matrix<int> matA(3, 2, 1);
    *      
    * Result:
    *          [ 1 1 ]
    *   matA = | 1 1 |
    *          [ 1 1 ]
    * 
    * Note: the dimensions argument are of type int.
    * Providing a negative value would result in an
    * implicit conversion to a very large size_t
    * value, potentially causing excessive memory
    * allocation and program termination.
    * 
    * No 3: The copy constructor is the default one  
    * and requires no explicit implementation.
    ***********************************************/

    Matrix( int rows = 1, int columns = 1, T value = T{} ){
      
      if ( (rows <= 0) || (columns <= 0) ){
        throw std::invalid_argument( "Matrix dimensions have to be greater than 0! ");
      }

      m_rows = static_cast<size_t>( rows );   // explicit type conversion and dimension setup
      m_columns = static_cast<size_t>( columns );  
      m_matrix = TMatrix( m_rows, TVector( m_columns, value ) ); // creating a matrix 
    }




    /***********************************************
    *               GETTERS                       *
    * No 1: Getting the value of a matrix element 
    * specified by its row and column indices.
    * 
    * Ex.:
    *  Given Matrix:
    *          [ 1 2 ]
    *   matA = | 3 4 |
    *          [ 5 6 ]
    * 
    *  Calling:
    *    matA.getValue(2, 1);
    * 
    *  rerurns T-type value (here int): 6
    ***********************************************/

    T getValue( size_t row, size_t column ) const {

      //Range validation: indices must be within matrix bounds

      if ( (row >= m_rows) || (column >= m_columns) ){
        throw std::out_of_range( "Out of range!" );
      }

      return m_matrix[row][column];      
    }

    /***********************************************
    * No 2: Getting an entire row.
    * 
    * Ex.:
    * Given Matrix:
    *          [ 1 2 ]
    *   matA = | 3 4 |
    *          [ 5 6 ]
    * 
    *  Calling:
    *    matA.getRow(1);
    * 
    *  returns T-type vector (here vector<int>):
    *  {3, 4}
    ***********************************************/

    TVector getRow( size_t row ) const {

      if ( row >= m_rows ){
        throw std::out_of_range( "Out of range!" );  
      }

      return m_matrix[row];
    }

    /***********************************************
    * No 3: Getting an  entire column.
    * 
    * Ex.:
    *  Given matrix:
    *          [ 1 2 ]
    *   matA = | 3 4 |
    *          [ 5 6 ]
    * 
    *  Calling
    *    matA.getColumn(1);
    * 
    *  returns T-type vector (int):
    *  {2, 4, 6}
    ***********************************************/

    TVector getColumn( size_t column ) const {

      if ( column >= m_columns ){
        throw std::out_of_range( "Out of range!" );  
      }

      TVector temp {};    // creating an empty TVector

      for ( size_t row = 0; row < m_rows; row++ ){
        temp.push_back( m_matrix[row][column] ); //filling it with values
      }

      return temp;  
    }

    /***********************************************
    *         Getting matrix dimensions           *
    * Fuction returns a 2-elements vector containing
    * number of rows and columns of matrix.
    * 
    * Two others return only one dimension.
    ***********************************************/

    std::vector<size_t> getDimensions() const {
      return std::vector<size_t> { m_rows, m_columns };
    }

    size_t getRowDim() const {
      return m_rows;
    }
    size_t getColumnDim() const {
      return m_columns;
    }

    /***********************************************
    *               SETTERS                       *
    * No. 1: Setting the value of a matrix element
    * specified by its row and column indices.
    * 
    * Ex.:
    * Given matrix:
    *          [ 1.0 2.0 ]
    *   matA = | 3.0 4.0 |
    *          [ 5.0 6.0 ]
    * 
    * Calling:
    *   matA.setValue(2, 1, 5.7);
    * 
    *  result in:
    *          [ 1.0 2.0 ]
    *   matA = | 3.0 4.0 |
    *          [ 5.0 5.7 ]
    ***********************************************/

    void setValue( size_t row, size_t column, T value ){

      //Validation
      if ( (row >= m_rows) || (column >= m_columns) ){
        throw std::out_of_range( "Out of range!" );
      }

      m_matrix[row][column] = value;
    }

    /***********************************************
    * No 2: Setting an entire row.
    * 
    * Ex.:
    * Given matrix:
    *          [ 1.0 2.0 ]
    *   matA = | 3.0 4.0 |
    *          [ 5.0 6.0 ]
    * 
    *  Executing:
    *    std::vector<double> v1{1.1, 2.2};
    *    matA.setRow(1, v1);
    * 
    *  result in 
    *          [ 1.0 2.0 ]
    *   matA = | 1.1 2.2 |
    *          [ 5.0 5.7 ]
    ***********************************************/

    void setRow( size_t row, const TVector & vector ){

      if ( row >= m_rows ){
        throw std::out_of_range( "Out of range!" );  
      }

      if ( vector.size() != m_columns ){
        throw std::invalid_argument( "Rows dimensions mismatch!" );
      }

      m_matrix[row] = vector;
    }

    /***********************************************
    * No 3: Setting whole column.
    * 
    * Ex.:
    * Given matrix:
    *          [ 1.0 2.0 ]
    *   matA = | 3.0 4.0 |
    *          [ 5.0 6.0 ]
    * 
    *  Executiong code:
    *    std::vector<double> v1{1.1, 2.2, 3.3};
    *    matA.setColumn(0, v1);
    * 
    *  results in
    *          [ 1.1 2.0 ]
    *   matA = | 2.2 4.0 |
    *          [ 3.3 6.0 ]
    ***********************************************/
  
    void setColumn( size_t column, const TVector & vector ){

      if ( column >= m_columns ){
        throw std::out_of_range( "Out of range!" );  
      }
      if ( vector.size() != m_rows ){
        throw std::invalid_argument( "Column dimensions missmatch!" );
      }

      for ( size_t row = 0; row < m_rows; row++ ){
        m_matrix[row][column] = vector[row];
      }
    }

    /***********************************************
    *      Setting display width of elements      *
    * 
    * It provides better-looking and configurable
    * display option for matrix using << operator.
    * User can set with it how wide each element is.
    * That provides alignment of elements.
    * Default value is 6.
    ***********************************************/  

    void setElementDisplayWidth( int value ){

      if ( value <= 0 ){
        throw std::invalid_argument( "Display width has to be bigger than 0!" );
      }

      elementDisplayWidth = static_cast<unsigned int>(value);
    }


    /***********************************************
    *            Matrix operations                *
    * 
    * The following comments provide short examples
    * of each operation along with the mathematical
    * assumptions that must be satisfied.
    ***********************************************/
    /***********************************************
    *            Matrix transposition              *
    *              matA.transpose()
    * 
    * Any matrix can be transposed. 
    * If matA is m x n matrix then matA_trasposed 
    * is n x m matrix such that:
    * 
    * [matA_trasposed]_{ij} = [matA]_{ji}
    ***********************************************/
    Matrix transpose() const {

      Matrix result( m_columns, m_rows, T{0} ); 

      for ( size_t row = 0; row < result.m_rows; row++ ){
        for ( size_t column = 0; column < result.m_columns; column++ ){
          result( row, column ) = (*this)( column, row );
        }
      }

      return result; 
    }

    /***********************************************
    *                 Operators                   *
    * The following operators implement basic
    * matrix operations. Short examples and the
    * required mathematical assumptions are
    * provided in the comments.
    *
    * To allow operations between matrices of
    * different arithmetic types, some operators
    * are implemented as templates. For example,
    * it is possible to add a matrix of type double
    * to a matrix of type int.
    ***********************************************/

    /***********************************************
    *                operator +=                   *
    *               matA += matB                   *
    * 
    * Mathematical assumptions:
    *   1) dim(matA) = dim(matB)
    * 
    * Description of operation:
    * Elements of both martixcs with the same 
    * indices are added.
    * 
    * The diferent types of matrixes can be used,
    * but the result will be always cast to the 
    * left hand side matrix type.
    * 
    * Ex.:
    *          [ 1.0 2.0 ]
    *   matA = | 3.0 4.0 |
    *          [ 5.0 6.0 ]
    *
    *          [ 6.1 5.2 ]
    *   matB = | 4.3 3.4 |
    *          [ 2.5 1.6 ]
    *
    *          matA += matB;  
    *
    *          [ 7.1 7.2 ]
    *   matA = | 7.3 7.4 |
    *          [ 7.5 7.6 ]
    ***********************************************/  

    template <typename T2>
    requires std::is_arithmetic_v<T2>
    Matrix<T> & operator += ( const Matrix<T2> & other){

      if ( (m_rows != other.getRowDim()) || ( m_columns != other.getColumnDim()) ){
        throw std::invalid_argument("Matrix dimensions mismatch");
      }

      for ( size_t row = 0; row < m_rows; row++ ){
        for ( size_t column = 0; column < m_columns; column++ ){
            (*this)( row, column) += static_cast<T>( other( row, column ) );
        }
      }
      
      return *this;
    }

    /***********************************************
    *                operator -=                   *
    *               matA -= matB                   *
    * 
    * Mathematical assumptions:
    *   1) dim(matA) = dim(matB)
    * 
    * Description of operation:
    * Elements of both martices with the same 
    * indices are subtracted. 
    * 
    * The diferent types of matrixes can be used,
    * but the result will be always cast to the 
    * left hand side matrix type.
    * 
    * Ex.:
    *          [ 1.0 2.0 ]
    *   matA = | 3.0 4.0 |
    *          [ 5.0 6.0 ]
    *
    *          [ 6.1 5.2 ]
    *   matB = | 4.3 3.4 |
    *          [ 2.5 1.6 ]
    *
    *          matA -= matB;  
    *
    *          [ 5.1  3.2 ]
    *   matA = | 1.3 -1.4 |
    *          [-2.5 -5.6 ]
    ***********************************************/  

    template <typename T2>
    requires std::is_arithmetic_v<T2>
    Matrix<T> & operator -= ( const Matrix<T2> & other ) {

      if ( (m_rows != other.getRowDim()) || ( m_columns != other.getColumnDim()) ){
        throw std::invalid_argument(" Matrices dimensions mismatch!" );
      }
      
      for ( size_t row = 0; row < m_rows; row++ ){
        for ( size_t column = 0; column < m_columns; column++ ){
          (*this)( row, column) -= static_cast<T>( other( row, column ) );
        }
      }
      return *this;  
    }

    /***********************************************
    *                operator *=                   *
    * For matrix multiplication this operator has 
    * sense only for square matrix due to 
    * dimmensions restrictions. 
    * It has thou sense for scalar multiplication.       
    ***********************************************/  

    template <typename S>
    requires std::is_arithmetic_v<S>  
    Matrix & operator *= ( const S scalar ) {

      for ( auto & row : m_matrix ){
        for ( auto & element : row ){
          element *= static_cast<T>(scalar) ;
        }
      }

      return *this;  
    }

    /***********************************************
    *                operator ()                   *
    *                matA(r,c)
    * Operator () allows to get elements value
    * with indices r and c.
    * Two versions for read and write are necessary.     
    ***********************************************/  

    T operator () (size_t r, size_t c) const { 
      return m_matrix.at(r).at(c); 
    } 

    T& operator () (size_t r, size_t c){ 
      return m_matrix.at(r).at(c); 
    }

    /***********************************************
    *                operator <<                   *
    *                
    * Allows a matrix to be sent to an output stream.
    * The matrix is always printed with a leading
    * and trailing newline.
    *
    * To provide a readable and aligned output, the
    * user may adjust the display width of elements
    * using the setElementDisplayWidth() function.
    * The default width is 6 characters.     
    ***********************************************/  
    
    friend std::ostream& operator<<( std::ostream& os, const Matrix & matrix ){

      os << "\n"; // always begins from new line

      for ( size_t row = 0; row < matrix.getRowDim(); ++row ) {
        for ( size_t column = 0; column < matrix.getColumnDim(); ++column ) {
          os << std::setw(matrix.elementDisplayWidth) << matrix( row, column ) << " ";
        } 
        os << "\n";
      }

      return os;
    }

  };  //class
  
  
  /***********************************************
  *                operator +                   *
  *            matC = matA + matB               *
  * 
  * Mathematical assumptions:
  *   1) dim(matA) = dim(matB)
  * 
  * Description of operation:
  * Elements of both martices with the same 
  * indices are added.
  * 
  * The resulting matrix has a value type equal 
  * to the common arithmetictype of both 
  * operand matrices.
  * 
  * Ex.:
  *          [ 1.0 2.0 ]
  *   matA = | 3.0 4.0 |
  *          [ 5.0 6.0 ]
  *
  *          [ 6.1 5.2 ]
  *   matB = | 4.3 3.4 |
  *          [ 2.5 1.6 ]
  *
  *          [ 7.1 7.2 ]
  *   matC = | 7.3 7.4 |
  *          [ 7.5 7.6 ]
  ***********************************************/    
    
  template <typename T1, typename T2>
  requires std::is_arithmetic_v<T1> && std::is_arithmetic_v<T2>
  auto operator + ( const Matrix<T1> & lhs, const Matrix<T2> & rhs ){

    if ( (lhs.getRowDim() != rhs.getRowDim()) || ( lhs.getColumnDim() != rhs.getColumnDim() ) ){
      throw std::invalid_argument( "Matrices dimensions mismatch!" );
    }
    
    using ResultType = std::common_type_t<T1, T2>;  //choose 'stronger' type, eg. double over int
    Matrix<ResultType> result( lhs.getRowDim(), lhs.getColumnDim(), ResultType{0} );
    
    for ( size_t row = 0; row < lhs.getRowDim(); row++ ){
      for ( size_t column = 0; column < lhs.getColumnDim(); column++ ){
        result( row, column ) = static_cast<ResultType>( lhs( row, column ))
                              + static_cast<ResultType>( rhs( row, column ));
      }
    }

    return result;  
  }

   /***********************************************
  *                operator -                   *
  *            matC = matB - matA               *
  * 
  * Mathematical assumptions:
  *   1) dim(matA) = dim(matB)
  * 
  * Description of operation:
  * Elements of the right-hand side matrix are 
  * subtracted from the corresponding elements 
  * of the left-hand side matrix.
  * 
  * The resulting matrix has a value type equal 
  * to the common arithmetictype of both 
  * operand matrices.
  * 
  * Ex.:
  *          [ 1.0 2.0 ]
  *   matA = | 3.0 4.0 |
  *          [ 5.0 6.0 ]
  *
  *          [ 6.1 5.2 ]
  *   matB = | 4.3 3.4 |
  *          [ 2.5 1.6 ]
  *
  *          [ 5.1  3.2 ]
  *   matC = | 1.3 -1.4 |
  *          [-2.5 -5.6 ]
  ***********************************************/    

  template <typename T1, typename T2>
  requires std::is_arithmetic_v<T1> && std::is_arithmetic_v<T2>
  auto operator - ( const Matrix<T1> & lhs, const Matrix<T2> & rhs ){

    if ( (lhs.getRowDim() != rhs.getRowDim()) || ( lhs.getColumnDim() != rhs.getColumnDim() ) ){
      throw std::invalid_argument( "Matrices dimensions mismatch!" );
    }
    
    using ResultType = std::common_type_t<T1, T2>;  //choose 'stronger' type, eg. double over int
    Matrix<ResultType> result( lhs.getRowDim(), lhs.getColumnDim(), ResultType{0} );

    for ( size_t row = 0; row < lhs.getRowDim(); row++ ){
      for ( size_t column = 0; column < lhs.getColumnDim(); column++ ){
        result( row, column ) = static_cast<ResultType>( lhs( row, column ))
                              - static_cast<ResultType>( rhs( row, column ));
      }
    }

    return result;  
  }
  
  /***********************************************
  *                operator *                   *
  *           matC = matA * matB                *
  * 
  * Mathematical assumptions:
  *   1) Dimensions:
  *      matA is an m x n matrix
         matB is an n x p matrix

  * Description of operation:
  *      matC is an m x p matrix such that
  *      c_ij = sum_{k=0}^(n-1) a_{ik} * b_{kj}
  *      for i = 0, ... , m-1
  *          j = 0, ... , p-1
  * Ex.:
  *          [ 1.0 2.0 ]
  *   matA = | 3.0 4.0 |
  *          [ 5.0 6.0 ]
  *
  *          [ 1.1 1.2  1.3]
  *   matB = [ 2.1 2.2  2.3]
  *          
  *
  *          matC = matA * matB;  
  *
  *          [ 5.3   5.6   5.9 ]
  *   matC = | 11.7  12.4  13.1|
  *          [ 18.1  19.2  20.3]
  ***********************************************/  

  template <typename T1, typename T2>
  requires std::is_arithmetic_v<T1> && std::is_arithmetic_v<T2>
  auto operator * ( const Matrix<T1> & lhs, const Matrix<T2> & rhs ){

    if ( lhs.getColumnDim() != rhs.getRowDim() ){
      throw std::invalid_argument( "Matrices dimensions mismatch!" );
    }

    using ResultType = std::common_type_t<T1, T2>;  //choose 'stronger' type, eg. double over int
    Matrix<ResultType> result( lhs.getRowDim(), rhs.getColumnDim(), ResultType{0} );

    for ( size_t row = 0; row < lhs.getRowDim(); row++ ){
      for ( size_t column = 0; column < rhs.getColumnDim(); column++ ){
        for ( size_t k = 0; k < lhs.getColumnDim(); k++ ){          
          result( row, column) += static_cast<ResultType>( lhs ( row, k ) )
                               *  static_cast<ResultType>( rhs( k, column) );
        }
      }
    }

    return result;  
  }

  /***********************************************
  *                operator *                   *
  *           matC = matA * scalar              *
  *           matC = scalar * matA              *
  * 
  * No mathematical assumptions.
  * 
  * Description of operation:
  *    matC is a matrix where every element is
  *    product of matA elements with scalar           
  * Operation is (unlike matrix multiplication) 
  * commutative and scalar can be of any arithmetic
  * type.
  * 
  * Ex.:
  *          [ 1.0 2.0 ]
  *   matA = | 3.0 4.0 |
  *          [ 5.0 6.0 ]
  *
  *   scalar = 2.0;
  *          
  *          matC = scalar * matA;  
  *
  *          [ 2.0   4.0 ]
  *   matC = | 6.0   8.0 |
  *          [ 10.0  12.0]
  ***********************************************/
  template <typename T, typename S>
  requires std::is_arithmetic_v<T> && std::is_arithmetic_v<S>
  auto operator * ( const Matrix<T> & matrix, const S scalar ){

    using ResultType = std::common_type_t<T, S>;  //choose 'stronger' type, eg. double over int
    Matrix<ResultType> result( matrix.getRowDim(), matrix.getColumnDim(), ResultType{0} );

    for ( size_t row = 0; row < matrix.getRowDim(); row++ ){
      for ( size_t column = 0; column < matrix.getColumnDim(); column++ ){
        result ( row, column ) = static_cast<ResultType>( matrix( row, column )) 
                               * static_cast<ResultType>( scalar ) ;
      }
    }
    return result;  
  }

  template <typename T, typename S>
  requires std::is_arithmetic_v<T> && std::is_arithmetic_v<S>
  auto operator * ( const S scalar, const Matrix<T> & matrix ){

    return matrix * scalar;

  }

} // namespace sd

