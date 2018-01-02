# ifndef __DENSE_MATRIX_H__
# define __DENSE_MATRIX_H__

# include "Matrix.H"


namespace mg {
                namespace numeric {
                                    namespace algebra {

// forward declaration 
template <typename Type>
class DenseMatrix ;

//--

template<typename U>
std::ostream& operator<<(std::ostream& os, const DenseMatrix<U>& m );

template<typename U>
std::vector<U> operator*(const DenseMatrix<U>& , const std::vector<U> ) ;

template<typename U>
DenseMatrix<U> operator*(const DenseMatrix<U>& , const U& rhs) noexcept ;

template<typename U>
DenseMatrix<U> operator/(const DenseMatrix<U>& , const U& rhs) noexcept ;


//using the default (ijk) algorithm time-complexity = O(N^3) 
template<typename U> 
DenseMatrix<U> operator* (const DenseMatrix<U>& , const DenseMatrix<U>&) ;  

template<typename U>
DenseMatrix<U> operator+(const DenseMatrix<U>& , const DenseMatrix<U>&) ;

template<typename U>
DenseMatrix<U> operator-(const DenseMatrix<U>& , const DenseMatrix<U>&) ;


template<typename U>
auto _det(const DenseMatrix<U>& a ) -> U ;

//perform mat-mat product using  Strassen algorithm time-complexity = O(N^(2.807))
template <typename U>
void strassen(const DenseMatrix<U>& , const DenseMatrix<U>&, 
                    DenseMatrix<U>&,  const std::size_t tam ) ;

//-----------------------------------------------------------------------------
// utility function for strassen algorithm
template <typename U>
auto sum(const DenseMatrix<U>& A , const DenseMatrix<U>& B ,
               DenseMatrix<U>& C, std::size_t tam ) noexcept;

template <typename U>
auto subtract(const DenseMatrix<U>& A , const DenseMatrix<U>& B ,
                    DenseMatrix<U>& C, std::size_t tam ) noexcept;



/**------------------------------------------------------------------------------
 * \class DenseMatrix
 * @brief Dense Matrix Class 
 *   
 *    using std::vector<T> for storing the whole matrix data
 *    
 *
 *
 *    
 *  @Marco Ghiani Dec 2017 , Glasgow UK
 * 
 ------------------------------------------------------------------------------*/



template <typename Type>
class DenseMatrix 
                           :     public Matrix<Type> 
{
      

//-- Non member function (friend)
//
       template<typename U>
       friend std::ostream& operator<<( std::ostream& os, const DenseMatrix<U>& m   );

       template<typename U>
       friend std::vector<U> operator*(const DenseMatrix<U>& , const std::vector<U> ) ;
 
       template<typename U>
       friend DenseMatrix<U> operator*(const DenseMatrix<U>& , const DenseMatrix<U>&) ;  

       template <typename U>
       friend void strassen(const DenseMatrix<U>& , const DenseMatrix<U>&, 
                                  DenseMatrix<U>& , const std::size_t tam ) ;

       template<typename U>
       friend DenseMatrix<U> operator+(const DenseMatrix<U>& , const DenseMatrix<U>&) ;
      
       template<typename U>
       friend DenseMatrix<U> operator-(const DenseMatrix<U>& , const DenseMatrix<U>&) ;

       template<typename U>
       friend DenseMatrix<U> operator*(const DenseMatrix<U>& , const U& rhs) noexcept ;

       template<typename U>
       friend DenseMatrix<U> operator/(const DenseMatrix<U>& , const U& rhs) noexcept ;
      
       template<typename U>
       friend auto _det(const DenseMatrix<U>& a ) -> U ;     
      

//--
//
   public:
     
    //-- constructor  
       constexpr DenseMatrix (std::initializer_list<std::vector<Type>> ) noexcept ;
       
       constexpr DenseMatrix (const std::string& );
       
       constexpr DenseMatrix (std::size_t , std::size_t) noexcept ;
      
       constexpr DenseMatrix(std::size_t) noexcept ;

       virtual  ~DenseMatrix() = default ;
      

    // - method   
       void constexpr print () const noexcept override ;

       auto constexpr size1() const noexcept { return Rows ; }

       auto constexpr size2() const noexcept { return Cols ; }
       
       auto constexpr isSquare() const noexcept { return (Rows == Cols) ; } 

       Type constexpr findValue(const std::size_t , const std::size_t ) const noexcept ;

       std::vector<Type> diag() const noexcept ;      
       
       constexpr DenseMatrix<Type> exctractMinor(std::size_t r, std::size_t c) ;

       constexpr auto det()-> Type ;
   
       constexpr DenseMatrix<Type> Minor(std::size_t r, std::size_t c);
       
       constexpr DenseMatrix<Type> minors(std::size_t r, std::size_t c) const;

   //-  operators 
   //
   //
       Type& operator()(const std::size_t , const std::size_t) noexcept override;

       const Type& operator()(const std::size_t , const std::size_t ) const noexcept override;
 
       DenseMatrix<Type>& operator=(const DenseMatrix<Type>& that) noexcept ; 
       
       DenseMatrix<Type>& operator+=(const DenseMatrix<Type>& rhs ) ;
       
       DenseMatrix<Type>& operator-=(const DenseMatrix<Type>& rhs ) ;
       
       DenseMatrix<Type>& operator*=(const Type& rhs ) ;

       DenseMatrix<Type>& operator/=(const Type& rhs ) ;     


//---
   protected:
      
       std::vector<Type> data ;
       
       std::size_t Rows ;
       std::size_t Cols ;
      
       mutable Type dummy ;
       std::size_t nnz    ; // number of non zero elem (for eval. degree of density)
       const std::size_t leafSize = 1;    // variable used in strassen alghorithm 
      
       Type zero = 0.0 ;
} ;

// template aliasing 
template <typename T> 
using matrix = DenseMatrix<T> ;
     


//-------------------------------        Implementation      -----------------------------------------



template<typename T>
constexpr DenseMatrix<T>::DenseMatrix (std::initializer_list<std::vector<T>> dense ) noexcept : nnz{0}
{
    Rows = dense.size()  ;
    Cols = (*dense.begin()).size() ;  

    data.resize(Rows*Cols);
    
    auto i=0; 
    auto j=0;  
    for(auto& row : dense )
    {
        for(auto& col : row)  
        {
            data.at(i) = col ;
            i++;
            if(col !=0 ) nnz++ ;
        }
    }
      
}

//default constructor
template <typename T>
constexpr DenseMatrix<T>::DenseMatrix(const std::size_t row,
                                      const std::size_t col) noexcept  
                                                                        : Rows{row}, Cols{col}, nnz{0}
{
      data.resize(Rows*Cols);
}


// construct square Identity Matrix

template<typename T>
constexpr DenseMatrix<T>::DenseMatrix(std::size_t n) noexcept : Rows{n},
                                                                Cols{n},
                                                                nnz{n}
{          
      data.resize(n*n) ;
      for(auto i=1; i<=Rows ; i++)
          data.at((i-1)*Cols + (i-1)) = 1.0;  

}


template <typename T>
constexpr DenseMatrix<T>::DenseMatrix(const std::string& filename)  
{
    std::ifstream f(filename , std::ios::in);  
    
    if(!f)
    {
       std::string mess = "Error opening file \'" + filename + 
                          "\'\n>> Exception Thrown in DenseMatrix constrsuctor <<" ; 
       throw OpeningFileException(mess.c_str());    
    }
    if( filename.find(".mtx") != std::string::npos )
    {
       std::string line ;
       T elem = 0 ;
       auto r = 0 , c = 0 ;       
       auto i=0 , j=0 ;
       getline(f,line) ; // jump out the header
       
       while(getline(f,line))
       {
          std::istringstream ss(line) ;  
          if(i==0)
          {
             ss >> Rows ; 
             ss >> Cols ;
             ss >> nnz  ;
             
             data.resize(Rows*Cols);
          }
          else
          {
             ss >> r ; 
             ss >> c ;
             ss >> elem ;
             
             this->operator()(r,c) = elem ;

          }
          i++;  
       }
    }
    else
    {
       std::string line, temp ;
       T elem = 0;
       auto i=0 , j=0;
       while(getline (f,line))
       {
          j=0;  
          std::istringstream ss(line);  
          while( ss >> elem){
             data.push_back(elem);
          j++;
             if(elem != 0) nnz++ ;
          }          
          i++ ;
       }
       Rows = i;     
       Cols = j;     
    }
}


// assignament operator
//
template<typename T>
DenseMatrix<T>& DenseMatrix<T>::operator=(const DenseMatrix<T>& that) noexcept 
{                                                                         
   data.resize(that.data.size()); 
   if(this != &that) 
   {
      auto i=0 ;
      for(auto x : that.data )
      {
            //std::cout << i << std::endl;
            data.at(i) = x ;
            i++;
      }
       Rows = that.Rows  ;
       Cols = that.Cols  ;
       nnz = that.nnz   ;                   
    }
    else
    {
      std::cerr << "Attempt to perform Self assignement occured in DenseMatrix<T>::operator= " << std::endl;
    }
    
    return *this;
}




template <typename T>
void constexpr DenseMatrix<T>::print() const noexcept 
{
    for(auto i=1; i<= size1() ; i++){
      for(auto j=1 ; j<= size2() ; j++ ){
         std::cout << std::setw(6) << findValue(i,j) << " " ;   
      }
      std::cout << std::endl;
    }  


}

template <typename T>
T constexpr DenseMatrix<T>::findValue(const std::size_t row , const std::size_t col) const noexcept 
{
    assert(row > 0      &&
           row <= Rows  && 
           col > 0      &&
           col <= Cols     );  

   return data.at((row-1)*Cols + (col -1));   
}



template <typename T>
T& DenseMatrix<T>::operator()(const std::size_t row, const std::size_t col) noexcept 
{
    assert(row > 0      &&
           row <= Rows  && 
           col > 0      &&
           col <= Cols     );  

     return data.at((row-1)*Cols + (col -1)); 
}



template <typename T>
const T& DenseMatrix<T>::operator()(const std::size_t row, const std::size_t col) const noexcept 
{
    assert(row > 0      &&
           row <= Rows  && 
           col > 0      &&
           col <= Cols     );  

     return data.at((row-1)*Cols + (col -1)); 
}


template <typename T>       
std::vector<T> DenseMatrix<T>::diag() const noexcept 
{
      
      if(! isSquare())
      {
         std::cerr << "Impssible get diagonal from non-square matrix!" << std::endl;
         exit(-1); // avoid to thrown exception 
      }
      else
      {
            
         std::vector<T> d(Rows);
         for(auto i=1 ; i<= Rows ; i++){
               d.at(i-1) =  this->operator()(i,i) ;         
         }
      
         return d;  
      }
}


// product matrix * scalar 

template<typename T>
DenseMatrix<T>& DenseMatrix<T>::operator*=(const T& rhs )
{     
# pragma omp parallel for 
   for(auto i=1 ; i<= Rows ; i++){
      for(auto j=1 ; j <= Cols ; j++){
            this->operator()(i,j) *= rhs ;
      }
   }   
}

template<typename T>
DenseMatrix<T>& DenseMatrix<T>::operator/=(const T& rhs )
{     
# pragma omp parallel for 
   for(auto i=1 ; i<= Rows ; i++){
      for(auto j=1 ; j <= Cols ; j++){
            this->operator()(i,j) /= rhs ;
      }
   }   
}




//--- 
//
template <typename T>
DenseMatrix<T>& DenseMatrix<T>::operator+=(const DenseMatrix<T>& rhs )
{
      if( this->size1() != rhs.size1() ||
          this->size2() != rhs.size2()    ) 
      {
          throw InvalidSizeException("Matrix dimension doesn't match in operator += !");
      }
      else
      {
# pragma omp parallel for     
         for(auto i=1 ; i<= this->Rows ; i++){
            for(auto j=1 ; j<= this->Cols ; j++){
                  this->operator()(i,j) += rhs(i,j);
                  if(this->operator()(i,j) == 0 && rhs(i,j) !=0 ) nnz++;
            }
         }
      return *this;   
      }
}


//---
//
template <typename T>
DenseMatrix<T>& DenseMatrix<T>::operator-=(const DenseMatrix<T>& rhs )
{
      if( this->size1() != rhs.size1() ||
          this->size2() != rhs.size2()    ) 
      {
          throw InvalidSizeException("Matrix dimension doesn't match in operator += !");
      }
      else
      {
# pragma omp parallel for     
         for(auto i=1 ; i<= this->Rows ; i++){
            for(auto j=1 ; j<= this->Cols ; j++){
                  this->operator()(i,j) -= rhs(i,j);
                  if(this->operator()(i,j) == 0 && rhs(i,j) !=0 ) nnz++;
            }
         }
       return *this;   
      }
}

/*  @fun Extract a minor (from row and column to exclude) 
 *       only for square matrix !
 */ 

template <typename T>
constexpr DenseMatrix<T> DenseMatrix<T>::exctractMinor(std::size_t r, std::size_t c)
{  
   if( ! this->isSquare() )
   {
       throw InvalidSizeException("Error in Minor method >>Matrix must be square<< ");   
       
   }
   else
   {
     auto colCount = 1 , rowCount = 1 ;    
     DenseMatrix<T> res{Rows-1 , Cols-1} ;

       for(auto i=1 ; i <= Rows ; i++ )
       {
         if(i != r ) // when r is not the row choose (r) 
         {
           colCount = 1;  
            
           for(auto j=1 ; j <= Cols ; j++ )
           {
                // when j is not the colums choose (c)
                if(j != c)
                {
                  res(rowCount,colCount) = this->operator()(i,j) ;
                  colCount++;
                }
            
           }
           rowCount ++ ; 
         }
      }
      return res ;
  }    
}


template <typename T>
constexpr DenseMatrix<T> DenseMatrix<T>::minors(const std::size_t r , const std::size_t c ) const 
{
      //if(r > 0 && r <= row && c > 0 && c <= columns)
            
      if(r > 0 && r <= Rows && c > 0 && c <= Cols) // 1-based
      {
          DenseMatrix<T> res{Rows-1 , Cols-1} ;

          for(std::size_t ri = 1 ; ri <= ( Rows - (r >= Rows)); ri++ )
          {
             for(std::size_t ci = 1 ; ci <= (Cols - (c >= Cols)); ci++ )
             {
                  res(ri - (ri > r), ci - (ci > c)) = this->data[(ri-1) * Cols + (ci -1)];
             }
          }  
          return res;  
      }
      else
      {
          throw InvalidSizeException("Index for minor out of range");  
      }
       
}







/// @fun compute the determinant of square matrix  
//  with the power way (the minors way) just for small matrices 
//  
//  Method 
template <typename T>
constexpr auto DenseMatrix<T>::det() -> T  
{
   if(! this->isSquare())
   {
     throw InvalidSizeException("Matrix must be SQUARE for compute the DETERMINANT WITH MINORS-method");     
   }
   
   DenseMatrix<T> a{*this};   // matrix used to compute the step 

   T d = 0 ; // value of determinant   
   auto rows = a.size1();
   auto cols = a.size2();    
   
   if(rows==cols) // a is a square Matrix 
   {
      if(rows==1)  // 1x1 matrix
      {
         d = a(1,1);   
      }
      else if(rows==2)
      {
         d = a(1,1)*a(2,2) - a(2,1)*a(1,2) ;    
      }
      else
      {
            // matrix equal or larger than 3x3
            for(auto c=1 ; c <= cols ; c++)
            {
                DenseMatrix<T> M = a.minors(1,c);

                d += pow(-1,1+c) * a(1,c) * M.det();
            }
      
      }
    } 
    return d;
}



// Non Member Function
template <class U>
auto _det(const DenseMatrix<U>& a ) -> U
{
   U d = 0;   // value of determinant ;   
   std::size_t row  = a.size1();   
   std::size_t cols = a.size2();   
   
   if(row == cols)      
   {
      // square Matrix
      if(row == 1)
      {     
         d = a(1,1);    // 1x1 Matrix
      }
      else if(row == 2)
      {
         d = a(1,1) * a(2,2) - a(2,1) * a(1,2);
      }
      else
      {
         // this is a 3x3 matrix or larger    
         for(std::size_t c = 1 ; c <= cols ; c++ )
         {
            DenseMatrix<U> M = a.minors(1,c);
            
            d += pow(-1,1+c) * a(1,c) * _det(M);
         }
      }
   }   
   else
   {
      throw InvalidSizeException(">>Martrix must be square<<");
   }  
   return d;   
}










//-- nom member function 
//
template<typename U>
std::ostream& operator<<(std::ostream& os, const DenseMatrix<U>& m )
{
    for(auto i=1 ; i <= m.size1() ; i++ ){
      for(auto j=1 ; j <= m.size2() ; j++){
         std::cout << std::setw(8) << (fabs(m(i,j)) > 1.0e-14 ? m(i,j) : 0) << "  " ; 
      }
      std::cout << std::endl ;
    }  
}

       


template<typename U>
DenseMatrix<U> operator+(const DenseMatrix<U>& m1 , const DenseMatrix<U>& m2) 
{
      if(m1.size1() != m2.size1() || 
         m1.size2() != m2.size2()    )
      {
           throw InvalidSizeException(">>> Matrix dimension doesn't match in operator+  <<<"); 

      }
      else
      {
          DenseMatrix<U> res(m1.size1(), m1.size2());  
# pragma omp parallel for        
          for(auto i=1 ; i <= res.Rows ; i++){
              for(auto j=1; j<= res.Cols ; j++){
                 res(i,j) = m1(i,j) + m2(i,j) ; 
                 if(res(i,j)) res.nnz++; 
              }
          }    
         return res;   
      }
}



template<typename U>
DenseMatrix<U> operator-(const DenseMatrix<U>& m1 , const DenseMatrix<U>& m2) 
{
      if(m1.size1() != m2.size1() || 
         m1.size2() != m2.size2()    )
      {
           throw InvalidSizeException(">>> Matrix dimension doesn't match in operator-  <<<"); 

      }
      else
      {
          DenseMatrix<U> res(m1.size1(), m1.size2());  
# pragma omp parallel for         
          for(auto i=1 ; i <= res.Rows ; i++){
              for(auto j=1; j<= res.Cols ; j++){
                 res(i,j) = m1(i,j) - m2(i,j) ; 
                 if(res(i,j) != 0 ) res.nnz++ ; 
              }
          }    
         return res;   
      }
}



// MvP (Matrix Vector Product)
// ijk - alghorithm 
template<typename T>
std::vector<T> operator*(const DenseMatrix<T>& A, const std::vector<T> x) 
{
      

      if(A.size2() != x.size()) 
      { 
          std::string to = "x" ;
          std::string mess = "Error occured in operator* attempt to perfor productor between\n>> op1: "
                        + std::to_string(A.size1()) + to + std::to_string(A.size2()) +
                        " and op2: " + std::to_string(x.size()) + " <<";
          throw InvalidSizeException(mess.c_str());
      }
      else
      {
         std::vector<T> b(A.size1(),0);
   # pragma omp parallel for 
            for(auto i=0; i< A.size1(); i++)
            {     
                for(auto j=0 ; j < A.size2() ; j++ )
                {  
                        b.at(i) += A(i+1,j+1) * x.at(j);
                }
            }

        return b ;           
      }
}

//perform mat-mat product using  Strassen algorithm time-complexity = O(N^(2.807))
template<typename U> 
DenseMatrix<U> operator* (const DenseMatrix<U>& m1, const DenseMatrix<U>& m2) 
{
      
      if( m1.size2() != m2.size1() )
      {
         std::string to = "x" ;
         std::string mess = "Error occured in operator* attempt to perfor productor between op1: "
                        + std::to_string(m1.size1()) + to + std::to_string(m1.size2()) +
                        " and op2: " + std::to_string(m2.size1()) + to + std::to_string(m2.size2()) ;
         throw InvalidSizeException(mess.c_str());
      }
      else
      {
         DenseMatrix<U> res(m1.size1(), m2.size2() );       
 
#pragma omp parallel for

         for(auto i=1 ; i<= res.size1() ; i++){
             for(auto j=1 ; j <= res.size2(); j++){
                for(auto k=1 ; k <= m1.size2() ; k++){
                   res(i,j) += m1(i,k)*m2(k,j);
                  // if(res(i,j) != 0) res.nnz++;     
                }
             }   
         }         
         return res;
    }
}


template<typename U>
DenseMatrix<U> operator*(const DenseMatrix<U>& m, const U& rhs) noexcept 
{
   
   DenseMatrix<U> res(m.Rows, m.Cols );

# pragma omp parallel for 

   for(auto i=1 ; i <= m.Rows ; i++ ){
      for(auto j=1 ; j<= m.Cols ; j++){
            res(i,j) = m(i,j) * rhs ;
      }
   }
   return res;
}


template<typename U>
DenseMatrix<U> operator/(const DenseMatrix<U>& m, const U& rhs) noexcept 
{
   
   DenseMatrix<U> res(m.Rows, m.Cols );

# pragma omp parallel for 

   for(auto i=1 ; i <= m.Rows ; i++ ){
      for(auto j=1 ; j<= m.Cols ; j++){
            res(i,j) = m(i,j) / rhs ;
      }
   }
   return res;
}


//  perform < Mat * Mat > product using Strassen alghorithm (square Matrix) 
// 
//  uses only for large square matrices ! 
//
template <typename T>
void strassen(const DenseMatrix<T>& A ,  const DenseMatrix<T>& B ,
                    DenseMatrix<T>& C  , const std::size_t tam ) 
{
    
    if( !A.isSquare() || !B.isSquare() )
    {
       throw InvalidSizeException(">>> Matrix must be square ! <<< \nException thrown in strassen product"); 
    }    
    else if (A.size1() != A.size1() )
    {
      std::string to = "x" ;
      std::string mess = "Error occured in strassen function: attempt to perfor productor between op1: "
                        + std::to_string(A.size1()) + to + std::to_string(A.size2()) +
                        " and op2: " + std::to_string(B.size1()) + to + std::to_string(B.size2()) ;

      throw InvalidSizeException(mess.c_str());
    }
    else
    {
      //  DenseMatrix<T> c(m1.size1(), m2.size2());
        if ( tam <= A.leafSize )
        {
           C = A*B ;
        }      
        else
        { 
            std::size_t newTam = tam/2 ;   
          
            DenseMatrix<T> a11(newTam,newTam); 
            DenseMatrix<T> a12(newTam,newTam); 
            DenseMatrix<T> a21(newTam,newTam); 
            DenseMatrix<T> a22(newTam,newTam); 
            
            DenseMatrix<T> b11(newTam,newTam); 
            DenseMatrix<T> b12(newTam,newTam); 
            DenseMatrix<T> b21(newTam,newTam); 
            DenseMatrix<T> b22(newTam,newTam); 
           
            DenseMatrix<T> c11(newTam,newTam); 
            DenseMatrix<T> c12(newTam,newTam); 
            DenseMatrix<T> c21(newTam,newTam); 
            DenseMatrix<T> c22(newTam,newTam); 
            
            DenseMatrix<T> p1(newTam,newTam);   
            DenseMatrix<T> p2(newTam,newTam); 
            DenseMatrix<T> p3(newTam,newTam); 
            DenseMatrix<T> p4(newTam,newTam); 
            DenseMatrix<T> p5(newTam,newTam); 
            DenseMatrix<T> p6(newTam,newTam); 
            DenseMatrix<T> p7(newTam,newTam); 
               
            DenseMatrix<T> AA(newTam,newTam);   // matrix of a-result
            DenseMatrix<T> BB(newTam,newTam);   // matrix of b-result
              
            // dividing the matrices in 4 sub-matrices
            for(auto i=1 ; i <= newTam ; i++){
               for(auto j=1 ; j <= newTam ; j++ ){   
                  
                    a11(i,j)  = A( i          , j           );
                    a12(i,j)  = A( i          , j + newTam  );
                    a21(i,j)  = A( i + newTam , j           );
                    a22(i,j)  = A( i + newTam , j + newTam  );
                  
                    b11(i,j)  = B( i          , j           );
                    b12(i,j)  = B( i          , j + newTam  );
                    b21(i,j)  = B( i + newTam , j           );
                    b22(i,j)  = B( i + newTam , j + newTam  );
                   
               }   
            }
            
            sum(a11,a22,AA,newTam);
            sum(b11,b22,BB,newTam);      
            strassen(AA,BB,p1,newTam);// p1 = (a11+a22) *(b11+b22)
            
            sum(a21,a22,AA,newTam);
            strassen(AA,b11,p2,newTam ); // p2 = (a21+ a22) * b11

            subtract(b12,b22,BB,newTam);
            strassen(a11,BB,p3,newTam ); // p3 = a11 * (b12-b22)

            subtract(b21,b11,BB,newTam);
            strassen(a22,BB,p4,newTam);  // p4 = a22 * (b21-b11)
            
            sum(a11,a12,AA,newTam);
            strassen(AA,b22,p5,newTam);  // p5 = (a11+a12) * b22 

            subtract(a21,a11,AA,newTam);
            sum(b11,b12,BB,newTam);
            strassen(AA,BB,p6,newTam);   // p6 = (a21 - a11) * (b11+b12)
            
            subtract(a12,a22,AA,newTam);
            sum(b21,b22,BB,newTam);      
            strassen(AA,BB,p7,newTam);   // p7 = (a12- a22) * (b21 + b22)

            sum(p1,p4,AA,newTam);
            subtract(p7,p5,BB,newTam);
            sum(AA,BB,c11,newTam);  // c11 = p1 + p4 - p5 + p7
            
            sum(p3,p5,c12,newTam);  // c12 = p3+p5

            sum(p2,p4,c21,newTam); // c21 = p2+p4

            subtract(p1,p2,AA,newTam);
            sum(p3,p6,BB,newTam);
            sum(AA,BB,c22,newTam); // c22 = p1-p2+p3+p6       



         // Grouping the result

# pragma omp parallel for            
            for(auto i=1 ; i <= newTam ; i++){
                for(auto j=1 ; j <= newTam ; j++ ){
                  
                   C( i          , j          ) = c11(i,j);
                   C( i          , j + newTam ) = c12(i,j);
                   C( i + newTam , j          ) = c21(i,j);
                   C( i + newTam , j + newTam ) = c22(i,j);
                }
            }
        }
    }
}


template <typename U>
auto sum(const DenseMatrix<U>& A , const DenseMatrix<U>& B , DenseMatrix<U>& C, std::size_t tam ) noexcept
{  
# pragma omp parallel for 

      for(auto i=1 ; i <=tam ; i++){
         for(auto j=1 ; j <= tam ; j++){
            C(i,j) = A(i,j) + B(i,j) ;
         }
      }   
}

template <typename U>
auto subtract(const DenseMatrix<U>& A , const DenseMatrix<U>& B ,
                              DenseMatrix<U>& C, std::size_t tam ) noexcept
{

# pragma omp parallel for 
      for(auto i=1 ; i <=tam ; i++){
         for(auto j=1 ; j <= tam ; j++){
            
            C(i,j) = A(i,j) - B(i,j) ;

         }
      }   
}







  }//algebra 
 }//numeric
}//mg
# endif


