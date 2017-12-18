# ifndef __CRS_MATRIX_H__
# define __CRS_MATRIX_H__


# define  __DEBUG__
# define __TESTING__

# include "../CompressedMatrix.H"

namespace mg { namespace numeric { namespace algebra {

// foward declarations 
template <typename U>
class CRSmatrix ;


template <typename U>
std::ostream& operator<<(std::ostream& os , const CRSmatrix<U>& m ); 

template <typename U>
std::vector<U> operator*(const CRSmatrix<U>& , const std::vector<U>& x);

template<typename U>
CRSmatrix<U> operator*(const CRSmatrix<U>& m1, const CRSmatrix<U>& m2) ;



/*------------------------------------------------------------
 *
 *   Class for Compressed Row Storage (sparse) Matrix  
 *
 *    @Marco Ghiani Nov 2017 :
 *    v 1.01 testing 
 *
 ------------------------------------------------------------*/

template <typename Type>
class CRSmatrix :
                             public CompressedMatrix<Type>
{

         template <typename U>
         friend std::ostream& operator<<(std::ostream& os , const CRSmatrix<U>& m ); 

         template <typename U>
         friend std::vector<U> operator*(const CRSmatrix<U>& , const std::vector<U>& x);

         template<typename U>
         friend CRSmatrix<U> operator*(const CRSmatrix<U>& m1, const CRSmatrix<U>& m2) ;


      public:
         
         constexpr CRSmatrix(std::initializer_list<std::initializer_list<Type>> row ) noexcept;
         
         constexpr CRSmatrix(std::size_t i, std::size_t j) noexcept ;
         
         constexpr CRSmatrix(const std::string& );
         
         virtual  ~CRSmatrix() = default ;
         
         virtual Type& operator()(const std::size_t , const std::size_t) noexcept override final;

         virtual const Type& operator()(const std::size_t , const std::size_t) const noexcept override final ;
         void constexpr print() const noexcept override final;
      
         auto constexpr printCRS() const noexcept;
         
         using CompressedMatrix<Type>::printCompressed;

      
      private:
      
        using SparseMatrix<Type>::aa_ ;
        using SparseMatrix<Type>::ia_ ;
        using SparseMatrix<Type>::ja_ ;
      
        using SparseMatrix<Type>::denseRows ;
        using SparseMatrix<Type>::denseCols ;

        using SparseMatrix<Type>::nnz       ;
        using SparseMatrix<Type>::zero      ;     

        std::size_t constexpr findIndex(std::size_t row, std::size_t col) const noexcept override final ;
        
        Type constexpr findValue(const std::size_t , const std::size_t ) const noexcept override  final ;

        void insertAt(const std::size_t row, const std::size_t col,const Type val) noexcept override final;
 
 };

//---------------------------------     Implementation 

template<typename T>
inline constexpr CRSmatrix<T>::CRSmatrix(std::initializer_list<std::initializer_list<T>> row ) noexcept
{
    this->denseRows = row.size();
    this->denseCols =(*row.begin()).size() ;

    auto i=0, j=0, RowCount=0;  
    
    ia_.resize(denseRows+1);
    ia_[0] = 0;
    for(auto & r : row)
    {
        j=0 ;  RowCount=0 ;    
        for(auto & c : r)
        {
           if( c != 0.0 )
           {  
             aa_.push_back(c);
             ja_.push_back(j);
             RowCount++;
           }
           j++; 
        }
        i++;
        ia_[i] = ia_[i-1] + RowCount ;
    }
    nnz = aa_.size() ;
#ifdef __TESTING__
   printCompressed();   
# endif
}

//
//
template <typename T>
inline constexpr CRSmatrix<T>::CRSmatrix(std::size_t i, std::size_t j) noexcept 
                                                                         
{
      this->denseRows= i;
      this->denseCols= j; 
      
      aa_.resize(denseRows);
      ja_.resize(denseRows);
      ia_.resize(denseRows+1);
} 


// -- construct from file 
//
template <typename T>
constexpr CRSmatrix<T>::CRSmatrix(const std::string& filename )
{
      
      std::ifstream f( filename , std::ios::in );

      if(!f)
      {  
            std::string mess = "Error opening file  " + filename +
                               "\n Exception thrown in COOmatrix constructor" ;
            throw OpeningFileException(mess);
      }

      if( filename.find(".mtx") != std::string::npos )    // Coo format 
      {
          
          std::string line;
          T elem = 0 ;
          
          getline(f,line); // jump out the header line
          
          line = " " ;

          auto i=0; 
          auto i1=0 , j1=0 ;    
          while(getline(f,line))
          {
             std::istringstream ss(line);
             if(i==0)
             {
                 ss >> denseRows ; 
                 ss >> denseCols ;
                 ss >> nnz ;  
                
                ia_.resize(denseRows + 1); 
                ia_.at(0)=0;
             }
             else
             {
                ss >> i1 >> j1 >> elem ;

                for(auto i=(i1) ; i <= denseRows ; i++)  
                  ia_.at(i)++ ;
                
                ja_.insert(ja_.begin() + (ia_.at(i1)-1),  j1-1  ); 
                aa_.insert(aa_.begin() + (ia_.at(i1)-1),  elem  );  
             }                  
             i++ ; 

          }
          nnz = aa_.size();
      }
      else
      {
           std::string line;
           std::string tmp ;
           T elem = 0; 
           auto RowCount=0 ;           
           
           ia_.push_back(0) ;
           
           auto i=0 , j=0 , k = 0; 
           while(getline(f,line))
           {
              std::istringstream ss(line);
              j=0; RowCount=0;
              
              while (ss >> elem)
              {      
                  if(elem != 0)  
                  {
                     aa_.push_back(elem);
                     ja_.push_back(j);
                     RowCount++;
                  }
                  j++;
              this->denseCols = j;
              }
              ia_.push_back(ia_.at(i)+RowCount);
              i++ ;
           }
           this->denseRows = i;
      nnz = aa_.size() ; 
     }  

#  ifdef __TESTING__
     printCompressed(); 
#  endif      
}

// print out the CRS storage 
//
template <typename T>
inline auto constexpr CRSmatrix<T>::printCRS() const noexcept 
{
   std::cout << "aa_ :   " ;   
   for(auto &x : aa_)
      std::cout << x << ' ' ;
   std::cout << std::endl;   
   
   std::cout << "ja_ :   " ;   
   for(auto &x : ja_)
      std::cout << x << ' ' ;
   std::cout << std::endl;   
   
   std::cout << "ia_ :   " ;   
   for(auto &x : ia_)
      std::cout << x << ' ' ;
   std::cout << std::endl;   
}

// print out the whole matrix
//
template<typename T>
inline void constexpr CRSmatrix<T>::print() const noexcept
{     
      
      for(std::size_t i=1 ; i <= denseRows ; i++)
      {
         for(std::size_t j=1 ; j <= denseCols ; j++)
                  std::cout << std::setw(8) << this->operator()(i,j) << ' ' ;        
         std::cout << std::endl;
      }
}


//- - private utility function 
//
template<typename T>
inline std::size_t constexpr CRSmatrix<T>::findIndex(std::size_t row, std::size_t col) const noexcept
{
    
    assert( row >= 0 && row < denseRows 
         && col >= 0 && col < denseCols    );

    auto jit = std::find(ja_.begin()+ia_.at(row) , ja_.begin()+ia_.at(row+1),col );  
    
    return static_cast<std::size_t>(std::distance(ja_.begin(), jit) );

}

template<typename T>
T constexpr CRSmatrix<T>::findValue(const std::size_t row, const std::size_t col) const noexcept 
{
      assert( row > 0 && row <= denseRows 
           && col > 0 && col <= denseCols    );
 
      const auto j = findIndex(row-1,col-1);

      if(j < ia_.at(row))
            return aa_.at(j);
      else
            return zero ;
    
}

//
//
template <typename T>
inline void CRSmatrix<T>::insertAt(const std::size_t row, const std::size_t col,const T val) noexcept 
{
   if(val != 0)
   {
      auto j = findIndex(row,col);
      
      if( j < ia_.at(row+1) )
      {
         aa_.at(j) = val;    // change non zero   
      }
      else                  // remove non-zero 
      {
          for(auto i=(row+1) ; i <= this->denseRows ; i++ ) 
             ia_.at(i)++ ;

          ja_.insert(ja_.begin() + (ia_.at(row+1)-1) , col);
          aa_.insert(aa_.begin() + (ia_.at(row+1)-1) , val);
      }
   
   }
}


//--
template<typename T>
inline const T& CRSmatrix<T>::operator()(const std::size_t row, const std::size_t col) const noexcept 
{
      
    assert( row > 0 && row <= denseRows 
         && col > 0 && col <= denseCols );


      const auto j = findIndex(row-1,col-1);

      if(j < ia_.at(row))
            return aa_.at(j);
      else
            return zero ;
}

//--
template <typename T>
inline T& CRSmatrix<T>::operator()(const std::size_t row, const std::size_t col) noexcept  
{
    assert( row > 0 && row <= denseRows 
         && col > 0 && col <= denseCols );

      const auto j = findIndex(row-1,col-1);

      if(j < ia_.at(row))
            return aa_.at(j);
      else
            return zero ;
}


// ------ non member function 

template <typename T>
std::ostream& operator<<(std::ostream& os , const CRSmatrix<T>& m )
{
      for(auto i=1 ; i <= m.size1() ; i++ ){
            for(auto j=1 ; j <= m.size2() ; j++){
                std::cout << std::setw(8) << m.findValue(i,j) << " " ;  
            }
      std::cout << std::endl;       
      }
}


// ----- perform product

template <typename U>
std::vector<U> operator*(const CRSmatrix<U>& m, const std::vector<U>& x)
{
    if(m.size2() != x.size() )
    {
       std::string to = "x" ;
       std::string mess = "Error occured in operator* attempt to perfor productor between op1: "
                        + std::to_string(m.size1()) + to + std::to_string(m.size2()) +
                        " and op2: " + std::to_string(x.size());
       throw InvalidSizeException(mess.c_str());
    }
    std::vector<U> y(m.size1());
    for(auto i=0 ; i < m.size1() ; i++){
      for(auto j= m.ia_.at(i) ; j < m.ia_.at(i+1) ; j++){
          y.at(i) += m.aa_.at(j) * x.at(m.ja_.at(j));
      }
    }
    return y;
}


//--- Perform CRS * CRS 
//

template<typename T>
CRSmatrix<T> operator*(const CRSmatrix<T>& m1, const CRSmatrix<T>& m2) 
{

      if( m1.size2() != m2.size1() )
      {
         std::string to = "x" ;
         std::string mess = "Error occured in operator* attempt to perfor productor between op1: "
                        + std::to_string(m1.size1()) + to + std::to_string(m1.size2()) +
                        " and op2: " + std::to_string(m2.size1()) + to + std::to_string(m2.size2()) ;
         throw InvalidSizeException(mess.c_str());
      }
      
      CRSmatrix<T> res(m1.size1(), m2.size2());
      double sum ;
      for(std::size_t i=1 ; i <= res.size1() ; i++ )    
      {            
            for(std::size_t j=1; j<= res.size2() ; j++ )
            {     
                  sum = 0;
                  for(std::size_t k=1; k<= m1.size2() ; k++ )
                  {
                     sum += m1(i,k)*m2(k,j) ;
                  }
                  res.insertAt(i-1,j-1,sum);
            }
      }      
      return res;
      
}

/*
template<typename T>
CRSmatrix<T> operator*(const CRSmatrix<T>& m1, const CRSmatrix<T>& m2) 
{
      if( m1.size2() != m2.size1() )
      {
          throw std::runtime_error("Exception in dot product :  Matrix dimension must agree !"); 
      }
      
      std::vector<T> col ;
      CRSmatrix<T> res(m1.size1(), m2.size2());
      //ret.a.clear();
      col.resize(res.size2());

      for(std::size_t j=0 ; j < res.size2() ; j++ )   // loop over columns
      {      
            col.clear();      // 
            for(std::size_t i=0; i< res.size1() ; i++ )
                  col.push_back(static_cast<T> (m2(i+1,j+1) ) );
                         
            // perform multiplication
            col = m1*col;
            for(std::size_t i=0; i< col.size() ; i++ )
            {
                  res.insertAt(i,j,col.at(i));
            }
      }
      return res;
}
*/





  }//algebra
 }//numeric 
}// mg
# endif

