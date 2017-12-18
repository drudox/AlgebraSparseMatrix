# ifndef __CCS_MATRIX_H__
# define __CCS_MATRIX_H__


# include "../CompressedMatrix.H"

# define __TESTING__

namespace mg {
               namespace numeric {
                                    namespace algebra {

// forward declaration 
template <typename U>
class CCSmatrix ;


template <typename U>
std::ostream& operator<<(std::ostream& out ,const CCSmatrix<U>& m );

template <typename U>
std::vector<U> operator*(const CCSmatrix<U>& m, const std::vector<U>& v);
 
template <typename U> 
CCSmatrix<U> operator*(const CCSmatrix<U>& m1, const CCSmatrix<U>& m2); 


 
// - CCSmatrix Class 
//
template < typename Type>
class CCSmatrix :                                       
                   public CompressedMatrix<Type>   //  base Class     
{

    public:  
    
      template <typename U>
      friend std::ostream& operator<<(std::ostream& out ,const CCSmatrix<U>& m );
      
      template <typename U>
      friend std::vector<U> operator*(const CCSmatrix<U>& m, const std::vector<U>& v);
 
      template <typename U> 
      friend CCSmatrix<U> operator*(const CCSmatrix<U>& m1, const CCSmatrix<U>& m2);



   public:
     
     constexpr CCSmatrix( std::initializer_list<std::vector<Type>> il) noexcept;
     
     constexpr CCSmatrix(std::size_t, std::size_t) noexcept ; 
     
     constexpr CCSmatrix(const std::string& );

     virtual  ~CCSmatrix() = default ; 

     virtual Type& operator()(const std::size_t i,const std::size_t j) noexcept override final;
     
     virtual const Type& operator()(const std::size_t i,const std::size_t j)const noexcept override final;
     
     void constexpr print()const noexcept override final; 

     auto constexpr printCCS()const noexcept;

     using CompressedMatrix<Type>::printCompressed;
      
   private:
   
     
        using SparseMatrix<Type>::aa_ ;
        using SparseMatrix<Type>::ia_ ;
        using SparseMatrix<Type>::ja_ ;
      
        using SparseMatrix<Type>::denseRows ;
        using SparseMatrix<Type>::denseCols ;

        using SparseMatrix<Type>::nnz       ;
        using SparseMatrix<Type>::zero      ;     

        std::size_t constexpr findIndex(std::size_t row, std::size_t col) const noexcept override final;
     
        void insertAt(const std::size_t, const std::size_t, const Type) noexcept override final;
        
        Type constexpr findValue(const std::size_t i, const std::size_t j) const noexcept override final;
};

/***************************************************************************************************/
//                                  Implementation                                                 //


template< typename T> 
constexpr CCSmatrix<T>::CCSmatrix( std::initializer_list<std::vector<T>> row) noexcept
{
     this->denseRows = row.size();
     this->denseCols = (*row.begin()).size();
     
      std::size_t i=1 , j=0 , colCount=0; 
      ja_.resize(denseCols+1);
      ja_[0] = 0;
      std::vector<std::vector<T>> temp( row );
      //w =0 ;
                
        for(size_t c = 0; c < temp[0].size(); c++)
        {
             colCount = 0;
             for(size_t r = 0; r < temp.size(); r++)
             {                
                 if(temp[r][c] != 0)
                 {
                   aa_.push_back( temp[r][c] ) ;
                   ia_.push_back(r); 
                   colCount++ ; 
                 } 
            }
            ja_[i] = ja_[i-1] + colCount ;
            i++ ;
     } 
     #ifdef __TESTING__
     printCompressed(); 
     # endif
}     

//
// 
template<typename T> 
constexpr CCSmatrix<T>::CCSmatrix(std::size_t row, std::size_t col) noexcept 
{
    this->denseRows = row ;
    this->denseCols = col ;

    aa_.resize(denseCols);
    ia_.resize(denseCols);
    ja_.resize(denseCols+1);  
}



//  construct from file - matrix data 
// 
template <typename T>
constexpr CCSmatrix<T>::CCSmatrix(const std::string& filename ) 
{
    std::ifstream f(filename , std::ios::in);

    if(!f)
    {
        std::string mess = "Error opening file  " + filename +
                               "\n Exception thrown in COOmatrix constructor" ;
        throw OpeningFileException(mess);
    }
    
    if( filename.find(".mtx") != std::string::npos )
    {
        std::string line;
        T elem =0; 
        
        getline(f,line); // jump out the header
        
        line = " ";
        
        auto i=0;
        auto i1 = 0, j1 = 0;
        
        while(getline(f,line))
        {
           std::istringstream ss(line);

           if(i==0)
           {
               ss >> denseRows ;
               ss >> denseCols ;
               ss >> nnz ;

               ja_.resize(denseCols+1);
               ja_.at(0) = 0; 
           }
           else
           {
               ss >> i1 >> j1 >> elem ;
               for(auto j=(j1); j<= denseCols ; j++ )
                   ja_.at(j)++ ;

               ia_.insert( ia_.begin() + (ja_.at(j1)-1) , i1-1 );
               aa_.insert( aa_.begin() + (ja_.at(j1)-1) , elem );
           }
           i++ ; 
        }
    }
    else
    {
        std::string line , tmp;
        std::vector<std::vector<std::string>> temp ;
        while( getline(f,line) )  
        {
           std::istringstream ss(line);
           std::vector<std::string> vtmp;
           vtmp.clear();
           while(ss >> tmp )
           {     
             vtmp.push_back (tmp);
           }   
           temp.push_back(vtmp) ;
         }

         this->denseRows = temp.size(); 
         this->denseCols = temp[0].size(); 
   
         std::size_t i=1 , j=0 , w=0; 
         ja_.resize(denseCols+1);
         ja_[0] = 0;
         for(size_t c = 0; c < temp[0].size(); c++)
         {     
            w =0 ;
            for(size_t r = 0; r < temp.size(); r++)
            {
                 if( std::stof(temp[r][c]) != 0)
                 {
                   //T val =      
                   aa_.push_back( std::stof(temp[r][c]) ) ;
                   ia_.push_back(r); 
                   w++ ; 
                 } 
            }
            ja_[i] = ja_[i-1] + w ;
            i++ ;
         }
    } 
    #ifdef __TESTING__ 
     printCompressed(); 
    #endif  
}     


// --- utility function
//
template<typename T>
inline std::size_t constexpr CCSmatrix<T>::findIndex(const std::size_t row, const std::size_t col) const noexcept
{
    auto ijt = std::find(ia_.begin()+ja_.at(col) , ia_.begin()+ja_.at(col+1), row );  
    
    return static_cast<std::size_t>(std::distance(ia_.begin(), ijt ) );

}

template<typename T>
T constexpr CCSmatrix<T>::findValue(const std::size_t row, const std::size_t col) const noexcept {
    const auto i = findIndex(row,col);  
    
    if(i < ja_.at(col+1))
      return aa_.at(i) ;
    else 
      return zero ;
}

//
// 
template<typename T>
inline const T& CCSmatrix<T>::operator()(const std::size_t row,const std::size_t col)const noexcept  
{
    const auto i = findIndex(row,col);  
    
    if(i < ja_.at(col+1))
      return aa_.at(i) ;
    else 
      return zero ;
}

//
//
template<typename T>
inline T& CCSmatrix<T>::operator()(const std::size_t row, const std::size_t col) noexcept  
{
    const auto i = findIndex(row,col);  
    
    if(i < ja_.at(col+1))
      return aa_.at(i) ;
    else 
      return zero ;
}


// insert element in the matrix at i,j 
//
template <typename T>
inline void CCSmatrix<T>::insertAt(const std::size_t row, const std::size_t col, const T val) noexcept 
{
   if(val != 0)
   {
      auto i = findIndex(row,col);
      
      if( i < ja_.at(col+1) )
      {
         aa_.at(i) = val;    // substitute if already exist in this position   
      }
      else                   // update index and pointer vector , insert value
      {
          for(auto j=(col+1) ; j <= this->denseCols ; j++ ) 
             ja_.at(j)++ ;

          ia_.insert(ia_.begin() + (ja_.at(col+1)-1) , row);
          aa_.insert(aa_.begin() + (ja_.at(col+1)-1) , val);
      }
   
   }
}

// print the whole matrix
//
template<typename T>
inline void constexpr CCSmatrix<T>::print()const noexcept 
{
      for(std::size_t i=0 ; i < denseRows ; i++)
      {
         for(std::size_t j=0 ; j < denseCols ; j++)
                  std::cout << std::setw(8) << this->operator()(i,j) << ' ' ;        
         std::cout << std::endl;
      }
}



//
template<typename T >
auto constexpr CCSmatrix<T>::printCCS()const noexcept
{    

    std::cout << "aa_ :   " ;   
    for(auto& i : aa_ )
      std::cout << i << ' ' ;
    std::cout << std::endl ;  
    

    std::cout << "ia_ :   " ;   
    for(auto& i : ia_ )
      std::cout << i << ' ' ;
    std::cout << std::endl ;  
    
    
    std::cout << "ja_ :   " ;   
    for(auto& i : ja_ )
      std::cout << i << ' ' ;
    std::cout << std::endl ;  

}

//   -----    non member function ----
//

template <typename T>
std::ostream& operator<<(std::ostream& os ,const CCSmatrix<T>& m )
{
    for(auto i=0; i < m.denseRows ; i++)
    {
       for(auto j=0 ; j < m.denseCols ; j++)
       {
          os << std::setw(8) << m.findValue(i,j) <<  ' ' ;
       }
       os << std::endl;
    }
      
    return os ;
}



// -- perform matrix times vector 
//
//
template <typename T>
std::vector<T> operator*(const CCSmatrix<T>& m, const std::vector<T>& x)
{
      
      if( m.size2() != x.size() )
      {
       std::string to = "x" ;
       std::string mess = "Error occured in operator* attempt to perfor productor between op1: "
                        + std::to_string(m.size1()) + to + std::to_string(m.size2()) +
                        " and op2: " + std::to_string(x.size());
       throw InvalidSizeException(mess.c_str());
      }
      
      std::vector<T> y(m.size1()) ;
      
      for(std::size_t i=0 ; i < x.size() ; i++)
      {
            for(std::size_t j = m.ja_.at(i) ; j <= m.ja_.at(i+1)-1 ; j++  )
            {     
                  y.at(m.ia_.at(j)) += m.aa_.at(j) * x.at(i);
            }
      }
      
      return y;
}




template<typename T>
CCSmatrix<T> operator*(const CCSmatrix<T>& m1, const CCSmatrix<T>& m2) 
{

      if( m1.size2() != m2.size1() )
      {
         std::string to = "x" ;
         std::string mess = "Error occured in operator* attempt to perfor productor between op1: "
                        + std::to_string(m1.size1()) + to + std::to_string(m1.size2()) +
                        " and op2: " + std::to_string(m2.size1()) + to + std::to_string(m2.size2()) ;
         throw InvalidSizeException(mess.c_str());
      }
      
      CCSmatrix<T> res(m1.size1(), m2.size2());
      double sum ;
      for(std::size_t i=0 ; i < res.size1() ; i++ )    
      {            
            for(std::size_t j=0; j< res.size2() ; j++ )
            {     
                  sum = 0;
                  for(std::size_t k=0; k< m1.size2() ; k++ )
                  {
                     sum += m1(i,k)*m2(k,j) ;
                  }
                  res.insertAt(i,j,sum);
            }
      }      
      return res;
      
}
//



  }//algebra    
 }//numeric
}// mg
# endif
