# ifndef ___COO_MATRIX_H___
# define ___COO_MATRIX_H___

# include "../../SparseMatrix.H"

# define __DEBUG__


namespace mg{ 
              namespace numeric { 
                                  namespace algebra {


// forward declararion
template <typename T> class COOmatrix ;

template <typename T>
std::ostream& operator<<(std::ostream& os , const COOmatrix<T>& m ) ; 

template <typename T>
std::vector<T> operator*(const COOmatrix<T>& , const std::vector<T>& ) ;



/*
 *    COO Matrix ( coordinate storage for Sparse Matrix )
 *    
 *    @ Marco Ghiani Nov 2017 Glasgow 
 *    
 */


template <typename data_type>
class COOmatrix :
                   public SparseMatrix<data_type>   
{                   


      template <typename T>
      friend std::ostream& operator<<(std::ostream& os , const COOmatrix<T>& m );

      template <typename T>
      friend std::vector<T> operator*(const COOmatrix<T>& , const std::vector<T>& ) ;



//--      
   public:   

      constexpr COOmatrix (std::initializer_list<std::initializer_list<data_type >> ) noexcept ;
      
      constexpr COOmatrix (const std::string& );
      
      virtual data_type& operator()(const std::size_t , const std::size_t) noexcept override ;

      virtual const data_type& operator()(const std::size_t , const std::size_t) const noexcept override;
      
      void constexpr print() const noexcept override ;
     
      auto constexpr printCOO() const noexcept ;      

//--
   private:

      using SparseMatrix<data_type>::aa_ ;
      using SparseMatrix<data_type>::ja_ ;
      using SparseMatrix<data_type>::ia_ ;

      data_type constexpr findValue(const std::size_t , const std::size_t ) const noexcept override;

      mutable data_type dummy ;
      using SparseMatrix<data_type>::zero  ;

};
 
/*
 *   ----------       Implementation    ---------
 *         -->   to be stored in a *.cpp file
 */
      
template <typename T>      
constexpr COOmatrix<T>::COOmatrix(std::initializer_list<std::initializer_list<T>> rows) noexcept 
{
      this->denseRows = rows.size();
      this->denseCols = (*rows.begin()).size();

      //std::cout << this->denseRows << " x " << this->denseCols << std::endl;      
      auto i=0 , j=0;
      for(auto& row : rows)
      {
        j=0;    
        for(auto& col : row )
        {
            if(col != 0)
            {     
               aa_.push_back(col);
               ja_.push_back(j);
               ia_.push_back(i);   
            }
            j++;
        }    
       i++;
      }
# ifdef __DEBUG__
      printCOO();
# endif

}


template <typename T>     
constexpr COOmatrix<T>::COOmatrix (const std::string& fname)
{
      
    std::ifstream f(fname , std::ios::in);  

      if(fname.find(".mtx") != std::string::npos )    // Coo format 
      {
          if(!f) 
          {
              std::string mess = "Error opening file  " + fname +
                                 "\n>>> Exception thrown in COOmatrix constructor <<<" ;
              throw OpeningFileException(mess);
          }
   
          std::string line;
          getline(f,line);
          std::size_t i1=0 , j1=0 ;
          T elem =0;
         // getline(f,line)
          auto i=0 ;
          while(getline(f,line))
          {
            std::istringstream ss(line);
            if(i==0)
            {
                  ss >> this->denseRows >> this->denseCols >> this->nnz ;
            }
            else
            {
                 ss >> i1 >> j1 >> elem ;
                 aa_.push_back(elem);
                 ia_.push_back(i1-1);
                 ja_.push_back(j1-1);
            }
            i++;      
          }  
      }     
      else
      {
          if(!f) 
          {
              std::string mess = "Error opening file  " + fname +
                                 "\n>>> Exception thrown in COOmatrix constructor <<<" ;
              throw OpeningFileException(mess);
          }
          else
          {
              std::string line ;    
              T elem = 0;    
              auto i=0, j=0 ;

              while(getline(f,line))
              {
                 std::istringstream ss(line);
                 j=0; 
                 while(ss >> elem)
                 {
                    if(elem != 0)    
                    {   
                        aa_.push_back(elem);
                        ja_.push_back(j);
                        ia_.push_back(i);
                    }   
                    j++ ;
                 }
                 ++i;
              }
              this->denseRows = i ;
              this->denseCols = j ;    
          }    
      }
#  ifdef __DEBUG__
      printCOO();
#  endif
}





template <typename T>
inline auto constexpr COOmatrix<T>::printCOO() const noexcept 
{
   std::cout << "aa_ :  " ;
   for(auto& x : aa_ )
       std::cout << x << ' ' ;
   std::cout << std::endl;   
 
   std::cout << "ia_ :  " ;
   for(auto& x : ia_ )
       std::cout << x << ' ' ;
   std::cout << std::endl;   
  
   std::cout << "ja_ :  " ;
   for(auto& x : ja_ )
       std::cout << x << ' ' ;
   std::cout << std::endl;   
}


template <typename T>
void constexpr COOmatrix<T>::print() const noexcept {
      
      for(auto i=1; i<= this->denseRows ; i++){
         for(auto j=1; j <= this->denseCols ; j++){
            std::cout << std::setw(8) << this->operator()(i,j) << ' ';
          //  findValue(i,j);
         }
         std::cout << std::endl;   
      }   
}



template <typename T>
T constexpr COOmatrix<T>::findValue(const std::size_t r, const std::size_t c) const noexcept
{
   
   auto i1 = ia_.begin();
   while( (i1 = std::find(i1, ia_.end(), r) )  != ia_.end() )   
   {     
           auto i = static_cast<std::size_t >( std::distance(ia_.begin() , i1) );
             if( ja_.at(i) == c )
                 return aa_.at(i)  ;
      ++i1;
    }    
    return T(0);   
    
}

//-----   operator Overloading 

template <typename T>
T& COOmatrix<T>::operator()(const std::size_t i, const std::size_t j) noexcept
{
      assert(i > 0 && i <= this->denseRows &&
             j > 0 && j <= this->denseCols    );

      dummy = findValue(i-1,j-1);
      return dummy ;
}

template <typename T>
const T& COOmatrix<T>::operator()(const std::size_t i, const std::size_t j) const noexcept
{
      assert(i > 0 && i <= this->denseRows &&
             j > 0 && j <= this->denseCols    );

      dummy =findValue(i-1,j-1);
      return dummy ;
}


/*
 *    Non member function 
 */ 


template <typename T>
std::ostream& operator<<(std::ostream& os, const COOmatrix<T>& m)  
{
      for(auto i=1; i <= m.denseRows ; i++){
            for(auto j=1 ; j <= m.denseCols ; j++){
                  os << std::setw(8) << m(i,j) << ' ';
            }      
            os << std::endl;
      }
}


template <typename T>
std::vector<T> operator*(const COOmatrix<T>& m, const std::vector<T>& x)
{
      
      std::vector<T> y(x.size(),0);
      
      if(m.size1() != x.size())
      {
          std::string to = "x" ;
          std::string mess = "Error occured in operator* attempt to perfor productor between op1: "
                        + std::to_string(m.size1()) + to + std::to_string(m.size2()) +
                        " and op2: " + std::to_string(x.size());
          throw InvalidSizeException(mess.c_str());
      }
      else
      {
        for(auto i=0 ; i < m.aa_.size() ; i++ )      
            y.at(m.ia_.at(i)) += m.aa_.at(i)* x.at(m.ja_.at(i));   
      }
      return y;

}



    } // algebra  
  } // namespace numeric
} // namespace mg
# endif
