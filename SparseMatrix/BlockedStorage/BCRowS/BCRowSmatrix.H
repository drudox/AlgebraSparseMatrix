# ifndef __BLOCK_COMPRESSED_SINGLEROW_MATRIX_H__
# define __BLOCK_COMPRESSED_SINGLEROW_MATRIX_H__

# include "../BlockCompressedMatrix.H"

# define __TESTING__

namespace mg {
                namespace numeric {
                                    namespace algebra {

//
//forward Declaration
template <typename Type,std::size_t S>
class BCRowSmatrix ;


template <typename T,std::size_t S>
std::vector<T> operator*(const BCRowSmatrix<T,S>& A , const std::vector<T>& x ) noexcept ;

template <typename T, std::size_t S>
std::ostream& operator<<(std::ostream& os , const BCRowSmatrix<T,S>& m ) noexcept ;


/*-------------------------------------------------------------------------------------------*
 *
 *    Variable (row-size) Block (Size) Compressed Row Storage Class 
 *    
 *    this way of storage, store non zero "block" only along the row 
 *    therefore the block will be 1 x (NNZ - consecutives cols ) 
 *
 *    @Marco Ghiani  Dec. 2017 , Glasgow
 *
 --------------------------------------------------------------------------------------------*/


template <typename Type, std::size_t BS=2>
class BCRowSmatrix  : 
                                public BlockCompressedMatrix<Type,BS,BS> 
                        
{

      template <typename T,std::size_t S>
      friend std::vector<T> operator*(const BCRowSmatrix<T,S>& A , const std::vector<T>& x ) noexcept ;

      template <typename T,std::size_t S>
      friend std::ostream& operator<<(std::ostream& os , const BCRowSmatrix<T,S>& m ) noexcept ;


//--
//
   public:
      
      constexpr BCRowSmatrix(std::initializer_list<std::vector<Type>>) noexcept;
      
      constexpr BCRowSmatrix(const std::string& ); 

      auto constexpr printBCRS() const noexcept ;
   
      void constexpr print() const noexcept override final ;

      Type& operator()(const std::size_t , const std::size_t ) noexcept override final;
      
      const Type& operator()(const std::size_t , const std::size_t )const noexcept override final;

      using SparseMatrix<Type>::size1;
      
      using SparseMatrix<Type>::size2;

//--
//
   private:

     using SparseMatrix<Type>::denseRows ;
     using SparseMatrix<Type>::denseCols ;

     using SparseMatrix<Type>::aa_ ;
     using SparseMatrix<Type>::ja_ ;
     using SparseMatrix<Type>::ia_ ;
     
     std::vector<std::size_t>  nz_ ;
     
     using SparseMatrix<Type>::dummy ;

     Type constexpr findValue(const std::size_t r , const std::size_t c) const noexcept override final; 
      
     std::size_t constexpr findBlockIndex(const std::size_t , const std::size_t ) const noexcept override final; 

};


template <typename T, std::size_t S>
inline constexpr BCRowSmatrix<T,S>::BCRowSmatrix(std::initializer_list<std::vector<T>> rows) noexcept
{
   this->denseRows = rows.size();
   this->denseCols =(*rows.begin()).size() ;
 

   
   ia_.push_back(1);
   nz_.push_back(1);
   auto i=0 , blockCount=0, elemCount=0 ;   
   for(auto& row : rows )
   {
      blockCount = 0;
      for(auto j=0 ; j< row.size() ; j++)
      {
          if(row.at(j) != 0 )  
          {
             ja_.push_back(j+1);
             blockCount++ ;
             for(auto k=j ;k < row.size() && row.at(k) != 0 && j < row.size() ; k++ , ++j)
             {
                   aa_.push_back(row.at(k));
                   elemCount++ ;
             }
             nz_.push_back(elemCount+1);
          }
            
      }
      ia_.push_back(ia_.at(i)+blockCount);
      i++ ;
   }
   # ifdef __TESTING__
   printBCRS();   
   # endif
}

// --- construc from file 
//
template<typename T, std::size_t S>
inline constexpr BCRowSmatrix<T,S>::BCRowSmatrix(const std::string& fname) 
{
     std::ifstream f(fname , std::ios::in);
     
     if(!f)
     {
         throw OpeningFileException("EXCEPTION TROWN:\n >>> Error opening file in BCRowSmatrix class constructor <<<");   
     }       
     else
     {
         std::vector<std::vector<T>> rows;
         std::vector<T> r ;
         std::string line ;
         T elem = 0;

         auto i=0,  j=0;

         while(getline(f,line))
         {  
            r.clear();
            std::istringstream ss(line);      
            
            if(i==0)
            {
                while(ss >> elem){
                   r.push_back(elem);  
                   j++ ;  
                }   
            }
            else
            {
               while(ss >> elem)
                  r.push_back(elem);   
                                 
            }
            rows.push_back(r);
            i++;
         }
         this->denseRows = i;
         this->denseCols = j;
         
         i=0; j=0;
         auto elemCount =0;
         auto blockCount=0;
    
         ia_.push_back(1);
         nz_.push_back(1);
        
        for(auto& row : rows )
         {
            blockCount = 0;
            for(auto j=0 ; j< row.size() ; j++)
            {
                if(row.at(j) != 0 )  
                {
                    ja_.push_back(j+1);
                    blockCount++ ;
                    for(auto k=j ;k < row.size() && row.at(k) != 0 && j < row.size() ; k++ , ++j)
                    {
                         aa_.push_back(row.at(k));
                         elemCount++ ;
                    }
                    nz_.push_back(elemCount+1);
                }
            
            }
            ia_.push_back(ia_.at(i)+blockCount);
            i++ ;
         }
        # ifdef __TESTING__ 
         printBCRS();   
        # endif     
     } 
}


// print the storage 
//
template<typename T,std::size_t S>
inline auto constexpr BCRowSmatrix<T,S>::printBCRS() const noexcept 
{
   std::cout << "Af:     " ;   
   for(auto& x : aa_ )
       std::cout << x << ' ' ;
   std::cout << std::endl;
   
   std::cout << "Colind: " ;
   for(auto& x : ja_ )
       std::cout << x << ' ' ;
   std::cout << std::endl; 
   
   std::cout << "RowPtr: " ;
   for(auto& x : ia_ )
       std::cout << x << ' ' ;
   std::cout << std::endl; 
   
   std::cout << "Nzptr:  " ;
   for(auto& x : nz_ )
       std::cout << x << ' ' ;
   std::cout << std::endl; 
}


// find value into original sparse matrix
// 
template <typename T,std::size_t S>
inline T constexpr BCRowSmatrix<T,S>::findValue(const std::size_t row, const std::size_t col) const noexcept
{
    for(auto i = ia_.at(row-1)-1 ; i < ia_.at(row)-1 ; i++ )
    {  
      auto w = i ;
      auto z = 0 ;
      for(auto k= nz_.at(i)-1 ; k < nz_.at(i+1)-1 ; k++ ) 
      {
            
        if(ja_.at(w)+z  == col ) 
        {
          return aa_.at(k) ; 
        }
        z++ ;
      }
    }   
}

template<typename T, std::size_t S>
std::size_t constexpr BCRowSmatrix<T,S>::findBlockIndex(const std::size_t row,
                                                        const std::size_t col  ) const noexcept  
{
    for(auto i = ia_.at(row-1)-1 ; i < ia_.at(row)-1 ; i++ )
    {  
      auto w = i ;
      auto z = 0 ;
      for(std::size_t k= nz_.at(i)-1 ; k < nz_.at(i+1)-1 ; k++ ) 
      {
            
        if(ja_.at(w)+z  == col ) 
        {
          return k ; 
        }
        z++ ;
      }
    } 
}


template <typename T, std::size_t S>
void constexpr BCRowSmatrix<T,S>::print() const noexcept 
{
    
    for(std::size_t i=1 ; i <= size1() ; i++)
    {
       for(std::size_t j=1 ; j <= size2() ; j++)
       {
            std::cout << std::setw(8) << this->operator()(i,j) << ' ';   //findValue(i,j) << ' ';
       }
       std::cout << std::endl;
    }
   
}

//
// operator overload

template<typename T, std::size_t S>
T& BCRowSmatrix<T,S>::operator()(const std::size_t r, const std::size_t c) noexcept 
{
      assert(r > 0 && r <= size1() && c > 0 && c <= size2());

      dummy = findValue(r,c);  
      return dummy ;
}


template<typename T, std::size_t S>
const T& BCRowSmatrix<T,S>::operator()(const std::size_t r, const std::size_t c)const noexcept 
{
      assert(r > 0 && r <= size1() && c > 0 && c <= size2());

      dummy = findValue(r,c);
      return dummy ;
}




//------ non member function 

// SpMV -- perform matrix \times vector
//
template <typename T, std::size_t S>
std::vector<T> operator*(const BCRowSmatrix<T,S>& A , const std::vector<T>& x ) noexcept 
{
    
    if(A.size2() != x.size() )
    {
       std::cerr << "Matrix * vector error ! Dimension doesn't match!!" << std::endl;                  
       exit(-1);
    }

    std::vector<T> y(x.size());
    auto startCol = 0 , t = 0;  
    for(auto i=0 ; i < A.size2() ; i++ )
    {
       y.at(i) = 0.;    
       for(auto j= A.ia_.at(i)-1 ; j < A.ia_.at(i+1)-1 ; j++)
       {    
            startCol = A.ja_.at(j);
            t = 0;
            for(auto k = A.nz_.at(j)-1 ; k < A.nz_.at(j+1)-1 ; k++ )
            {   
                y.at(i) += A.aa_.at(k) * x.at(startCol-1 + t) ;
                ++t ;
            }
       }
    }  
    return y;  
}



//-- output 
//
template <typename T, std::size_t S>
std::ostream& operator<<(std::ostream& os , const BCRowSmatrix<T,S>& m ) noexcept 
{
      for(auto i=1; i <= m.size1() ; ++i)
      {
         for(auto j=1 ; j <= m.size2() ; ++j)
             os << std::setw(8) << m(i,j) << ' ' ;
         os << std::endl;  
      }
      return os;
 
}





  }//algebra 
 }//numeric
}//mg
# endif
