# ifndef __COMPRESSED_MATRIX_STORAGE_H__ 
# define __COMPRESSED_MATRIX_STORAGE_H__


# include "../SparseMatrix.H" 

namespace mg { namespace numeric { namespace  algebra {

template <typename T >
class CompressedMatrix :
                        public SparseMatrix<T>
{

   public:
      
      virtual ~CompressedMatrix() = default ;

      virtual T& operator()(const std::size_t , const std::size_t) noexcept override = 0;

      virtual const T& operator()(const std::size_t , const std::size_t) const noexcept override = 0;
    
      void constexpr printCompressed() const noexcept ;
      
      virtual void insertAt(const std::size_t, const std::size_t, const T) noexcept = 0 ;
            
      virtual void print() const noexcept override = 0;
      
      auto constexpr nnz_() const noexcept { return this->nnz ; }

   protected:
     
      virtual std::size_t findIndex(std::size_t , std::size_t ) const noexcept =0;  
      
      virtual T findValue(const std::size_t , const std::size_t ) const noexcept override =0 ;

      std::size_t dim ;
};



template <typename T>
inline void constexpr CompressedMatrix<T>::printCompressed() const noexcept 
{
   std::cout << "aa_ :   " ;   
   for(auto &x : SparseMatrix<T>::aa_)
      std::cout << x << ' ' ;
   std::cout << std::endl;   
   
   std::cout << "ia_ :   " ;   
   for(auto &x : SparseMatrix<T>::ia_)
      std::cout << x << ' ' ;
   std::cout << std::endl;   
   
   std::cout << "ja_ :   " ;   
   for(auto &x : SparseMatrix<T>::ja_)
      std::cout << x << ' ' ;
   std::cout << std::endl;   
}





  }//algebra    
 }//numeric
}//mg
# endif




