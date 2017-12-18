# ifndef __MODIFIED_COMPRESSED_MATRIX_H__
# define __MODIFIED_COMPRESSED_MATRIX_H__


# include "../SparseMatrix.H"

namespace mg { namespace numeric { namespace algebra {


template <typename T>
class ModifiedCompressedMatrix : 
                                    public SparseMatrix<T>
{                                    

   public:

      virtual ~ModifiedCompressedMatrix() = default ;
      
      virtual T& operator()(const std::size_t , const std::size_t ) noexcept override = 0; 
      
      virtual const T& operator()(const std::size_t , const std::size_t) const noexcept override = 0;

      void constexpr printModCompressed() const noexcept ;
      
      virtual void print() const noexcept override = 0;

      auto constexpr size() const noexcept {return dim;}
      
      auto constexpr nnz() const noexcept { return this->aa_.size(); }

      std::vector<T> constexpr  diag() noexcept ;

   protected:

      virtual std::size_t findIndex(std::size_t , std::size_t ) const noexcept = 0;
      
      virtual T findValue(const std::size_t , const std::size_t ) const noexcept override = 0 ;


      std::size_t dim ;
  
};
  
template <typename T>
void constexpr ModifiedCompressedMatrix<T>::printModCompressed() const noexcept 
{
      std::cout << "AA :   " ;     
      for(auto& x : this->aa_ )
          std::cout << x << ' ' ;
      std::cout << std::endl;

      std::cout << "JA :   " ;     
      for(auto& x : this->ja_ )
          std::cout << x << ' ' ;
      std::cout << std::endl;
}

template <typename T>
std::vector<T> constexpr ModifiedCompressedMatrix<T>::diag() noexcept 
{
    std::vector<T> d(dim);  
    for(auto i=0; i <dim ; i++)
      d.at(i) = this->aa_.at(i) ;
      
    return d;
}



  }//algebra
 }//numeric
}//mg 
#endif
