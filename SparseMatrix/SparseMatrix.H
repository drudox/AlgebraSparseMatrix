# ifndef ___SPARSE_MATRIX_H___
# define ___SPARSE_MATRIX_H___

# include "../Matrix.H"

namespace mg{ namespace numeric { namespace algebra {

template <typename T> 
class SparseMatrix : 
                     public Matrix<T> 
{
     public: 
      
      virtual ~SparseMatrix() =default ;

      virtual T& operator()(const std::size_t , const std::size_t) noexcept override = 0 ;

      virtual const T& operator()(const std::size_t , const std::size_t) const noexcept override = 0 ;
      
      virtual void print() const noexcept override = 0;
      
      auto constexpr size1() const noexcept { return denseRows ;}

      auto constexpr size2() const noexcept { return denseCols ;}

      
    protected:
      
      std::vector<T> aa_ ;            // vectror of non zero elem 
      std::vector<std::size_t> ia_ ;  // vector row index / pointer 
      std::vector<std::size_t> ja_ ;  // vector col index / pointer 
      

      std::size_t denseRows ;
      
      std::size_t denseCols ;  
      
      std::size_t nnz ;       // non zero member

      T zero  = static_cast<T>(0.) ;
      
      virtual T findValue(const std::size_t , const std::size_t ) const noexcept = 0;
      
      mutable T dummy ;
};

//virtual void print() const noexcept override ;



  } // algebra
 } // numeric
} // mg

# endif
