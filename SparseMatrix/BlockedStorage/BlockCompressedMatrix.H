# ifndef __BLOCK_COMPRESSED_MATRIX_H__
# define __BLOCK_COMPRESSED_MATRIX_H__

# include "../SparseMatrix.H"

namespace mg { 
               namespace numeric {
                                    namespace algebra {
  
 /*----------------------------------------------------
  *   Abstract class for Block Compressed 
  *   Sparse Matrix 
  *
  *   @ Marco Ghiani Dec. 2017 Glasgow
  ---------------------------------------------------*/

template <typename T, std::size_t S1, std::size_t S2>

class BlockCompressedMatrix :
                              public SparseMatrix<T>
{


   public: 


      virtual T& operator()(const std::size_t , const std::size_t) noexcept override = 0 ;

      virtual const T& operator()(const std::size_t , const std::size_t) const noexcept override = 0 ;
      
      virtual void print() const noexcept override = 0;
  
      virtual void printBlockMatrix() const noexcept ;


   protected:
        
     virtual std::size_t findBlockIndex(const std::size_t, const std::size_t ) const noexcept = 0 ;

     virtual T findValue(const std::size_t , const std::size_t ) const noexcept = 0;



};

template<typename T, std::size_t S1, std::size_t S2>
inline void BlockCompressedMatrix<T,S1,S2>::printBlockMatrix() const noexcept  
{
      for(auto i=0; i < this->denseRows / S1 ; i++)
      {
          for(auto j=0 ; j < this->denseCols /S2 ; j++ )
             std::cout <<  findBlockIndex(i,j) << ' ' ;
          std::cout << std::endl;  
      }

}




  }//algebra 
 }//numeric
}//mg
# endif
