# ifndef __MATRIX_H__
# define __MATRIX_H__


# include <iosfwd>
# include <initializer_list>
# include <iomanip>
# include <fstream>
# include <cassert>
# include <cmath>
# include <vector>
# include <string>
# include "MatrixException.H"
# include <sstream>
# include <set>
# include <algorithm>
# include <utility>
# include <memory>
 

namespace mg { namespace numeric { namespace algebra {

template <typename data_type>
class Matrix {
   
   public:
      
      
      virtual ~Matrix(){};
      
      virtual data_type& operator()(const std::size_t , const std::size_t) noexcept = 0; 
      
      virtual const data_type& operator()(const std::size_t , const std::size_t) const noexcept = 0 ;
      
      virtual void print()const noexcept =0;
      
      auto constexpr printMatrix(std::ostream&) noexcept {}
};

/*
template<typename T>
virtual mg::numeric::linearAlgebra::Matrix<T>::~Matrix()
                                          {}
*/
  } // linearAlgebra
 }  // numeric 
}   // mg 

# endif









