#ifndef __SPARSE_BLOCK_CR_STORAGE_H__
#define __SPARSE_BLOCK_CR_STORAGE_H__

# include "../BlockCompressedMatrix.H"

namespace mg {
               namespace numeric {
                                    namespace algebra {


# define __TESTING__

// forward declaration
template <typename T, std::size_t S>
class SBCRSmatrix ;



//- prototype of non member function 

template <typename T ,std::size_t S>
std::ostream& operator<<(std::ostream& os , const SBCRSmatrix<T,S>& m ) noexcept ;

template <typename T, std::size_t S>
std::vector<T> operator*(const SBCRSmatrix<T,S>& m , const std::vector<T>& );
 


/*------------------------------------------------------------------------------------
 *   
 *   Square Block Compressed Row storage Matrix Class
 *   
 *   This storing consist to storing only non zero element inside each square Block
 *
 *
 *
 *   @Marco Ghiani Dec 2017 , Glasgow 
 *  
 -----------------------------------------------------------------------------------*/


template <typename Type, std::size_t Size>
class SBCRSmatrix 
                   : public BlockCompressedMatrix<Type,Size,Size>
{
   

   public:

      template <typename T ,std::size_t S>
      friend std::ostream& operator<<(std::ostream& os , const SBCRSmatrix<T,S>& m ) noexcept ;

      template <typename T , std::size_t S>
      friend std::vector<T> operator*(const SBCRSmatrix<T,S>& m , const std::vector<T>& );
      

   public:
      
     constexpr SBCRSmatrix(std::initializer_list<std::vector<Type>> );
     
     constexpr SBCRSmatrix(const std::string& ) ;
 
     virtual ~SBCRSmatrix() = default ; 


     auto constexpr validate_block(const std::vector<std::vector<Type>>& dense,
                                                    std::size_t i, std::size_t j) const noexcept ;
     
     auto constexpr insert_valueBlock(const std::vector<std::vector<Type>>& dense,
                                                    std::size_t i, std::size_t j) noexcept ;
     
     auto constexpr printSBCRS() const noexcept ;  
     
     using BlockCompressedMatrix<Type,Size,Size>::printBlockMatrix; 

     void constexpr print() const noexcept override final;       


     Type& operator()(const std::size_t , const std::size_t) noexcept; 

     const Type& operator()(const std::size_t , const std::size_t) const noexcept ; 

//--
   private:

     //std::size_t denseRows ;
     //std::size_t denseCols ;
 
     using SparseMatrix<Type>::denseRows ;
     using SparseMatrix<Type>::denseCols ;

     using SparseMatrix<Type>::dummy     ; 

     struct Block {
        std::size_t i_ ;
        std::size_t j_ ;
        Type        v_ ; 
     };
    
      

     std::vector< Block >       ba_ ; // container (contain 2 index : position inside a block +  value ) 
     
     std::vector<std::size_t>   an_ ; // block ptr
     
     std::vector<std::size_t>   aj_ ; // bcol ind 

     std::vector<std::size_t>   ai_ ; // brow ind

     std::size_t index_ = 0;          // used for blocks indexing

     std::size_t constexpr findBlockIndex(const std::size_t r, const std::size_t c) const noexcept override final; 


//  
     Type constexpr findValue(const std::size_t r, const std::size_t c ) const noexcept override final;

};


/*
 * ------------------------------------   method definition    ----------------------------------------
*                                   to be puit into *.cpp file later 
*/

template <typename T,std::size_t S>
inline constexpr SBCRSmatrix<T,S>::SBCRSmatrix(std::initializer_list<std::vector<T>> dense_ )
{
    
    this->denseRows = dense_.size();
    this->denseCols = (*dense_.begin()).size();

    if(
        denseRows % S != 0 ||
        denseCols % S != 0    
      )
    {
        throw InvalidSizeException("Error in SB_CSR_Matrix constructor:\n" 
                                    "Matrix dimension and Block size doesn't match" );   
    }
     
    std::vector<std::vector<T>> dense(dense_);

    ai_.resize(denseRows/S + 1);
    ai_.at(0) = 1;

    for(auto i=0 ; i < dense.size() / S ; i++ )
    {
         auto rowCount = 0;   
         
         for(auto j=0 ; j < dense[i].size() / S ; j++)
         {
            if(validate_block(dense,i,j))
            {
                insert_valueBlock(dense,i,j);  
                aj_.push_back(j+1);  
                rowCount++;  
            }     

         }
         ai_.at(i+1) = ai_.at(i) + rowCount ;

    }
    an_.push_back(index_+1) ;

#ifdef __TESTING__    
    printSBCRS();   
#endif 
}

// construct matrix by file
//
template <typename T , std::size_t S>
inline constexpr SBCRSmatrix<T,S>::SBCRSmatrix(const std::string& fname) 
{
    std::ifstream f(fname, std::ios::in);
   
    if(!f)
    {
       std::string mess = "Error operning file " + fname + 
                          ":  EXCEPTION THROWN in SBCRSmatrix constructor! " ; 
       throw OpeningFileException(mess.c_str());     
    }
    else
    {
       std::vector<std::vector<T>> dense ;     
       std::string line ;
       T elem =0 ;
       std::vector<T> row ;
       
       auto i=0 , j=0 ;
       
       while(getline(f,line))
       {    
          std::istringstream ss(line);
          j=0 ;  
          row.clear();
          while(ss >> elem)
          {
            row.push_back(elem);
            j++;
          }
          dense.push_back(row);
          i++ ;  
       }
 
       this->denseRows = i;
       this->denseCols = j;

 //      std::cout << i << "x" << j << std::endl;

       ai_.resize(denseRows/S + 1);
       ai_.at(0) = 1;

       for(auto i=0 ; i < dense.size() / S ; i++ )
       {
           auto rowCount = 0;   
         
           for(auto j=0 ; j < dense[i].size() / S ; j++)
           {
              if(validate_block(dense,i,j))
              {
                  insert_valueBlock(dense,i,j);  
                  aj_.push_back(j+1);  
                  rowCount++;  
              }     

           }
           ai_.at(i+1) = ai_.at(i) + rowCount ;
        }
        an_.push_back(index_+1) ;
# ifdef __TESTING__    
       printSBCRS();   
# endif     
    }  

}


template <typename T,std::size_t S>
inline auto constexpr SBCRSmatrix<T,S>::validate_block(const std::vector<std::vector<T>>& dense,
                                                       std::size_t i, std::size_t j) const noexcept
{   
   bool nonzero = false ;
   for(std::size_t m = i * S ; m < S * (i + 1); ++m)
   {
      for(std::size_t n = j * S ; n < S * (j + 1); ++n)
      {
            if(dense[m][n] != 0) nonzero = true;
      }
   }
   return nonzero ;
}

template <typename T,std::size_t S>
inline auto constexpr SBCRSmatrix<T,S>::insert_valueBlock(const std::vector<std::vector<T>>& dense,
                                                       std::size_t i, std::size_t j) noexcept
{  
   //std::size_t value = index;   
   bool firstElem = true ;
   for(std::size_t m = i * S , o = 1; m < S * (i + 1); ++m ,++o)
   {
      Block b; 
      //mp.clear();
      for(std::size_t n = j * S , l=1 ; n < S * (j + 1); ++n, ++l)
      {    
         if(dense[m][n] != 0)
         {
            if(firstElem)
            {
                  an_.push_back(index_+1);
                  firstElem = false ;
            }
            index_ ++ ;
            b.i_ = o ;
            b.j_ = l ; 
            b.v_ = dense[m][n] ;  
            
            ba_.push_back(b);
         }   
      }
   }
}   

  
  
template <typename T, std::size_t S>
inline void constexpr SBCRSmatrix<T,S>::print() const noexcept
{
     for(auto i =1 ; i <= denseRows ; i++){
         for(auto j=1; j <= denseCols ; j++){
             
              std::cout << std::setw(8) << findValue(i,j) << ' ' ;    
         }
         std::cout << std::endl;
    }  

}


template <typename T , std::size_t S>
inline auto constexpr SBCRSmatrix<T,S>::printSBCRS()const noexcept 
{
  std::cout << "      (    " ;
  for(auto i=0 ; i < ba_.size() ; i++ )
  { 
      std::cout << ba_.at(i).v_ << ' ' ;
      
  }
  std::cout << std::endl; 
  
  std::cout << "ba_   (    " ;
  for(auto i=0 ; i < ba_.size() ; i++ )
  { 
      std::cout << ba_.at(i).i_ << ' ' ;
      
  }
  std::cout << std::endl; 
  
  std::cout << "      (    " ;
  for(auto i=0 ; i < ba_.size() ; i++ )
  { 
      std::cout << ba_.at(i).j_ << ' ' ;
      
  }
  std::cout << std::endl; 


  std::cout << "an_ :   " ;
  for(auto &x : an_ ) 
      std::cout <<  x << ' ' ;
    std::cout << std::endl; 
 
  std::cout << "aj_ :   " ;
  for(auto &x : aj_ ) 
      std::cout <<  x << ' ' ;
    std::cout << std::endl; 
   
   std::cout << "ai_ :   " ; 
   for(auto &x : ai_ ) 
      std::cout << x << ' ' ;
    std::cout << std::endl; 
  

}


// 
//
template <typename T, std::size_t S> 
inline std::size_t constexpr SBCRSmatrix<T,S>::findBlockIndex(const std::size_t r,
                                                       const std::size_t c) const noexcept 
{
      
      for(auto j= ai_.at(r) ; j < ai_.at(r+1) ; j++ )
      {   
         if( aj_.at(j-1) == c  )
         {
            return static_cast<std::size_t>(j) ;
         }
      }
}

//
template <typename T, std::size_t S>
inline T constexpr SBCRSmatrix<T,S>::findValue(const std::size_t r, const std::size_t c )const noexcept  {
      
    assert(r > 0 && r <= denseRows && c > 0 && c <= denseCols );
      
    for(auto i=0; i < denseRows/S ; ++i){
        for(auto j=ai_.at(i) ; j < ai_.at(i+1) ; ++j){ 
           for(auto k=an_.at(j-1) ; k <= an_.at(j)-1 ; ++k )
           {    
              if (ba_.at(k-1).i_+(i*2) == r && (ba_.at(k-1).j_+ (aj_.at(j-1)-1)*2) == c )
              {
                  return ba_.at(k-1).v_  ;
              }
                  
          } 
       }    
    }
    return 0;
}

//-------   Operator overloading


template <typename T, std::size_t S>
T& SBCRSmatrix<T,S>::operator()(const std::size_t r, const std::size_t c) noexcept
{
      assert(r > 0 && r <= denseRows && c > 0 && c <= denseCols);
      dummy = findValue(r,c);
      return dummy ;
}


template <typename T, std::size_t S>
const T& SBCRSmatrix<T,S>::operator()(const std::size_t r, const std::size_t c) const noexcept  
{
      assert(r > 0 && r <= denseRows && c > 0 && c <= denseCols);
      dummy = findValue(r,c);
      return dummy ;
}


//--------------------------- non member function 

template <typename T ,std::size_t S>
std::ostream& operator<<(std::ostream& os , const SBCRSmatrix<T,S>& m ) noexcept 
{
      for(auto i=1 ; i <= m.denseRows ; i++ ){ 
         for(auto j=1 ; j<= m.denseCols ; j++ ){  
            os << std::setw(8) << m(i,j) << ' ' ;
         }
         os << std::endl;   
      } 
     return os; 

}

//    perform matrix times vector product 
//
template <typename T, std::size_t S>
std::vector<T> operator*(const SBCRSmatrix<T,S>& m , const std::vector<T>& x )
{
    std::vector<T> y(x.size());
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
                       
        auto brows = m.denseRows/S ;
           
           
       for(auto b=0; b < brows ; b++){
          for(auto j = m.ai_.at(b) ; j <= m.ai_.at(b+1)-1 ; j++){
             for(auto z= m.an_.at(j-1) ; z <= m.an_.at(j)-1 ; z++ ){
               y.at(S*b+(m.ba_.at(z-1).i_ - 1)) += m.ba_.at(z-1).v_ * 
                                                   x.at(S*(m.aj_.at(j-1)-1) + m.ba_.at(z-1).j_-1);
             }
          }
       }
    }
    return y;
}
 
  }//agebra
 }//numeric
}//mg
# endif
