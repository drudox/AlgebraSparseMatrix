# ifndef __SQUARE_BLOCK_COMPRESSED_STORAGE_H__
# define __SQUARE_BLOCK_COMPRESSED_STORAGE_H__


# include "../BlockCompressedMatrix.H"
 

namespace mg {
               namespace numeric {
                                    namespace algebra {

// forward declarations
template <typename T, std::size_t S >
class SqBCSmatrix ;



template <typename T, std::size_t S>
std::ostream& operator<<(std::ostream& os , const SqBCSmatrix<T,S>& m ) noexcept ;


template <typename T, std::size_t S>
std::vector<T> operator*(const SqBCSmatrix<T,S>& m, const std::vector<T>& x );

/*---------------------------------------------------------------------------------
 *
 *      Square Block Compressed Storage (sparse) Matrix 
 *          
 *      Store square block into single vector in top-bottom , left-right order  
 *
 *     @Marco Ghiani Dec 2017, Glasgow 
 *
 - -------------------------------------------------------------------------------*/



template <typename Type, std::size_t BS>
class SqBCSmatrix 
                   : public BlockCompressedMatrix<Type,BS,BS>
{


      template <typename T, std::size_t S>
      friend std::ostream& operator<<(std::ostream& os , const SqBCSmatrix<T,S>& m ) noexcept ;

      template <typename T, std::size_t S>
      friend std::vector<T> operator*(const SqBCSmatrix<T,S>& m, const std::vector<T>& x );


   
   public:

     constexpr SqBCSmatrix(std::initializer_list<std::vector<Type>> dense );  
     
     constexpr SqBCSmatrix(const std::string& );  

     virtual ~SqBCSmatrix() = default ; 


     auto constexpr print_block(const std::vector<std::vector<Type>>& dense,
                                  std::size_t i, std::size_t j) const noexcept ; 
     
     auto constexpr validate_block(const std::vector<std::vector<Type>>& dense,
                                  std::size_t i, std::size_t j) const noexcept ; 
 
     auto constexpr insert_block(const std::vector<std::vector<Type>>& dense,
                                                       std::size_t i, std::size_t j) noexcept ;
      
     auto constexpr printSqBCS() const noexcept ; 
  
     using BlockCompressedMatrix<Type,BS,BS>::printBlockMatrix ; 

     Type& operator()(const std::size_t , const std::size_t ) noexcept override final; 

     const Type& operator()(const std::size_t , const std::size_t )const  noexcept override final; 


     auto constexpr printBlock(std::size_t i) const noexcept ;
  
     void constexpr print() const noexcept override final; 
    
   
   private:
    
    std::size_t bn  ;
    std::size_t bBS ;
    std::size_t nnz ;
    using SparseMatrix<Type>::denseRows ;
    using SparseMatrix<Type>::denseCols ;
    
    using SparseMatrix<Type>::dummy ;

    std::vector<Type>    ba_ ; 
    std::vector<std::size_t>  an_ ;
    std::vector<std::size_t>  ai_ ;
    std::vector<std::size_t>  aj_ ;

      
    std::size_t index =0 ;

    std::size_t constexpr findBlockIndex(const std::size_t r, const std::size_t c) const noexcept override final;  
    
    auto constexpr recomposeMatrix() const noexcept ;
    
     
    Type constexpr findValue( const std::size_t , const std::size_t ) const noexcept override final;    
      
};

//---------------------------      IMPLEMENTATION      ------------------------------

// 
//
template <typename T, std::size_t BS>
constexpr SqBCSmatrix<T,BS>::SqBCSmatrix(std::initializer_list<std::vector<T>> dense_ )
{
      this->denseRows = dense_.size();   
      auto it         = *(dense_.begin());
      this->denseCols = it.size();
      
      if( (denseRows*denseCols) % BS != 0 )
      {
            throw InvalidSizeException("Error block size is not multiple of dense matrix size");
      }
      
     std::vector<std::vector<T>> dense(dense_);
     bBS = BS*BS ;  
     bn  = denseRows*denseCols/(BS*BS) ;
      
    ai_.resize(denseRows/BS +1);
    ai_[0] = 1;
      
    for(std::size_t i = 0; i < dense.size() / BS ; i++)
    {    
        auto rowCount =0;
        for(std::size_t j = 0; j < dense[i].size() / BS ; j++)
        {
            if(validate_block(dense,i,j))
            {     
                  aj_.push_back(j+1);
                  insert_block(dense, i, j);
                  rowCount ++ ;
            }      
            
        }
        ai_[i+1] = ai_[i] + rowCount ;
     }
     printSqBCS();
}

//-- read dense matrix from file 
//
template <typename T, std::size_t BS>
constexpr SqBCSmatrix<T,BS>::SqBCSmatrix(const std::string& fname)
{
    std::ifstream f(fname , std::ios::in);
    if(!f)
    {
       throw OpeningFileException("error opening file in constructor !");
    }
    else
    {
       std::vector<std::vector<T>> dense;
       std::string line, tmp;
       T elem = 0 ;
       std::vector<T> row;
       std::size_t i=0, j=0 ;     
       
       while(getline(f, line))
       {
          row.clear();  
          std::istringstream ss(line);
          if(i==0) 
          {
            while(ss >> elem)
            {
              row.push_back(elem);
              j++;
            }
          }
          else
          {
            while(ss >> elem) 
              row.push_back(elem);
          }
          dense.push_back(row);  
          i++;  
       }
       
       this->denseRows = i;
       this->denseCols = j;
       
       bBS = BS*BS ;  
       bn  = denseRows*denseCols/(BS*BS) ;
      
       ai_.resize(denseRows/BS +1);
       ai_[0] = 1;

      
       for(std::size_t i = 0; i < dense.size() / BS ; i++)
       {    
          auto rowCount =0;
          for(std::size_t j = 0; j < dense[i].size() / BS ; j++)
          {
              if(validate_block(dense,i,j))
              {     
                    aj_.push_back(j+1);
                    insert_block(dense, i, j);
                    rowCount ++ ;
              }      
            
          }
          ai_[i+1] = ai_[i] + rowCount ;
       }
 
    }
    printSqBCS();

}




template <typename T,std::size_t BS>
inline auto constexpr SqBCSmatrix<T,BS>::printBlock(std::size_t i) const noexcept
{  
   auto w = i-1 ;
   auto k = 0;   
   for(std::size_t i = 0 ; i < BS ; ++i) 
   {
      for(std::size_t j=0 ; j < BS ; ++j )
      {
          std::cout << std::setw(8) << ba_.at(an_.at(w)-1+k) << ' ';
          k++;        
      }            
   }
}





// --- print out all the blocks 
//
template <typename T,std::size_t BS>
inline auto constexpr SqBCSmatrix<T,BS>::print_block(const std::vector<std::vector<T>>& dense,
                                                       std::size_t i, std::size_t j) const noexcept
{   
   for(std::size_t m = i * BS ; m < BS * (i + 1); ++m) 
   {
      for(std::size_t n = j * BS ; n < BS * (j + 1); ++n)
                  std::cout << dense[m][n] << ' ';
      std::cout << '\n';
   }
}


//   --- check for non zero block 
//
template <typename T,std::size_t BS>
inline auto constexpr SqBCSmatrix<T,BS>::validate_block(const std::vector<std::vector<T>>& dense,
                                                       std::size_t i, std::size_t j) const noexcept
{   
   bool nonzero = false ;
   for(std::size_t m = i * BS ; m < BS * (i + 1); ++m)
   {
      for(std::size_t n = j * BS ; n < BS * (j + 1); ++n)
      {
            if(dense[m][n] != 0) nonzero = true;
      }
   }
   return nonzero ;
}


//-----  insert block into Ba vector 
//
template <typename T,std::size_t BS>
inline auto constexpr SqBCSmatrix<T,BS>::insert_block(const std::vector<std::vector<T>>& dense,
                                                       std::size_t i, std::size_t j) noexcept
{   
   //std::size_t value = index;   
   bool firstElem = true ;
   for(std::size_t m = i * BS ; m < BS * (i + 1); ++m)
   {
      for(std::size_t n = j * BS ; n < BS * (j + 1); ++n)
      {    
            if(firstElem)
            {
                  an_.push_back(index+1);
                  firstElem = false ;
            }
            ba_.push_back(dense[m][n]);
            index ++ ;
      }
   }
}   
 
 
template <typename T, std::size_t BS> 
std::size_t constexpr SqBCSmatrix<T,BS>::findBlockIndex(const std::size_t r, const std::size_t c) const noexcept 
{
      for(auto j= ai_.at(r) ; j < ai_.at(r+1) ; j++ )
      {   
         if( aj_.at(j-1) == c+1  )
         {
            return static_cast<std::size_t>(j) ;
         }
      }
}

// --- print the four vector of SqBCS format
//
template <typename T, std::size_t BS>
auto constexpr SqBCSmatrix<T,BS>::printSqBCS() const noexcept 
{ 

  std::cout << "ba_ :   " ;
  for(auto &x : ba_ ) 
      std::cout << x << ' ' ;
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

template <typename T, std::size_t BS> 
void constexpr SqBCSmatrix<T,BS>::print() const noexcept 
{      
    for(auto i=0 ; i < denseRows  ; i++){
       for(auto j = 0; j < denseCols ; j++){
               std::cout<< findValue(i, j) <<'\t';
       }
       std::cout << std::endl;
    }
}


template <typename T, std::size_t BS>
T constexpr SqBCSmatrix<T,BS>::findValue(const std::size_t i, const std::size_t j) const noexcept
{
    auto index = findBlockIndex(i/BS, j/BS);
    if(index != 0)
      return ba_.at(an_.at(index-1)-1 + j%BS + (i%BS)*BS);
    else
      return T(0);
}



template<typename T, std::size_t BS>
T& SqBCSmatrix<T,BS>::operator()(const std::size_t r, const std::size_t c) noexcept
{
      dummy = findValue(r-1,c-1);
      return dummy ;
}

template <typename T, std::size_t BS>
const T& SqBCSmatrix<T,BS>::operator()(const std::size_t r, const std::size_t c)const  noexcept  
{
      dummy = findValue(r-1,c-1);
      return dummy ;
}




template <typename T, std::size_t BS> 
auto constexpr SqBCSmatrix<T,BS>::recomposeMatrix() const noexcept
{

    std::vector<std::vector<T>> sparseMat(denseRows, std::vector<T>(denseCols, 0));
    auto BA_i = 0, AJ_i = 0;
    //for each BCSR row
    for(auto r = 0; r < denseRows/BS; r++){
        //for each Block in row
        for(auto nBlock = 0; nBlock < ai_.at(r+1)-ai_.at(r); nBlock++){  
            //for each subMatrix (Block)
            for(auto rBlock = 0; rBlock < BS; rBlock++){
                for(auto cBlock = 0; cBlock < BS; cBlock++){
                    //insert value
                    sparseMat.at(rBlock + r*BS).at(cBlock + (aj_.at(AJ_i)-1)*BS) = ba_.at(BA_i);
                    ++BA_i;
                }
            }
        ++AJ_i;
        }
    }
    return sparseMat;
}


//----------------------   non member functions 

template <typename T, std::size_t BS>
std::ostream& operator<<(std::ostream& os , const SqBCSmatrix<T,BS>& m ) noexcept 
{
    for(auto i=0 ; i < m.denseRows ; i++)
    {
       for(auto j = 0; j < m.denseCols ; j++)
       {
         os << m.findValue(i, j) <<'\t';
       }
       os << std::endl;
    }
    return os;  
}


template <typename T, std::size_t BS>
std::vector<T> operator*(const SqBCSmatrix<T,BS>& m, const std::vector<T>& x )
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
            auto brows = m.denseRows/BS ;  
            auto bnze  = m.an_.size()   ;

            auto z=0;

            for(auto b=0 ; b < brows ; b++)
            {     
               for(auto j= m.ai_.at(b) ; j <= m.ai_.at(b+1)-1; j++ )
               {      
                  for(auto k=0 ; k < BS ; k++ )
                  {
                     for(auto t=0 ; t < BS ; t++)
                     {
                         y.at(BS*b+k) += m.ba_.at(z) * x.at(BS*(m.aj_.at(j-1)-1)+t) ;          
                         z++ ;
                     }     
                  }
               }   
            }

      }
      return y;      

}

      
  }//algebra
 }//numeric
}//mg 
# endif
