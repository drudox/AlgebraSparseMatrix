# ifndef __BLOCK_COMPRESSED_ROW_STORAGE_H__
# define __BLOCK_COMPRESSED_ROW_STORAGE_H__

# define __TESTING__

# include "../BlockCompressedMatrix.H"

namespace mg { 
                namespace numeric {
                                    namespace algebra {


// forward declarations
template <typename T, std::size_t R, std::size_t C>
class BCRSmatrix ;


template <typename T, std::size_t R, std::size_t C>
std::ostream& operator<<(std::ostream& os , const BCRSmatrix<T,R,C>& m );


template <typename T, std::size_t Br, std::size_t Bc >
std::vector<T> operator*(const BCRSmatrix<T,Br,Bc>& m, const std::vector<T>& x );


/*-------------------------------------------------------------------------------------------
 *
 *    ( Rectangular ) Blocks Compressed Rows Storage   
 *    
 *    Stored rectangular non-zero block into vector using :
 *
 *    - aa : vector containg blocks (top-bottom, left-right order)
 *    - an : vector of index to the first element of the block in the vector aa 
 *    - ia : (pointers row) stored the first block index of each row in the blocked matrix.
 *    - ja : stored the index of the block columns in the blocked matrix.
 *
 *    @Marco Ghiani Dec. 2017 
 *
 *
 --------------------------------------------------------------------------------------------*/


template < typename    Type,
           std::size_t BR  ,
           std::size_t BC   
         >
class BCRSmatrix 
                  :     public BlockCompressedMatrix<Type,BR,BC>
{


      template <typename T, std::size_t R, std::size_t C>
      friend std::ostream& operator<<(std::ostream& os , const BCRSmatrix<T,R,C>& m );

      template <typename T, std::size_t Br,std::size_t Bc>
      friend std::vector<T> operator*(const BCRSmatrix<T,Br,Bc>& m, const std::vector<T>& x );
 
//-
//
   public:

     constexpr BCRSmatrix(std::initializer_list<std::vector<Type>> dense );  
     
     constexpr BCRSmatrix(const std::string& );  

     virtual ~BCRSmatrix() = default ; 

     auto constexpr print_block(const std::vector<std::vector<Type>>& dense,
                                  std::size_t i, std::size_t j) const noexcept ; 
     
     auto constexpr validate_block(const std::vector<std::vector<Type>>& dense,
                                  std::size_t i, std::size_t j) const noexcept ; 
 
     auto constexpr insert_block(const std::vector<std::vector<Type>>& dense,
                                                  std::size_t i, std::size_t j) noexcept ;
      
     auto constexpr printBCRS() const noexcept ; 
  
     auto constexpr printBlock(std::size_t i) const noexcept ;
  
     void constexpr print() const noexcept override final; 
     
     Type& operator()(const std::size_t , const std::size_t) noexcept override final ;
     
     const Type& operator()(const std::size_t , const std::size_t) const noexcept override final ;
   

   
   private:
    
    std::size_t bn  ;
    std::size_t bBR ;
    std::size_t nnz ;
    
    using SparseMatrix<Type>::denseRows ;
    using SparseMatrix<Type>::denseCols ;

    
    using SparseMatrix<Type>::aa_ ;
    std::vector<std::size_t>  an_ ;

    using SparseMatrix<Type>::ia_ ;
    using SparseMatrix<Type>::ja_ ;

    using SparseMatrix<Type>::dummy;   // mutable Type var.  

    std::size_t index =0 ;

    std::size_t constexpr findBlockIndex(const std::size_t r, const std::size_t c) const noexcept override final;  
    
    auto constexpr recomposeMatrix() const noexcept ;
    
    Type constexpr findValue(const std::size_t i, const std::size_t j) const noexcept override final;
  

};

//---------------------------      IMPLEMENTATION      ------------------------------

// 
//
template <typename T, std::size_t BR, std::size_t BC>
constexpr BCRSmatrix<T,BR,BC>::BCRSmatrix(std::initializer_list<std::vector<T>> dense_ )
{
      this->denseRows = dense_.size();   
      auto it         = *(dense_.begin());
      this->denseCols = it.size();
      
      if( (denseRows*denseCols) % BR != 0 )
      {
            throw InvalidSizeException("Error block size is not multiple of dense matrix size");
      }
      
     std::vector<std::vector<T>> dense(dense_);
     bBR = BR*BC ;  
     bn  = denseRows*denseCols/(BR*BC) ;
      
    ia_.resize(denseRows/BR +1);
    ia_[0] = 1;
      
    for(std::size_t i = 0; i < dense.size() / BR ; i++)
    {    
        auto rowCount =0;
        for(std::size_t j = 0; j < dense[i].size() / BC ; j++)
        {
            if(validate_block(dense,i,j))
            {     
                  ja_.push_back(j+1);
                  insert_block(dense, i, j);
                  rowCount ++ ;
            }      
            
        }
        ia_[i+1] = ia_[i] + rowCount ;
     }
     # ifdef __TESTING__
     printBCRS();
     # endif
}

//-- read dense matrix from file 
//
template <typename T, std::size_t BR, std::size_t BC>
constexpr BCRSmatrix<T,BR,BC>::BCRSmatrix(const std::string& fname)
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
       
       bBR = BR*BR ;  
       bn  = denseRows*denseCols/(BR*BC) ;
      
       ia_.resize(denseRows/BR +1);
       ia_[0] = 1;

      
       for(std::size_t i = 0; i < dense.size() / BR ; i++)
       {    
          auto rowCount =0;
          for(std::size_t j = 0; j < dense[i].size() / BC ; j++)
          {
              if(validate_block(dense,i,j))
              {     
                    ja_.push_back(j+1);
                    insert_block(dense, i, j);
                    rowCount ++ ;
              }      
            
          }
          ia_[i+1] = ia_[i] + rowCount ;
       }
 
    }
   # ifdef __TESTING__
   printBCRS();
   # endif    
}


//-------------------------- methods  
//


//    operator one-start idx based
//
template <typename T, std::size_t BR, std::size_t BC>
T& BCRSmatrix<T,BR,BC>::operator()(const std::size_t row, const std::size_t col) noexcept
{
      assert(row > 0 && row <= denseRows && col > 0 && col <= denseCols );

      dummy = findValue(row-1,col-1) ; 
      return dummy ;
}

//  const version - operator one-start idx based
//
template <typename T, std::size_t BR, std::size_t BC>
const T& BCRSmatrix<T,BR,BC>::operator()(const std::size_t row, const std::size_t col) const noexcept
{
      assert(row > 0 && row <= denseRows && col > 0 && col <= denseCols );

      dummy = findValue(row-1,col-1) ; 
      return dummy ;
}



template <typename T,std::size_t BR,std::size_t BC>
inline auto constexpr BCRSmatrix<T,BR,BC>::printBlock(std::size_t i) const noexcept
{  
   auto w = i-1 ;
   auto k = 0;   
   for(std::size_t i = 0 ; i < BR ; ++i) 
   {
      for(std::size_t j=0 ; j < BC ; ++j )
      {
          std::cout << std::setw(8) << aa_.at(an_.at(w)-1+k) << ' ';
          k++;        
      }            
   }
}





// --- print out all the blocks 
//
template <typename T,std::size_t BR, std::size_t BC>
inline auto constexpr BCRSmatrix<T,BR,BC>::print_block(const std::vector<std::vector<T>>& dense,
                                                       std::size_t i, std::size_t j) const noexcept
{   
   for(std::size_t m = i * BR ; m < BR * (i + 1); ++m) 
   {
      for(std::size_t n = j * BC ; n < BC * (j + 1); ++n)
                  std::cout << dense[m][n] << ' ';
      std::cout << '\n';
   }
}


//   --- check for non zero block 
//
template <typename T,std::size_t BR, std::size_t BC>
inline auto constexpr BCRSmatrix<T,BR,BC>::validate_block(const std::vector<std::vector<T>>& dense,
                                                       std::size_t i, std::size_t j) const noexcept
{   
   bool nonzero = false ;
   for(std::size_t m = i * BR ; m < BR * (i + 1); ++m)
   {
      for(std::size_t n = j * BC ; n < BC * (j + 1); ++n)
      {
            if(dense[m][n] != 0) nonzero = true;
      }
   }
   return nonzero ;
}


//-----  insert block into Ba vector 
//
template <typename T,std::size_t BR, std::size_t BC>
inline auto constexpr BCRSmatrix<T,BR,BC>::insert_block(const std::vector<std::vector<T>>& dense,
                                                       std::size_t i, std::size_t j) noexcept
{   
   //std::size_t value = index;   
   bool firstElem = true ;
   for(std::size_t m = i * BR ; m < BR * (i + 1); ++m)
   {
      for(std::size_t n = j * BC ; n < BC * (j + 1); ++n)
      {    
            if(firstElem)
            {
                  an_.push_back(index+1);
                  firstElem = false ;
            }
            aa_.push_back(dense[m][n]);
            index ++ ;
      }
   }
}   
 
 
template <typename T, std::size_t BR,std::size_t BC> 
std::size_t constexpr BCRSmatrix<T,BR,BC>::findBlockIndex(const std::size_t r, const std::size_t c) const noexcept 
{
      for(auto j= ia_.at(r) ; j < ia_.at(r+1) ; j++ )
      {   
         if( ja_.at(j-1) == c+1  )
         {
            return j ;
         }
      }
}

// --- print the four vector of BCRS format
//
template <typename T, std::size_t BR, std::size_t BC>
auto constexpr BCRSmatrix<T,BR,BC>::printBCRS() const noexcept 
{ 

  std::cout << "aa_ :   " ;
  for(auto &x : aa_ ) 
      std::cout << x << ' ' ;
    std::cout << std::endl; 
  
  std::cout << "an_ :   " ;
  for(auto &x : an_ ) 
      std::cout <<  x << ' ' ;
    std::cout << std::endl; 
 
  std::cout << "ja_ :   " ;
  for(auto &x : ja_ ) 
      std::cout <<  x << ' ' ;
    std::cout << std::endl; 
   
   std::cout << "ia_ :   " ; 
   for(auto &x : ia_ ) 
      std::cout << x << ' ' ;
    std::cout << std::endl; 
      
}

template <typename T, std::size_t BR, std::size_t BC> 
void constexpr BCRSmatrix<T,BR,BC>::print() const noexcept 
{      
    for(auto i=0 ; i < denseRows ; i++)
    {
        for(auto j = 0; j < denseCols ; j++)
        {
              std::cout<< findValue(i, j) <<'\t';
        }
    std::cout << std::endl;
    }
}


template <typename T, std::size_t BR,std::size_t BC> 
auto constexpr BCRSmatrix<T,BR,BC>::recomposeMatrix() const noexcept
{

    std::vector<std::vector<T>> sparseMat(denseRows, std::vector<T>(denseCols, 0));
    auto BA_i = 0, AJ_i = 0;
    //for each BCSR row
    for(auto r = 0; r < denseRows/BR; r++){
        //for each Block in row
        for(auto nBlock = 0; nBlock < ia_.at(r+1)-ia_.at(r); nBlock++){  
            //for each subMatrix (Block)
            for(auto rBlock = 0; rBlock < BR; rBlock++){
                for(auto cBlock = 0; cBlock < BC; cBlock++){
                    //insert value
                    sparseMat.at(rBlock + r*BR).at(cBlock + (ja_.at(AJ_i)-1)*BC) = aa_.at(BA_i);
                    ++BA_i;
                }
            }
        ++AJ_i;
        }
    }
    return sparseMat;
}


template <typename T, std::size_t BR,std::size_t BC>
T constexpr BCRSmatrix<T,BR,BC>::findValue(const std::size_t i, const std::size_t j) const noexcept 
{
    auto index = findBlockIndex(i/BR, j/BC);
    if(index != 0)
        return aa_.at(an_.at(index-1)-1 + j%BC + (i%BR)*BC);
    else
        return T(0);
}





//
//----------------------   non member functions 

template <typename T, std::size_t BR,std::size_t BC>
std::ostream& operator<<(std::ostream& os , const BCRSmatrix<T,BR,BC>& m )
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



//
// perform (SpMV) Sparse-Matrix Vector product
//
template <typename T, std::size_t BR, std::size_t BC>
std::vector<T> operator*(const BCRSmatrix<T,BR,BC>& m, const std::vector<T>& x )
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
            auto brows = m.denseRows/BR ;  
            auto bnze  = m.an_.size()   ;

            auto z=0;

//# pragma omp parallel for   // using openmp library

            for(auto b=0 ; b < brows ; b++)
            {     
               for(auto j= m.ia_.at(b) ; j <= m.ia_.at(b+1)-1; j++ )
               {      
                  for(auto k=0 ; k < BR ; k++ )
                  {
                     for(auto t=0 ; t < BC ; t++)
                     {
                         y.at(BC*b+k) += m.aa_.at(z) * x.at(BC*(m.ja_.at(j-1)-1)+t) ;          
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
