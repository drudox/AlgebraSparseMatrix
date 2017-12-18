# ifndef __ELL_MATRIX_H__
# define __ELL_MATRIX_H__

# include "../../SparseMatrix.H"

# define  __DEBUG__
/*  # define __TESTING__ */


namespace mg {
              namespace numeric {
                                   namespace algebra {

// forward declaration
template <typename Type>
class ELLmatrix;

template <typename T>
std::ostream& operator<<(std::ostream& os , const ELLmatrix<T>& m );

template <typename T>
std::vector<T> operator*( const ELLmatrix<T>& , const std::vector<T>& );




//---------------------------------------
/*
 *          ELL matrix Class 
 *
 *    @Marco Ghiani Dec 2017 Glasgow
 *    
 ---------------------------------------*/


template <typename Type>
class ELLmatrix : 
                  public SparseMatrix<Type> 
{

      template <typename T>
      friend std::ostream& operator<<(std::ostream& os , const ELLmatrix<T>& m );

      template <typename T>
      friend std::vector<T> operator*( const ELLmatrix<T>& , const std::vector<T>& );




   public:
     
     constexpr ELLmatrix(std::initializer_list<std::initializer_list<Type>>, std::size_t = 3);
     
     constexpr ELLmatrix(const std::string& , std::size_t = 3);
     

     virtual Type& operator()(const std::size_t , const std::size_t) noexcept override ;

     virtual const Type& operator()(const std::size_t , const std::size_t) const noexcept override;
      
     void constexpr print() const noexcept override ;
 
     auto constexpr printELL() const noexcept ; 

   private:   
     
     using SparseMatrix<Type>::denseRows ;
     using SparseMatrix<Type>::denseCols ;
 
    
     std::vector<std::vector<Type>>        val_ ;
     std::vector<std::vector<std::size_t>> col_ ;
     
     using SparseMatrix<Type>::dummy ;
     using SparseMatrix<Type>::nnz ;

     Type constexpr findValue(const std::size_t , const std::size_t ) const noexcept override ; 

     std::size_t maxCols ;
};

//    ----------------------    Implementation 

template <typename T>
constexpr ELLmatrix<T>::ELLmatrix( std::initializer_list<std::initializer_list<T>> rows,
                                   std::size_t mx_col )   
                                                            : maxCols{mx_col}  
{
      
      this->denseRows = rows.size();
      this->denseCols = (*rows.begin()).size();
      
      val_.resize(denseRows);
      col_.resize(denseRows);

      auto i=0 ; auto j=1, k=0;  // column 1-start index cased
      for(auto row : rows)
      {
       j=1; k=0;   
       for(auto col : row )
       {    
          if(col != 0)   
          {
             val_.at(i).push_back(col);     
             col_.at(i).push_back(j);
             k++ ;
             if(k > maxCols)
             {
                throw InvalidSizeException(" >>> Matrix not conformal for ELL format! <<<");  
             }
             
          }
          ++j;  
       }
       if(val_.at(i).size() < maxCols) 
       {
            j = val_.at(i).size() ;
            while( val_.at(i).size() < maxCols )
            {
               val_.at(i).push_back(0.0) ;    
               col_.at(i).push_back(0.0) ;   
            }
       
       }
       ++i;
      }
# ifdef __DEBUG__
     printELL();  
# endif
}

//---

template<typename T>
constexpr ELLmatrix<T>::ELLmatrix(const std::string& fname , std::size_t mxcol ) 
                                                                                    : maxCols{mxcol} 
{
    std::ifstream f(fname , std::ios::in);  

    if(!f)
    {
         std::string mess = "Error opening file  " + fname +
                               "\n >>> Exception thrown in COOmatrix constructor <<< " ;
         throw OpeningFileException(mess);
    }
    
    if( fname.find(".mtx") != std::string::npos )
    {
       std::string line;     
       T elem =0;
       std::size_t i1=0 , j1=0 ;     
       getline(f,line); // jump out the header 
       line = " " ;  
       auto i=0;
       while(getline(f,line))
       {
          std::istringstream ss(line);  
          if(i==0)
          {
             ss >> denseRows >> denseCols >> nnz ; 
             val_.resize(denseRows);
             col_.resize(denseCols);
             i++;
             //maxCols = denseCols / 3 ;
          }
          else
          {
            ss >> i1 >> j1 >> elem ;
            
            col_.at(i1-1).push_back(j1);
            val_.at(i1-1).push_back(elem);
            
            if(val_.at(i1).size() > maxCols)
            {
                throw InvalidSizeException(" >>> Matrix not conformal for ELL format! <<<");  
            }
            ++i;
          }
       }

       for(i=0; i < denseRows ; i++ ){
         if( val_.at(i).size() < maxCols ){
           while(val_.at(i).size() < maxCols )
           {
               val_.at(i).push_back(0.0);   
               col_.at(i).push_back(0.0);   
           }
         }   
       }
    }
    else
    {

       std::string line;     
       T elem =0;
       auto i=0, j=1 ;
       while(getline(f,line))
       {  
          val_.push_back(std::vector<T>() );       
          col_.push_back(std::vector<std::size_t>() );       
          std::istringstream ss(line);
          j=1; auto k=0 ;
          while(ss >> elem)
          {
             if(elem != 0)
             {
                val_.at(i).push_back(elem);
                col_.at(i).push_back(j);
                k++;
                if( k > maxCols)
                {
                     throw InvalidSizeException(" >>> Matrix not conformal for ELL format! <<<");  
                }
             }
             ++j;    
          }
          denseCols = j-1;
          if(val_.at(i).size() < maxCols) 
          {
               j = val_.at(i).size() ;
               while( val_.at(i).size() < maxCols )
               {
                   val_.at(i).push_back(0.0) ;    
                   col_.at(i).push_back(0.0) ;   
               }
       
          }
          i++;
       }
       denseRows = i;
    }     
# ifdef __DEBUG__
      printELL();
# endif 
}





//--
template <typename T>
auto constexpr ELLmatrix<T>::printELL() const noexcept 
{

    std::cout <<  "   --  VAL  --    "  << std::endl;  
    for(auto& row : val_ )
    {
      for(auto& x : row )
      {
            std::cout << std::setw(4) << x << " " ;   
      }
      std::cout << std::endl;
    }
    
    std::cout << "  -- Col idx --     "  << std::endl;  
    for(auto& row : col_ )
    {
      for(auto& x : row )
      {
        std::cout << std::setw(4) << x << " " ;   
      }
      std::cout << std::endl;
    } 
}

//--  
template<typename T>  
T constexpr ELLmatrix<T>::findValue(const std::size_t i, const std::size_t j) const noexcept  
{
      assert(i >= 0 && i < denseRows &&
             j >= 0 && j < denseCols    );

      auto iter = std::find(col_.at(i).begin(),col_.at(i).end() , j+1 ); 
      std::size_t index ;   
      if(iter != col_.at(i).end())   
      {
         index = static_cast<std::size_t>( std::distance(col_.at(i).begin(), iter) ) ;    
         return val_.at(i).at(index) ;  
      }
      else
        return 0;
          
}     

template <typename T>
void constexpr ELLmatrix<T>::print() const noexcept 
{
      for(auto i=0 ; i < denseRows ; i++ ){
         for(auto j=0 ; j < denseCols ; j++){
            std::cout << std::setw(6) << this->operator()(i,j) << ' ' ;      
         }
         std::cout << std::endl;   
      }       

}


//--
template <typename T>
T& ELLmatrix<T>::operator()(const std::size_t r, const std::size_t c) noexcept 
{
    assert(r >= 0 && r < denseRows &&
           c >= 0 && c < denseCols );  
    dummy = findValue(r,c);  
      return dummy ;

}

//---
template <typename T>
const T& ELLmatrix<T>::operator()(const std::size_t r, const std::size_t c)const noexcept 
{
    assert(r >= 0 && r < denseRows &&
           c >= 0 && c < denseCols );  
  dummy = findValue(r,c);  
      return dummy ;

}


// non member function 
template <typename T>
std::ostream& operator<<(std::ostream& os , const ELLmatrix<T>& m )
{
   for(auto i=0 ; i < m.denseRows ; ++i){
      for(auto j=0 ; j < m.denseCols ; ++j){
          os << std::setw(6) << m(i,j) << "  " ;  
      }
      os << std::endl;
  }
  return os ;         
}

// perform matrix vector product
template <typename T>
std::vector<T> operator*( const ELLmatrix<T>& m, const std::vector<T>& x)
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
    
    //auto jj=0 , c=0;
    
    for(auto i=0 ; i < m.size1() ; i++){   
       for(auto j=0 ; j < m.maxCols ; j++){
            if(m.col_.at(i).at(j)!=0)
               y.at(i) += m.val_.at(i).at(j) * x.at(m.col_.at(i).at(j)-1); 
       }
    }

    return y;  
}









  }//algebra
 }//numeric
}//mg 
# endif
