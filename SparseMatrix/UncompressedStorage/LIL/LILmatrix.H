# ifndef __LIL_MATRIX_H__
# define __LIL_MATRIX_H__

# include "linklist.H"

# include "../../SparseMatrix.H"

namespace mg {
               namespace numeric {
                                    namespace algebra {

// forward declaration 
template <typename Type>
class LILmatrix ;

template <typename T>
std::ostream& operator<<(std::ostream& os , LILmatrix<T> m) noexcept ; 

template <typename T>
LILmatrix<T> operator+(const LILmatrix<T>& m1 ,const LILmatrix<T>& m2) ;

template <typename T>
std::vector<T> operator*(const LILmatrix<T>& A , const std::vector<T>& x ) ;

/*---------------------------------------------------------------------------
 *    LInked-List matrix class 
 *    
 *    @Marco Ghiani Glasgow Nov 2017
 *
 -------------------------------------------------------------------------*/


template <typename Type>
class LILmatrix 
                 : public SparseMatrix<Type>        
{


    template <typename T>
    friend std::ostream& operator<<(std::ostream& os , LILmatrix<T> m) noexcept ; 

    template <typename T>
    friend LILmatrix<T> operator+(const LILmatrix<T>& m1 ,const LILmatrix<T>& m2) ;

    template <typename T>
    friend std::vector<T> operator*(const LILmatrix<T>& A , const std::vector<T>& x ) ;



  public:
    
    constexpr LILmatrix(std::initializer_list<std::vector<Type>> row) noexcept ;
    
    constexpr LILmatrix(const std::string& fname ) ;
      
    constexpr LILmatrix(const std::size_t r, std::size_t c) noexcept ;

    virtual ~LILmatrix() = default ;  

    void constexpr print() const noexcept override final;  

    const Type& operator()(const std::size_t , const std::size_t )const noexcept override final;
    
    Type& operator()(const std::size_t , const std::size_t ) noexcept override final;


  private:
   
   using SparseMatrix<Type>::denseRows ;
   using SparseMatrix<Type>::denseCols ;
   
   std::vector< std::shared_ptr< linkedList<Type> >> aa_ ;
   
   using SparseMatrix<Type>::nnz  ;
   using SparseMatrix<Type>::zero ;
   using SparseMatrix<Type>::dummy ;   // used for store local variable to be returned as reference
                                       // orrible !                                                                
      
   Type constexpr findValue(const std::size_t r, std::size_t c) const noexcept override final{
            
            return aa_.at(r)->find_Value(c) ;
   }

};


//-----             implementation 

template <typename T>
inline constexpr LILmatrix<T>::LILmatrix(std::initializer_list<std::vector<T>> row) noexcept 
{
   this->denseRows = row.size();   
   this->denseCols = (*row.begin()).size();

   aa_.resize(denseRows);
   
   for(auto i=0 ; i< aa_.size() ; i++)
   {
      aa_.at(i) = std::make_shared<linkedList<T>>() ;      
   }

   std::size_t i=0 , j=0 ;   
   for(auto col : row)
   {
      j=0;
      for(auto elem : col)
      {
         if(elem != 0.0)
         {
            aa_.at(i)->insertNode(j,elem);
         }
       j++;     
      }
    i++;
   }
}      

template<typename T>
inline constexpr LILmatrix<T>::LILmatrix(const std::string& fname )  
{
      std::ifstream f(fname , std::ios::in);

      if(!f)
      {
            throw OpeningFileException("Error opening file in constructor ! EXCEPTION THROWED"); 
      }

      if(fname.find(".mtx") != std::string::npos )
      {
          std::string line ;   
          getline(f,line);    // jump out the header   
          
          line = " "; 
             
          auto i=0 , i1 =0 , j1 =0 ;
          T elem = 0;
          while(getline(f,line))
          {
             std::istringstream ss(line);
             if(i==0)
             {
                  ss >> denseRows ;
                  ss >> denseCols ;
                  ss >> nnz ;
                  
                  aa_.resize(denseRows);
                  for(auto i=0 ; i < aa_.size() ; i++)
                     aa_.at(i) = std::make_shared< linkedList<T>>();
                         
             }
             else
             {
                 ss >> i1 >> j1 >> elem ; 
                 
                 aa_.at(i1-1)->insertNode(j1-1,elem);

             }
             i++; 
          }
      
      }
      else
      {
         T elem = 0.;   
         auto i=0, j=0;   
         std::string line;   
         std::string temp;   
         while(getline(f,line))
         {  
            aa_.resize(i+1);
            aa_.at(i) = std::make_shared<linkedList<T>>() ;

            j=0;
            std::istringstream ss(line);
            if(i==0)             
            { 
                  while(ss >> elem)
                  {
                    if(elem != 0.0)
                    {
                        aa_.at(i)->insertNode(j,elem) ;
                    }                    
                    j++ ;
                  }
            this->denseCols = j;
           }
           else
           {
               while(ss >> elem)
               {
                  if(elem !=0);
                  {
                    aa_.at(i)->insertNode(j,elem) ;    
                  }
                  j++;
               }
           }
           i++;
         }
         this->denseRows = i; 
            
      }

}

//
//
template <typename T>
inline constexpr LILmatrix<T>::LILmatrix(const std::size_t r, std::size_t c) noexcept 
{
      this->denseRows = r ;
      this->denseCols = c ;
      
      aa_.resize(denseRows);

      for(auto i=0 ; i < denseRows ; i++)
         aa_.at(i) = std::make_shared<linkedList<T>>() ;
      
}





template <typename T>
inline void constexpr LILmatrix<T>::print() const noexcept
{
   for(auto i=0; i < aa_.size() ; i++)
   { 
    for(auto j=0; j< denseCols ; j++) 
    { 
      std::cout << std::setw(8) <<  this->operator()(i,j) << ' ' ; //aa_.at(i)->find_Value(j) << ' ';
    }
    std::cout <<  std::endl;
   }   
}


template <typename T>
inline const T& LILmatrix<T>::operator()(const std::size_t i, const std::size_t j)const noexcept 
{
    dummy = aa_.at(i)->find_Value(j) ;
    return dummy;  
}

template <typename T>
inline T& LILmatrix<T>::operator()(const std::size_t i, const std::size_t j) noexcept 
{
    dummy = aa_.at(i)->find_Value(j) ;
    return dummy ;  
}




// ==   non member function  
//
template <typename T>
std::ostream& operator<<(std::ostream& os , LILmatrix<T> m) noexcept  
{
   for(auto i=0; i < m.aa_.size() ; i++)
   { 
    for(auto j=0; j< m.denseCols ; j++) 
    { 
      os << std::setw(8) << m.aa_.at(i)->find_Value(j) << ' ';
    }
    os <<  std::endl;
   }   
   return os;   
}



template <typename T>
LILmatrix<T> operator+(const LILmatrix<T>& m1 ,const LILmatrix<T>& m2) 
{
      
    if( m1.size1() != m2.size1() || m1.size2() != m2.size2() )
    {
         throw InvalidSizeException("Error in operator + ! Matrix dimension doesn't match! ");     
    }
    else
    {
       LILmatrix<T> res(m1.size1(), m2.size2());   // aa_ already initialized !    

        for(auto i=0 ; i < res.size1() ; i++ )
        {
           for (auto j=0 ; j < res.size2() ; j++ )
           {
              auto v1 = m1.aa_.at(i)->find_Value(j);
              auto v2 = m2.aa_.at(i)->find_Value(j);
              

              if( v1 && v2 ) // both ave the same index
              {
                  res.aa_.at(i)->insertNode(j, v1+v2);    
              }
              else if(v1 && !v2)
              {
                 res.aa_.at(i)->insertNode(j,v1);
              }
              else if(v2 && !v1)
              {
                 res.aa_.at(i)->insertNode(j,v2); 
              }
              else if(!v1 && !v2)
              {
                  // nothing
              }
              else
              {
                  std::cerr << "Error in operator + the end is reached without handling a column index"
                            << std::endl;
              }
           }
        }    

    return res ;  
    }
}


// perform matrix times vector product 

template <typename T>
std::vector<T> operator*(const LILmatrix<T>& A , const std::vector<T>& x )
{
      if(A.size1() != x.size())
      {     
            std::string to = "x" ;
            std::string mess = "Error occured in operator* attempt to perfor productor between op1: "
                                 + std::to_string(A.size1()) + to + std::to_string(A.size2()) +
                                 " and op2: " + std::to_string(x.size());
            throw InvalidSizeException(mess.c_str());
      }
      else
      {
         std::vector<T> y(x.size(),0);        
         
         for(auto i=0; i< A.size1() ; i++ )
         {
            for(auto j= A.aa_.at(i)->begin() ; j != A.aa_.at(i)->end() ; ++j )
              y.at(i) += (*j) * x.at(j.index()); 
           //std::cout <<  "  value :  " << *j << "    index :   " << j.index() << '\n' ;  
         }
         return y;   
      }

}


  }//algebra 
 }//numeric
}// mg 
# endif
