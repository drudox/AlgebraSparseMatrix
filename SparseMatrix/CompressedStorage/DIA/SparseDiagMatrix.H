# ifndef __COMPRESSED_DIAG_MATRIX_H__
# define __COMPRESSED_DIAG_MATRIX_H__

# include <map>
# include "../../SparseMatrix.H"

# define __TESTING__

namespace mg { 
                namespace numeric {
                                     namespace algebra {     


// forward declaration
template<typename T>
class DIAmatrix;


//
//
template <typename U>
std::ostream& operator<<(std::ostream& os, const DIAmatrix<U>& m );

//- SpMV
//

template<typename U>
std::vector<U> operator*(const DIAmatrix<U>& , const std::vector<U>& ) ;

/*---------------------------------------------------------------------------------------------
 *
 *    class Compressed Diagonal Storage varsion : 
 *       
 *     a) not fixed diag length (optimizing storage)
 *    
 *    @Marco Ghiani Dicembre 2017 Glasgow U.K
 *
 --------------------------------------------------------------------------------------------*/



template<typename Type >
class DIAmatrix 
                 : public SparseMatrix<Type>
{

      template <typename U>
      friend std::ostream& operator<<(std::ostream& os, const DIAmatrix<U>& m );

      template<typename U>
      friend std::vector<U> operator*(const DIAmatrix<U>& , const std::vector<U>& ) ;
 

 
   public:

      constexpr DIAmatrix(std::initializer_list<std::vector<Type>> ) ;
      
      constexpr DIAmatrix(const std::string& );
      
      virtual ~DIAmatrix() = default ;

      void constexpr print() const noexcept override final; 
   
      Type& operator()(const std::size_t , const std::size_t ) noexcept override ;
      
      const Type& operator()(const std::size_t , const std::size_t )const  noexcept override ;

      auto constexpr printDIA() const noexcept ; 
    
   private:
      
      using SparseMatrix<Type>::denseRows ;
      using SparseMatrix<Type>::denseCols ;
      
      using SparseMatrix<Type>::dummy     ;
      
      using SparseMatrix<Type>::nnz       ;

      std::size_t dim ;      

      std::map<int ,  std::vector<Type>> value ;   // value 
      std::set<int>                        dig ;   // off-diagonal distance

      Type constexpr findValue(const std::size_t , const std::size_t ) const noexcept override ;
};


template<typename T>
constexpr DIAmatrix<T>::DIAmatrix(std::initializer_list<std::vector<T>> rows )  
{
    denseRows = rows.size();
    denseCols = (*rows.begin()).size(); 
    
    if(denseRows != denseCols )
    {
       throw InvalidSizeException("DIA-Matrix , SIZE EXCEPTION THROWN :\n>>> Matrix Must be square! <<<");   
    }
    
    dim = denseRows ;
    
    value[0].resize(dim);
    value[1].resize(dim-1);     // diagonal 
    value[-1].resize(dim-1);      // tri-bands as default
    
    dig.insert(-1);
    dig.insert(0);            // tribands as default 
    dig.insert(1);
    

    auto i =0 , j=0 ;  
    for(auto row : rows )
    {
       j=0;     
       for(auto val : row )
       {
         if(val != 0)
         {
                if( i==j )   value[0].at(i) = val  ;
           else if(j==i-1)   value[-1].at(j) = val  ; 
           else if(j==i+1)   value[1].at(i) = val ;
           else if(j < i )
           {
              auto offset =  j-i ; 
              if( !(dig.find(offset) != dig.end()) )
              {
                  dig.insert(offset);
                  value[offset].resize(dim-(fabs(j-i)));
                  value[offset].at(j) = val ;
              }
              else
                  value[offset].at(j) = val ; 
           }
           else if(j > i)
           {
             auto offset = j-i ;
             if( !(dig.find(offset) != dig.end())) 
             {     
                 dig.insert(offset);
                 value[offset].resize(dim-(fabs(j-i))); 
                 value[offset].at(i) = val ; 
             }
             else
               value[offset].at(i) = val;
           }
         }
         j++;
       }
       i++;
    }
#  ifdef __TESTING__
      printDIA();
#  endif
}
      
template<typename T>      
constexpr DIAmatrix<T>::DIAmatrix(const std::string& fname)
{
      std::ifstream f(fname , std::ios::in);
      
      if(!f)
      {
         std::string mess = "Error opening file  " + fname +
                            "\n >>> Exception thrown in DIAmatrix constructor <<< " ;
         throw OpeningFileException(mess);
      }
      
      if( fname.find(".mtx") != std::string::npos)
      {

        std::string line;
        getline(f,line);      // jump through the header
        T val=0;

        auto i=0, j=0 , ii=0, jj=0;
        
        while(getline(f,line))
        {
            std::istringstream ss(line);
            
            if(ii == 0)
            {
                ss >> denseRows ;
                ss >> denseCols ;
                ss >> nnz;

                 if(denseRows != denseCols )
                 {
                      throw InvalidSizeException("DIA-Matrix , SIZE EXCEPTION THROWN :\n>>> Matrix Must be square! <<<");   
                 }
    
                 dim = denseRows ;
    
                 value[0].resize(dim);
                 value[1].resize(dim-1);     // diagonal 
                 value[-1].resize(dim-1);      // tri-bands as default
    
                 dig.insert(-1);
                 dig.insert(0);            // tribands as default 
                 dig.insert(1);

            }
            else         
            {
              ss >> i >> j >> val ;    
              
              i-- ;     
              j-- ;

                   if(i==j)   value[0].at(i) = val ;   
              else if(j==i-1) value[-1].at(j) = val  ; 
              else if(j==i+1) value[1].at(i) = val ;
              else if(j < i )
              {
                  auto offset =  j-i ; 
                  if( !(dig.find(offset) != dig.end()) )
                  {
                        dig.insert(offset);
                        value[offset].resize(dim-(fabs(j-i)));
                        value[offset].at(j) = val ;
                  }
                  else
                    value[offset].at(j) = val ; 
              }
              else if(j > i)
              {
                   auto offset = j-i ;
                   if( !(dig.find(offset) != dig.end())) 
                   {     
                       dig.insert(offset);
                       value[offset].resize(dim-(fabs(j-i))); 
                       value[offset].at(i) = val ; 
                   }
                   else
                      value[offset].at(i) = val;
              }

           jj++ ; 
           }
         ii++ ;
         }

      }
      else
      {
          std::string line, tmp;   
          T val =0 ;   
          auto i=0, j=0;
            
          while(getline(f,line))
          { 
            
            std::istringstream temp(line);
            if(i==0) while( temp >> tmp ) j++ ;
            i++;
          }

          denseCols = j;
          denseRows = i; 
          if(denseCols != denseRows)
          {
              throw InvalidSizeException("DIA-Matrix , SIZE EXCEPTION THROWN :\n>>> Matrix Must be square! <<<"); 
          }
          
          dim = denseCols ;
          f.clear();
          f.seekg(0, std::ios::beg);

          value[0].resize(dim);
          value[1].resize(dim-1);     // diagonal 
          value[-1].resize(dim-1);      // tri-bands as default
    
          dig.insert(-1);
          dig.insert(0);            // tribands as default 
          dig.insert(1);
             
          i=0 ; j=0; 
          while(getline(f,line))
          {
             std::istringstream ss(line);
             j=0;
             while(ss >> val ) 
             { 
                  if(val != 0)
                  {
                         if(i==j)   value[0].at(i) = val ;   
                    else if(j==i-1) value[-1].at(j) = val  ; 
                    else if(j==i+1) value[1].at(i) = val ;
                    else if(j < i )
                    {
                        auto offset =  j-i ; 
                        if( !(dig.find(offset) != dig.end()) )
                        {
                              dig.insert(offset);
                              value[offset].resize(dim) ; //-(fabs(j-i)));
                              value[offset].at(j) = val ;
                        }
                        else
                          value[offset].at(j) = val ; 
                    }
                    else if(j > i)
                    {
                         auto offset = j-i ;
                         if( !(dig.find(offset) != dig.end())) 
                         {     
                             dig.insert(offset);
                             value[offset].resize(dim-(fabs(j-i))); 
                             value[offset].at(i) = val ; 
                         }
                         else
                            value[offset].at(i) = val;
                    }

                  }
                  j++ ;
               }
               i++ ;
          }
            
      }
# ifdef __TESTING__
      printDIA();
# endif 
}

//-- private utility method
//
template<typename T>
T constexpr DIAmatrix<T>::findValue(const std::size_t r, const std::size_t c) const noexcept 
{
    T val = 0.0;  
    if( value.find(c-r) != value.end() )  
       val  = r < c ? value.at(c-r).at(r) : value.at(c-r).at(c) ;      
    return val ;     
}


//-- 
//
template<typename T>
void constexpr DIAmatrix<T>::print() const noexcept 
{
   for(auto i=0 ; i< denseRows ; i++)  
   {
      for(auto j=0 ; j< denseCols ; j++)
      {
         std::cout << std::setw(6) <<  findValue(i,j) << " " ;    
      }
      std::cout << std::endl;
   }      
}

//--
//
template<typename T>
auto constexpr DIAmatrix<T>::printDIA() const noexcept 
{
   
   auto it = dig.begin() ; 

   for( auto& n : dig)
    { 
      std::cout << std::setw(4) <<*it << " | :  " ;
      for(auto& x : value.at(n) )
      {
         std::cout << x << " " ;   
      }
      ++it ;
      std::cout << std::endl;
    }
}
 

//--
//
template <typename T>
T& DIAmatrix<T>::operator()(const std::size_t r, const std::size_t c) noexcept 
{
      assert( r >0 && r <= dim && c > 0 && c <= dim ) ;
      dummy = findValue(r-1,c-1);
      return dummy ;
}


//--
//
template<typename T> 
const T& DIAmatrix<T>::operator()(const std::size_t r, const std::size_t c) const noexcept 
{
      assert( r >0 && r <= dim && c > 0 && c <= dim ) ;
      dummy = findValue(r-1,c-1);
      return dummy ;
}



// -- nom member function
//

template <typename U>
std::ostream& operator<<(std::ostream& os, const DIAmatrix<U>& m )
{
      for(auto i=1; i <= m.denseRows ; i++ )
      {
          for(auto j=1 ; j <= m.denseCols ; j++) 
          {
               os << std::setw(6) << m(i,j) << " ";   
          }
          os << std::endl;  
      }
}



template<typename T>
std::vector<T> operator*(const DIAmatrix<T>& m, const std::vector<T>& x ) 
{
    std::vector<T> y(x.size());

    if(m.size2() != x.size())
    {
       std::string to = "x" ;
       std::string mess = "Error occured in operator* attempt to perfor productor between op1: "
                        + std::to_string(m.size1()) + to + std::to_string(m.size2()) +
                        " and op2: " + std::to_string(x.size());
       throw InvalidSizeException(mess.c_str());
    }
    else
    {
        for(auto j = m.dig.begin() ; j != m.dig.end() ; j++ )
        {   
            
             std::vector<T> vet(x.size());
             for(auto i=0; i< x.size() ; i++ )
             {
               vet.clear();

               if( *j < 0 )
               {
                  for(auto k=0 ; k < fabs(*j); k++ )
                     vet.push_back(0);
                  for(auto k=0 ; k< m.value.at(*j).size() ; k++ )
                     vet.push_back(m.value.at(*j).at(k));    
               }
               if(*j == 0)
               {
                 for(auto k=0 ; k< m.value.at(*j).size() ; k++ )
                     vet.push_back(m.value.at(*j).at(k));   
               }
               if(*j > 0)
               {
                 for(auto k=0 ; k< m.value.at(*j).size() ; k++ )
                     vet.push_back(m.value.at(*j).at(k));   
                 for(auto k=0 ; k < fabs(*j); k++ )
                     vet.push_back(0);
               }

               auto index = (i + *j) >= 0 ? (i + *j) : 0 ;   
                    y.at(i) += vet.at(i) * x[index] ;
               
            }
            //std::cout << std::endl;
        }
    }
    return y;
}




  }//algebra
 }//numeric
}// mg
# endif 
