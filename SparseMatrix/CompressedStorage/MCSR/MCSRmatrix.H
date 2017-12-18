# ifndef __MODIFIED_COMPRESSED_ROW_STORAGE__
# define __MODIFIED_COMPRESSED_ROW_STORAGE__

# define __TESTING__

# include "../ModifiedCompressedMatrix.H"


namespace mg {
              namespace numeric {
                                  namespace algebra {


template <typename Type>
class MCSRmatrix;             

template <typename T>
std::ostream& operator<<(std::ostream& os ,const MCSRmatrix<T>& m) noexcept ;

template <typename T>
MCSRmatrix<T> operator+(const MCSRmatrix<T>& m1, const MCSRmatrix<T>& m2 ) ;

template <typename T>
std::vector<T> operator*(const MCSRmatrix<T>& A, const std::vector<T>& x ) noexcept ;

/*-------------------------------------------------------------------------------------
 *    
 *    -- Modified Compressed Sparse Row Matrix --
 *    
 *    @ Marco Ghiani November 2017  Glasgow
 *    
 *
 -------------------------------------------------------------------------------------*/


template <typename Type>
class MCSRmatrix 
                  : public ModifiedCompressedMatrix<Type>
{

   public:

      template <typename T>
      friend std::ostream& operator<<(std::ostream& os ,const MCSRmatrix<T>& m) noexcept ;

      template <typename T>
      friend MCSRmatrix<T> operator+(const MCSRmatrix<T>& m1, const MCSRmatrix<T>& m2 ) ;

      template <typename T>
      friend std::vector<T> operator*(const MCSRmatrix<T>& A, const std::vector<T>& x ) noexcept; 


   public: 

      constexpr MCSRmatrix( std::initializer_list<std::vector<Type>> rows);

      constexpr MCSRmatrix(const std::string& fname);
      
      constexpr MCSRmatrix(const std::size_t& ) noexcept; 

      virtual  ~MCSRmatrix() = default ;

      const Type& operator()(const std::size_t r , const std::size_t c) const noexcept override final;

      Type& operator()(const std::size_t r , const std::size_t c) noexcept override final;

      auto constexpr printMCSR() const noexcept ;
      
      void constexpr print() const noexcept override final;
      
      using ModifiedCompressedMatrix<Type>::printModCompressed;


   private:

      using SparseMatrix<Type>::aa_  ;
      using SparseMatrix<Type>::ja_  ;
      
      using SparseMatrix<Type>::nnz ;
      
      using ModifiedCompressedMatrix<Type>::dim ;
      


      std::size_t constexpr findIndex(const std::size_t, const std::size_t ) const noexcept override final;
     
      Type constexpr findValue(const std::size_t i, const std::size_t j) const noexcept override final;


};

//-------------------------------           Implementation       -------------------------------------------

template <typename T>
inline constexpr MCSRmatrix<T>::MCSRmatrix( std::initializer_list<std::vector<T>> rows)
{
      this->dim  = rows.size();
      auto _rows = *(rows.begin());

      aa_.resize(dim+1);
      ja_.resize(dim+1);

      if(dim != _rows.size())
      {
          throw InvalidSizeException("Error in costructor! Mod CSR format require square matrix!");  
      }

      std::size_t w = 0 ;
      ja_.at(w) = dim+2 ;
      for(auto ii = rows.begin(), i=1; ii != rows.end() ; ++ii, i++)
      {
          for(auto ij = ii->begin(), j=1, elemCount = 0 ; ij != ii->end() ; ++ij, j++ )   
          {
              if(i==j)
                 aa_[i-1] = *ij ;
              else if( i != j && *ij != 0 )
              {   
                 ja_.push_back(j); 
                 aa_.push_back(*ij); 
                 elemCount++ ;
              }
              ja_[i] = ja_[i-1] + elemCount;           
          }
      }  
  #ifdef __TESTING__    
      //printMCSR();
      printModCompressed();
  #endif    
}


//
//
template <typename T>
inline constexpr MCSRmatrix<T>::MCSRmatrix(const std::string& fname) 
{
   std::ifstream f(fname , std::ios::in );
   
   if(!f)
   {
      std::string mess = "Error unable to open file " + fname + 
                 "\nError occurred in MCSRmatrix constructor\n EXCEPTION THROWN" ;
      throw OpeningFileException(mess.c_str());
   }

   if( fname.find(".mtx") != std::string::npos )
   {
      auto i=0 , i1= 0 , j1= 0;    
      T elem =0;
      std::string line;
      getline(f,line); // jump out the header
      
      line = " " ;
      
      while(getline(f,line)) 
      {
         std::istringstream ss(line) ;  
         if(i==0)   
         {
            ss >> i1 ;
            ss >> j1 ;
            ss >> nnz ;
            
              if( i1 != j1 )
              {
                std::string mess ="Error in MCSR Matrix constructor:\n MCRS Matrix must be square! EXCEPTION THROWN" ;    
                throw InvalidSizeException(mess.c_str());      
              }
              else
              {
                dim = i1;
                aa_.resize(i1+1);
                ja_.resize(i1+1);
                
                for(i=0; i < dim+1 ; i++ )
                {
                     ja_.at(i) = dim+2 ;
                }  
              }  
         } 
         else
         {
            ss >> i1 >> j1 >> elem ;

            if(i1 == j1)      // diagonal element
            {
                aa_.at(i1-1) = elem ;  
            }
            else if( i1 != j1 )
            {
               for(i=i1 ; i < dim+1 ; i++)   
               {
                    ja_.at(i)++ ;
               }
               
               ja_.insert(ja_.begin() + (ja_.at(i1-1)-1) , j1  );
               aa_.insert(aa_.begin() + (ja_.at(i1-1)-1) , elem);   

            }
            else
            {
                  // nothing
            }
         }
         i++;
      }
      printModCompressed();
   }
   else
   {
      std::string line;
      std::string temp;
      std::size_t i=0 , j=0 ;
      while(getline(f,line))
      {
        std::istringstream ss(line);
        if(i==0) while(ss >> temp) j++ ;
        i++ ;
      }
      if(i != j) 
      {
         std::string mess = "Error in MCSR Matrix constructor:\n MCRS Matrix must be square! EXCEPTION THROWN" ;    
         throw InvalidSizeException(mess.c_str());   
      }
      else
      {
         f.clear();
         f.seekg( 0, std::ios::beg );
         
         this->dim = i;
         aa_.resize(dim+1); 
         ja_.resize(dim+1);
     
         std::size_t w = 0 ;
         ja_.at(w) = dim+2 ;
         T elem = 0;
         for(auto i=1; i <= dim  ; i++)
         {
          getline(f,line);
          std::istringstream ss(line);
          for(auto j=1, elemCount = 0 ; j <= dim ; j++ )   
          {
              ss >> elem ;           
              if(i==j)
                 aa_[i-1] = elem ;
              else if( i != j && elem != 0 )
              {   
                 ja_.push_back(j); 
                 aa_.push_back(elem); 
                 elemCount++ ;
              }
              ja_[i] = ja_[i-1] + elemCount;           
          }
         }     
# ifdef __TESTING__
   // printMCSR(); 
   printModCompressed();
# endif
   }
 }
}






template <typename T>
inline constexpr MCSRmatrix<T>::MCSRmatrix(const std::size_t& n ) noexcept
{
         this->dim = n;
         aa_.resize(dim+1); 
         ja_.resize(dim+1);

      for(std::size_t i = 0; i < aa_.size()-1 ; i++)
        aa_.at(i) = 1 ;    

      for(std::size_t i = 0; i < ja_.size() ; i++)
        ja_.at(i) = aa_.size()+1 ;   


}     


template <typename T>
inline std::size_t constexpr MCSRmatrix<T>::findIndex(const std::size_t row ,  const std::size_t col) const noexcept 
{    
     assert( row > 0 && row <= dim && col > 0 && col <= dim ); 
     if(row == col)
     {
       return row-1;
     }
     int i = -1;

      for(int i = ja_.at(row-1)-1 ; i < ja_.at(row)-1 ; i++ )
      {
            if( ja_.at(i) == col ) 
            { 
              return i ;
            }
      }
      return -1;
}


template <typename T>
inline auto constexpr MCSRmatrix<T>::printMCSR() const noexcept 
{
      std::cout << "AA :   " ;
      for(auto& x : aa_ )
          std::cout << x << ' ' ;
      std::cout << std::endl;

      std::cout << "JA :   " ;
      for(auto& x : ja_ )
          std::cout << x << ' ' ;
      std::cout << std::endl;
}


template <typename T>
inline void constexpr MCSRmatrix<T>::print() const noexcept
{
     for(std::size_t i=1 ; i <= dim ; i++ )
     {       
         for(std::size_t j=1 ; j <= dim ; j++)   
            std::cout << std::setw(8) << findValue(i,j) << ' ' ;
         std::cout << std::endl;   
     } 

}

// utility function 
template <typename T>      
inline T constexpr MCSRmatrix<T>::findValue(const std::size_t r, const std::size_t c) const noexcept
{
   auto i = findIndex(r,c); 
     if( i != -1 && i < aa_.size() )
     {      
            return aa_.at(i) ;
     }
     else
     {
            return aa_.at(dim) ;
     }      
}


template <typename T>
inline const T& MCSRmatrix<T>::operator()(const std::size_t r , const std::size_t c) const noexcept 
{     
     auto i = findIndex(r,c); 
     if( i != -1 && i < aa_.size() )
     {      
            return aa_.at(i) ;
     }
     else
     {
            return aa_.at(dim) ;
     }       
}


template <typename T>
inline T& MCSRmatrix<T>::operator()(const std::size_t r , const std::size_t c) noexcept 
{     
     auto i = findIndex(r,c); 
     if( i != -1 && i < aa_.size() )
     {      
            return aa_.at(i) ;
     }
     else
     {
            return aa_.at(dim) ;
     }       
}

//
//
//   non member operator

template <typename T>
std::ostream& operator<<(std::ostream& os ,const MCSRmatrix<T>& m) noexcept 
{
     for(std::size_t i=1 ; i <= m.dim ; i++ )
     {       
         for(std::size_t j=1 ; j <= m.dim ; j++)   
            os << std::setw(8) << m(i,j) << ' ' ;
         os << std::endl;   
     } 
     return os;
}


// perform sum by 2 Mod Comp Sparse Row matrices 
template <typename T>
MCSRmatrix<T> operator+(const MCSRmatrix<T>& m1, const MCSRmatrix<T>& m2 )
{
      if(m1.dim != m2.dim)
      {
          throw std::runtime_error("Matrix Dimension doesn't match");  
      }
      else
      {
         
         MCSRmatrix<T> res(m1.dim); 
      
         for(auto i=0 ; i < res.dim ; i++)
            res.aa_.at(i) = m1.aa_.at(i) + m2.aa_.at(i) ;
         
      


         std::set<unsigned int> elemCount ;
         for(auto i=0 ; i < res.dim ; i++)
         {
           elemCount.clear(); 
           auto n1 = m1.ja_.at(i+1) - m1.ja_.at(i);  
           auto n2 = m2.ja_.at(i+1) - m2.ja_.at(i);      
           
           //std::cout << "row " << i << " : " << n1 << "  " <<  n2 << std::endl;
          if(n1 && n2)
          { 
             auto j1 = m1.ja_.at(i) ;
             auto j2 = m2.ja_.at(i) ;
             
             auto i1 = m1.ja_.at(i) ;
             auto i2 = m2.ja_.at(i) ;
             
             
             std::set<unsigned int> m1i , m2i ;
              m1i.clear();
              m2i.clear();
              for(auto k=0; k < std::max(n1,n2); k++)
              {  
                 m1i.insert(m1.ja_.at(i1-1)); 
                 m2i.insert(m2.ja_.at(i2-1)); 
                 
                 if(i1 < m1.ja_.at(i+1)-1) i1++ ;
                 if(i2 < m2.ja_.at(i+1)-1) i2++ ;
              }   
              
             // 
             //
              std::vector<int> v_inters ;
              std::set_intersection(m1i.begin(),m1i.end(),
                                    m2i.begin(),m2i.end(),
                                    std::back_inserter(v_inters) );
              
             
              std::vector<int> v_diff1 ;
              std::set_difference(m1i.begin(),m1i.end(),m2i.begin(),m2i.end(),
                                  std::inserter(v_diff1,v_diff1.begin()) );
             
              std::vector<int> v_diff2 ;
              std::set_difference(m2i.begin(),m2i.end(),m1i.begin(),m1i.end(),
                                  std::inserter(v_diff2,v_diff2.begin()) );
                        
              
           
            for(auto j = 0 ; j < std::max(n1,n2) ; j++)
            {  
                
           

              if( m1.ja_.at(j1-1) == m2.ja_.at(j2-1) )
              {
                         
                  if(std::find(v_inters.begin(), v_inters.end() , m1.ja_.at(j1-1) ) != v_inters.end())  
                  {  
                        res.aa_.push_back( m1.aa_.at(j1-1) + m2.aa_.at(j2-1) );
                        res.ja_.push_back(m1.ja_.at(j1-1));
                        elemCount.insert(m1.ja_.at(j1-1)) ;
                  }
              }   
              else if (m1.ja_.at(j1-1) !=  m2.ja_.at(j2-1))  
              { 
                   if( std::find(v_diff1.begin(),v_diff1.end(),m1.ja_.at(j1-1)) != v_diff1.end() )
                   {
                           res.aa_.push_back( m1.aa_.at(j1-1));
                           res.ja_.push_back( m1.ja_.at(j1-1));
                           elemCount.insert(m1.ja_.at(j1-1)) ;

                   }
                   if( std::find(v_diff2.begin(),v_diff2.end(),m2.ja_.at(j2-1)) != v_diff2.end() ) // )
                   {
                           res.aa_.push_back( m2.aa_.at(j2-1));
                           res.ja_.push_back( m2.ja_.at(j2-1));
                           elemCount.insert(m2.ja_.at(j2-1)) ;
               
                   }
                   else
                   {
                        // nothing   
                   }
              }
              
                if(j1 < m1.ja_.at(i+1)-1) j1++ ; 
                if(j2 < m2.ja_.at(i+1)-1) j2++ ;
            }
          }
          else if(n2 && !n1 ) 
          {
            for(auto j=m2.ja_.at(i) ; j < m2.ja_.at(i+1) ; j++)
            {
               //   std::cout << m2.ja_.at(j-1) << ' '; 
                  elemCount.insert(m2.ja_.at(j-1)) ;
                  
                  res.aa_.push_back( m2.aa_.at(j-1));
                  res.ja_.push_back( m2.ja_.at(j-1));
                  
            }      
          }
          else if(n1 && !n2)
          {
            for(auto j=m1.ja_.at(i) ; j < m1.ja_.at(i+1) ; j++)
            {
              //    std::cout << m1.ja_.at(j-1) << ' '; 
                  elemCount.insert(m1.ja_.at(j-1)) ;
                  
                  res.aa_.push_back( m1.aa_.at(j-1));
                  res.ja_.push_back( m1.ja_.at(j-1));
            }
            //std::cout << std::endl; 
          }
          res.ja_.at(i+1) = res.ja_.at(i) + elemCount.size() ;  
         }

      return res;
      }
}

// perform   `Matrix \times Vector`
//
template <typename T>
std::vector<T> operator*(const MCSRmatrix<T>& A, const std::vector<T>& x ) noexcept 
{     
     assert(A.dim == x.size()); 
     std::vector<T> b(A.dim); 
     for(auto i=0; i < A.dim ; i++ )
         b.at(i) = A.aa_.at(i)* x.at(i) ;   
     
 
     for(auto i=0; i< A.dim ; i++)
     {
         for(auto k=A.ja_.at(i)-1 ; k< A.ja_.at(i+1)-1 ; k++ )
         {    
              b.at(i) += A.aa_.at(k)* x.at(A.ja_.at(k)-1);
             
         } 
     }
     return b;
}



  }//algebra
 }//numeric
}//mg 
# endif
