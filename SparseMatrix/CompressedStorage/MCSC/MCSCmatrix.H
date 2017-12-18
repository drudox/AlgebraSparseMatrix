# ifndef __MODIFIED_COMPRESSED_COLUMNS_STORAGE__
# define __MODIFIED_COMPRESSED_COLUMNS_STORAGE__

# define __TESTING__ 

# include "../ModifiedCompressedMatrix.H"

namespace mg {
              namespace numeric {
                                 namespace algebra {

// forward declaration 
template <typename Type> class MCSCmatrix ;      


template <typename T>
std::ostream& operator<<(std::ostream& os ,const MCSCmatrix<T>& m ) noexcept ;

template <typename T>
MCSCmatrix<T> operator+(const MCSCmatrix<T>& m1 , const MCSCmatrix<T>& m2 );

template <typename T>
std::vector<T> operator*(const MCSCmatrix<T>& A ,const std::vector<T>& x)noexcept ;


/*----------------------------------------------------------------------------------
 *
 *      -- Modified Compressed sparse Columns Matrix -- 
 *    
 *       @ Marco Ghiani November 2017 Glasgow
 *
 *
 -----------------------------------------------------------------------------------*/


template <typename Type>
class MCSCmatrix 
                  : public ModifiedCompressedMatrix<Type> 
{

    public:

      template <typename T>
      friend std::ostream& operator<<(std::ostream& os ,const MCSCmatrix<T>& m ) noexcept ;

      template <typename T>
      friend MCSCmatrix<T> operator+(const MCSCmatrix<T>& m1 , const MCSCmatrix<T>& m2 );

      template <typename T>
      friend std::vector<T> operator*(const MCSCmatrix<T>& A ,const std::vector<T>& x) noexcept ;


    public:  
      
      constexpr MCSCmatrix(std::initializer_list<std::vector<Type>> row) ;
      
      constexpr MCSCmatrix(const std::string& ) ;
      
      constexpr MCSCmatrix(const std::size_t& n) noexcept ;
      
      virtual ~MCSCmatrix() = default;

      const Type& operator()(const std::size_t , const std::size_t ) const noexcept override final; 
      
      Type& operator()(std::size_t , std::size_t ) noexcept override final;
    
      auto constexpr printMCSC() const noexcept ;
      
      void constexpr print() const noexcept override final ;
      
      using ModifiedCompressedMatrix<Type>::printModCompressed;

    private:  
      
      using SparseMatrix<Type>::aa_  ;
      using SparseMatrix<Type>::ja_  ;
      
      using SparseMatrix<Type>::nnz ;
      
      using ModifiedCompressedMatrix<Type>::dim ;

     
      std::size_t constexpr findIndex(const std::size_t,const std::size_t) const noexcept override final;
      
      Type constexpr findValue(const std::size_t, const std::size_t ) const noexcept override final;

};


//----------------------------------   Implementation  --------------------------------------------//

template <typename T>
inline constexpr MCSCmatrix<T>::MCSCmatrix( std::initializer_list<std::vector<T>> row) 
{
      this->dim = row.size();
      auto il = *(row.begin());

      if(this-> dim != il.size())
      {
            throw InvalidSizeException("Matrix Must be square in Modified CSC format ");
      }
      
      std::vector<std::vector<T>> temp(row);
      
      aa_.resize(dim+1);
      ja_.resize(dim+1);

      ja_[0] = dim+2 ;
      auto elemCount = 0;
      for(auto c = 0 ; c < temp[0].size() ; c++ )
      {     
            elemCount =0 ;  
            for(auto r = 0  ; r < temp.size() ; r++) 
            {
                  if(c==r)
                  {
                     aa_[c] = temp[r][c] ;   
                  }
                  else if(c != r && temp[r][c] !=0)
                  {     
                      aa_.push_back(temp[r][c]);
                      ja_.push_back(r+1);
                      elemCount++ ;  
                  }
            }
            ja_[c+1] = ja_[c] + elemCount ;
      }
      
   # ifdef __TESTING__   
      printMCSC();
   #endif    
}


template <typename T>
inline constexpr MCSCmatrix<T>::MCSCmatrix(const std::string& fname ) 
{
    
    
   
    std::ifstream f(fname , std::ios::in);

    if(!f)
    {   
        std::string mess = "Error unable to open file " + fname + 
                 "\nError occurred in MCSRmatrix constructor\n EXCEPTION THROWN" ;

        throw OpeningFileException(mess.c_str());
    }
    
    if( fname.find(".mtx") != std::string::npos)
    {   
       auto i=0, i1=0 , j1 =0;
       T elem = 0;
       std::string line;
         
       getline(f,line); // jump out the header   
       
       line = " "; 
         
       while(getline(f,line))  
       {
         std::istringstream ss(line);   
         if(i==0)
         {
             ss >> i1 ;
             ss >> j1 ;
             ss >> nnz;
             
                 if(i1 != j1)
                 {
                   std::string mess ="Error in MCSR Matrix constructor:\n MCRS Matrix must be square! EXCEPTION THROWN" ;    
                   throw InvalidSizeException(mess.c_str());    
                 }
                 else
                 {
                      dim = j1 ;  
                      aa_.resize(j1+1);
                      ja_.resize(j1+1);
                 
                      for(auto j=0; j< dim+1 ; j++)
                      {
                         ja_.at(j) = dim+2;
                      }
                  }
         }
         else
         { 
              ss >> i1 >> j1 >> elem ;

              if(i1 == j1)    // diagonal element 
              {
                  aa_.at(j1-1) = elem ;
              }
              else if (i1 != j1)
              {
                 for(auto j=j1 ; j < dim+1 ; j++)
                 {
                     ja_.at(j)++ ;
                 } 
                  
                 ja_.insert(ja_.begin() + (ja_.at(j1-1)-1) , i1   ); 
                 aa_.insert(aa_.begin() + (ja_.at(j1-1)-1) , elem ); 
              }
              else
              {
                  // nothing
              }
          }
          i++ ;  
       }
    }
    else    // standard matrix input file
    {
       std::vector<std::vector<std::string>> temp ;
       std::string line, tmp ;

       while(getline(f,line))
       {
          std::istringstream ss(line);
          std::vector<std::string> vtmp ;
      
          vtmp.clear();
          while(ss >> tmp)
          {
              vtmp.push_back(tmp);   
          }
          temp.push_back(vtmp);
       }
       if(temp.size() != temp[0].size())
       {
           std::string mess = "Error in MCSR Matrix constructor:\n MCRS Matrix must be square! EXCEPTION THROWN" ;    
           throw InvalidSizeException(mess.c_str());   
       }
       else
       {
          this->dim = temp.size();

          aa_.resize(dim+1);
          ja_.resize(dim+1);  
     
          ja_[0] = dim+2 ;
          auto elemCount = 0;
    
          for(auto c = 0 ; c < temp[0].size() ; c++ )
          {     
              elemCount =0 ;  
               for(auto r = 0  ; r < temp.size() ; r++) 
               { 
                   if(c==r)
                   {
                      aa_[c] = std::stof(temp[r][c]) ;   
                   }
                   else if(c != r && std::stof(temp[r][c]) !=0)
                   {     
                      aa_.push_back(std::stof(temp[r][c]));
                      ja_.push_back(r+1);
                      elemCount++ ;  
                   }
               }
               ja_[c+1] = ja_[c] + elemCount ;
          }
        }
      }
   # ifdef __TESTING__
      printModCompressed();
   # endif   
}
//
//
template <typename T>
constexpr MCSCmatrix<T>::MCSCmatrix(const std::size_t& n) noexcept 
{
      dim = n ;      

      aa_.resize(dim+1);
      ja_.resize(dim+1);
      
      std::size_t i= 0; 
      for(auto i=0 ; i < dim ; i++)        // set to identity!
      {
         aa_.at(i) = 1. ;    
         ja_.at(i) = aa_.size()+1 ;   
      }
      aa_.at(dim) = 0;
      ja_.at(dim) = aa_.size()+1;
}


// utility function
//
template <typename T>
T constexpr MCSCmatrix<T>::findValue(const std::size_t r, const std::size_t c) const noexcept 
{
      auto i = findIndex(r,c);

      if(i != -1 && i < aa_.size())
      {
            return aa_.at(i);
      }
      else 
      {
           return aa_.at(dim); 
      }
}


//
//
template <typename T>
inline const T& MCSCmatrix<T>::operator()(const std::size_t r,const std::size_t c) const noexcept 
{
      auto i = findIndex(r,c);

      if(i != -1 && i < aa_.size())
      {
            return aa_.at(i);
      }
      else 
      {
           return aa_.at(dim); 
      }

}

// 
template <typename T>
inline T& MCSCmatrix<T>::operator()(std::size_t r, std::size_t c) noexcept 
{
      auto i = findIndex(r,c);

      if(i != -1 && i < aa_.size())
      {
            return aa_.at(i);
      }
      else 
      {
           return aa_.at(dim); 
      }

}

//
//
template <typename T>
inline std::size_t constexpr MCSCmatrix<T>::findIndex(const std::size_t row, const std::size_t col) const noexcept 
{
      assert(row > 0 && row <= dim && col > 0 && col <= dim);
      
      if(row == col)    // element of diag
      {     
          return col-1;         
      }
      int i = -1;

      for( auto i = ja_.at(col-1)-1 ; i < ja_.at(col)-1; i++ )
      {
            if( ja_.at(i) == row )
            {
               return i;   
            }
      }
      return -1 ;
}

// print function 
//
template <typename T>
inline void constexpr MCSCmatrix<T>::print() const noexcept 
{
      for(auto i=1 ; i <= dim ; i++ )
      {
            for(auto j=1 ; j <= dim ; j++)
                  std::cout << std::setw(8) << findValue(i,j) << ' ' ;
            std::cout << std::endl;
      }
}

// -- print vector aa -> vector of value and ja -> vector of index
template <typename T>
inline auto constexpr MCSCmatrix<T>::printMCSC() const noexcept 
{
      std::cout << "AA:  " ;     
      for(auto& x : aa_ ) 
            std::cout << x << ' ' ;
      std::cout << std::endl;
      
      std::cout << "JA:  " ;      
      for(auto& x : ja_ ) 
            std::cout << x << ' ' ;
      std::cout << std::endl;
 
}

template <typename T>
std::ostream& operator<<(std::ostream& os ,const MCSCmatrix<T>& m ) noexcept 
{
      for(auto i=1; i <= m.dim ; i++ )
      {     
            for(auto j=1; j <= m.dim ; j++)
            {
                  os << std::setw(8) << m(i,j) << ' ' ;
            }
            os << std::endl;
      }
      return os;
}

// perform sum by 2 Mod Comp Sparse Col matrices 
template <typename T>
MCSCmatrix<T> operator+(const MCSCmatrix<T>& m1, const MCSCmatrix<T>& m2 )
{
      if(m1.dim != m2.dim)
      {
          throw std::runtime_error("Matrix Dimension doesn't match");  
      }
      else
      {
         
         MCSCmatrix<T> res(m1.dim); 
      
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



template <typename T>
std::vector<T> operator*(const MCSCmatrix<T>& A ,const std::vector<T>& x) noexcept 
{
      assert(A.dim == x.size());
      std::vector<T> b(x.size());
      
      for(auto i=0 ; i < A.dim ; i++ )
            b.at(i) = A.aa_.at(i) * x.at(i) ; // diagonal value 

     for(auto i=0; i< A.dim ; i++)
     {
         for(auto k=A.ja_.at(i)-1 ; k < A.ja_.at(i+1)-1 ; k++ )
         {   
              b.at(A.ja_.at(k)-1) += A.aa_.at(k)* x.at(i);
         } 
     }
     return b;
}

  }//algebra
 }//numeric
}//mg
# endif 

