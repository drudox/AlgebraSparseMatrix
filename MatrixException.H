# ifndef __MATRIX_EXCEPTION_H__
# define __MATRIX_EXCEPTION_H__

# include <exception>
# include <iostream>
# include <string>

class MatrixException : public std::exception {

   public:   
    
    explicit MatrixException(std::string mess) : std::exception{} , message{mess}
           {} 

    virtual ~MatrixException() noexcept = default ;
      
    inline virtual const char* what() const noexcept { return this-> message.c_str() ;}  

   protected: 
     
     std::string message ;
      

};

//-------------------------------------------------------------------------------

class InvalidSizeException : public MatrixException { 

      public:
       
         InvalidSizeException(const std::string &mess) : MatrixException{mess}
                  {}

};

class InvalidCoordinateException : public MatrixException {
      
      public:
       InvalidCoordinateException(const std::string &mess) : MatrixException{mess}     
                  {}
};

class OpeningFileException : public MatrixException {

      public :
       OpeningFileException(const std::string &mess) : MatrixException{mess}
               {}
};

# endif
