# ifndef __LINK_LIST_H__
# define __LINK_LIST_H__


# include "listNode.H"
# include <memory>
# include <iosfwd>
# include <iomanip>
# include <cassert>
# include <iostream>


namespace mg { 
               namespace numeric {
                                   namespace algebra {


template<typename T> class LILmatrix ;


template <typename T>
class linkedList {

   // public:
      template <typename U>      
      friend class LILmatrix ; 

    public:

     
      
      constexpr linkedList() noexcept  : nNodes{0}, _first{nullptr} 
                  {} 
      auto constexpr insertNode(const std::size_t& index, const T& value ) ;
      
      auto constexpr isEmpty() { return _first == nullptr ; }
    
      auto constexpr print() const noexcept ;
      
      auto constexpr printIndex() const noexcept ;
      
      auto constexpr size() const noexcept { return nNodes ; }
      
      struct iterator ;       // iterator 

      auto begin() noexcept { return iterator(*this, _first); }

      const auto cbegin() const noexcept { return iterator(*this, _first); }

      auto end() noexcept { return iterator(*this, nullptr); }

      const auto cend() const noexcept { return iterator(*this, nullptr); }

    private:
      
      std::size_t nNodes ;

      std::shared_ptr<ListNode<T>> _first ;
    //  std::share_ptr<ListNode<T>> _last  ;

      std::shared_ptr<ListNode<T>> getNewNode(const std::size_t& i , const T& v) {
      return std::make_shared<ListNode<T>>(i,v);
      } 
    
    
    public:  
      T constexpr find_Value(const std::size_t i)const noexcept ;
};


template<typename T>
auto constexpr linkedList<T>::insertNode(const std::size_t& index , const T& value) 
{
    auto newNode = getNewNode(index,value);  
      
    if(isEmpty())
    {
       _first = newNode ;
       newNode->_next = nullptr ;
       nNodes++;
    }
    else
    {
      auto current = _first ;
      
      if( _first->_next == nullptr ) //only one node 
      {
          if(_first->index_ > newNode->index_ )
          { 
             newNode->_next = _first  ; 
              _first        = newNode ; 
             nNodes++ ;
          }
          else if(_first->index_ < newNode->index_ )  
          {
              _first->_next  = newNode ;
              newNode->_next = nullptr ; 
              nNodes++;
          }
      }
      else
      {
         
         current = _first ;
         
         while(current->_next != nullptr)
         {
           if(current->_next->index_ < newNode->index_ )
               current = current->_next ;
           else 
           {
              break ;     
           }
         }   
         newNode->_next = current->_next ;
         current->_next = newNode       ;
         nNodes++ ;  
      }



    }
}

template<typename T>
auto constexpr linkedList<T>::print() const noexcept 
{
      auto current = _first;
      while(current != nullptr)
      {
            std::cout << current->data_  << "  " ;
            current = current->_next ;
      }
      std::cout << std::endl;
}

template<typename T>
auto constexpr linkedList<T>::printIndex() const noexcept 
{
      auto current = _first;
      while(current != nullptr)
      {
            std::cout << current->index_  << "  " ;
            current = current->_next ;
      }
      std::cout << std::endl;
}

template <typename T>
inline T constexpr linkedList<T>::find_Value(const std::size_t i)const noexcept 
{
   auto current = _first ;    
    
   while(current != nullptr )
   {
      if (current->index_ == i) return current->data_ ;   
      current = current->_next;
   }
   return 0.0;
}

//-----------  iterator class

template<typename T>
struct linkedList<T>::iterator {
      
      
      constexpr iterator( 
                           linkedList<T>* l , ListNode<T>* elem 
                         ) noexcept : _lptr{std::make_shared<linkedList<T>>(l)} , 
                                     _curr{std::make_shared<ListNode<T>>(elem)}
                   {} // 
    
      constexpr iterator(
                           linkedList<T>& l , std::shared_ptr<ListNode<T>> elem
                         ) noexcept : _lptr{std::make_shared<linkedList<T>>(l)} ,
                                   _curr(elem)
                   {} //
                
      constexpr iterator(
                          linkedList<T>& l
                        ) noexcept : _lptr{std::make_shared<linkedList<T>>(&l) } ,
                                    _curr(l._first)
                   {} //
       constexpr iterator(
                           linkedList<T>& l , iterator i
                          ) noexcept : _lptr{std::make_shared<linkedList<T>>(&l) },
                                       _curr{i}
                   {} //

       T& operator*() { return _curr->data_ ; }

       iterator operator++() { 
            _curr = _curr->_next ; 
            return *this;
       }
      
       constexpr auto index() const noexcept { return _curr->index_ ; }

       bool operator!=(const iterator& rhs) { return ( _curr != rhs._curr ); }

   private:
      
      std::shared_ptr<linkedList<T>> _lptr ;
      std::shared_ptr<ListNode<T>>   _curr ;  

};




  }//algebra
 }//numeric
}// mg 
# endif












