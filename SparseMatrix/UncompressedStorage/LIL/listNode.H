# ifndef __LIST_NODE_H__
# define __LIST_NODE_H__

# include <memory>


namespace mg {
               namespace numeric { 
                                   namespace algebra {

template <typename T> class linkedList ;

template <typename T> 
class ListNode {

   public:
     
     friend class linkedList<T> ;
   
     constexpr ListNode(const std::size_t& i, const T& v ) noexcept : index_{i} , data_{v}, _next{nullptr} 
                  {}
     
     virtual ~ListNode() = default ;

     auto constexpr get() const noexcept {return data_;}
   
     
   private:

                      std::size_t index_ ; 
                                T  data_ ;
      
      std::shared_ptr<ListNode<T>> _next ; 

};

  } //algebra
 }// numeric
}//mg 
# endif






