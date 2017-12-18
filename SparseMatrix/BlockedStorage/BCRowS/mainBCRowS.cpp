# include "BCRowSmatrix.H"

using namespace std;
using namespace mg::numeric::algebra ;


int main() {

  BCRowSmatrix<int> vbcsr1 = {
                                {5,1,7,0,0},{0,1,0,2,3},{0,2,4,0,0},{0,0,1,3,0},{0,6,0,0,3}
                              };
  cout << "-----------------------------------------------------------------------------------" << endl;   

  BCRowSmatrix<int> vbcsr2("input19.dat");    
  cout << "-----------------------------------------------------------------------------------" << endl;   
  std::vector<int> x = {1,3,6,0,6};

  std::vector<int> b = vbcsr2 * x ;    

  for(auto& x : b )
      cout << x << ' ' ;
  cout << endl;    
  cout << "-----------------------------------------------------------------------------------" << endl;   
  cout << vbcsr1;
  
  BCRowSmatrix<int> vbcsr3("input17.dat");    
  cout << "-----------------------------------------------------------------------------------" << endl;   
  cout << vbcsr3 ;
  cout << "-----------------------------------------------------------------------------------" << endl;   
  
  std::vector<int> y = {3,4,0,1,6,8,1,19};

  std::vector<int> c = vbcsr3 * y ;    

  for(auto& x : c )
      cout << x << ' ' ;
  cout << endl;    
  cout << "-----------------------------------------------------------------------------------" << endl;   return 0;
}



