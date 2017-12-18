# include "ELLmatrix.H"

using namespace std;
using namespace mg::numeric::algebra ;

int main(){
  
  ELLmatrix<int> ell1={{1,2,3,0,0,0},{0,4,5,0,6,0},{7,0,8,0,9,0},{0,8,0,0,7,6},{0,0,5,0,0,0},{0,0,4,0,3,0}};
  
  cout << "-------------------------------------------------------------------------- " << std::endl;
  cout << ell1 ; //.print();
  cout << "-------------------------------------------------------------------------- " << std::endl;
  ELLmatrix<int> ell2 = {{11,12,0,14,0,0,0,0},{0,22,23,0,25,0,0,0},{31,0,33,34,0,0,0,0},{0,42,0,0,45,46,0,0},
       {0,0,0,0,55,0,0,0},{0,0,0,0,65,66,67,0},{0,0,0,0,75,0,77,78}, {0,0,0,0,0,0,87,88}} ;
  cout << "-------------------------------------------------------------------------- " << std::endl;
  cout << ell2 ;

  ELLmatrix<double> ell3("input23.dat");
  cout << "-------------------------------------------------------------------------- " << std::endl;
  cout << ell3 ;
  
  cout << "-------------------------------------------------------------------------- " << std::endl;
  ELLmatrix<double> ell4("coo_matrix.mtx");
  cout << "-------------------------------------------------------------------------- " << std::endl;
  cout << ell4;
  cout << "-------------------------------------------------------------------------- " << std::endl;
  
  ELLmatrix<int> ell5("input17.dat",4); 
  cout << "-------------------------------------------------------------------------- " << std::endl;
  cout << ell5;
  cout << "-------------------------------------------------------------------------- " << std::endl;
  std::vector<int> v1 =  {3,4,0,1,6,8,1,19};    
  std::vector<int> v2 =  ell5*v1 ;
   for(auto& x : v2 )
      cout << x << ' ' ;
   std::cout << endl;   
   cout << "--------------------------------------------------------------------------------" << endl;
  

  
  
  return 0;
}



