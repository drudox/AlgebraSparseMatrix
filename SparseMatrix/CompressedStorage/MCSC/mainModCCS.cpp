# include "MCSCmatrix.H"


using namespace std;
using namespace mg::numeric::algebra;

int main(){
  cout << "------------------------------------------------------------------------------" << endl;

  MCSCmatrix<int> m1 =  {{11,12,13,14,0,0},{0,22,23,0,0,0},{0,0,33,34,35,36},{0,0,0,44,45,0},
                           {0,0,0,0,0,56},{0,0,0,0,0,66}};
  MCSCmatrix<int> m111 =  {{11,12,13,14,0,0},{0,22,23,0,0,0},{0,0,33,34,35,36},{0,0,0,44,45,0},
                           {0,0,0,0,0,56},{0,0,0,0,0,66}};
  
 
  MCSCmatrix<int> m11 = {{0,44,9,3,0,0},{0,2,34,0,0,0},{0,0,11,24,0,8},{0,0,0,0,35,0},
                           {0,2,0,1,0,56},{0,0,4,3,0,11}};
    
      
      //{{}{}{}{}{}{}}
  cout << "------------------------------------------------------------------------------" << endl;
  MCSCmatrix<int> m2("input6.dat");
  cout << "------------------------------------------------------------------------------" << endl;
  MCSCmatrix<double> m3  = {{1.01, 0 , 0,0}, {0, 4.07, 1.14,0},{0,0,6.08,0},{1.06,0,2.2,9.9} };  //{0, 0, 0,  3.09}, {1.0, 2.4, 0,0}} ; 
  cout << "---------------------------------------------------------------------------------------" << endl; 
  MCSCmatrix<double> m4  = {{ 1.21, 1.06 , 0,2.15 }, {0, 3.31, 1.03,1.31},{0,1.06,-6.01,1.01},{4.12,1.1,1.06,-8.3}};// {
  cout << "---------------------------------------------------------------------------------------" << endl; 
   
  cout << m2;
  cout << "---------------------------------------------------------------------------------------" << endl; 
  cout << m3 ;
  cout << "---------------------------------------------------------------------------------------" << endl; 
  cout << m4 ;
  cout << "---------------------------------------------------------------------------------------" << endl; 
  //MCSCmatrix<int> m5 ("input5.dat");
  //cout << m5 ;
  cout << "---------------------------------------------------------------------------------------" << endl; 
  MCSCmatrix<double> m6 = m3+m4 ;   
  cout << m6 ;
  m6.printMCSC();
  cout << "---------------------------------------------------------------------------------------" << endl; 
  MCSCmatrix<int> m7 = m1+m11 ;  
  m7.print() ;
  m7.printMCSC();
  cout << "---------------------------------------------------------------------------------------" << endl; 
  
  MCSCmatrix<double> m100  = {{1.01, 0 , 2.34,0}, {0, 4.07, 0,0},{3.12,0,6.08,0},{1.06,0,2.2,9.9} };
  std::vector<double> v1={0,1.3,4.2,0.8};
  
  std::vector<double> v2 = m100*v1 ;
  cout << "---------------------------------------------------------------------------------------" << endl; 
  
  for(auto& x : v2)
    cout << x << ' ' ;
  cout << endl;
 
  MCSCmatrix<int> m12("input17.dat");
  std::vector<int> v3 =   {3,4,0,1,6,8,1,19};
  
 // std::cout << m12 ;
  m12.print(); 
  std::vector<int> v4 =  m12*v3 ;
  
  for (auto & x : v4 )
      cout << x << ' ' ;
  cout <<endl;   
 
  cout << "---------------------------------------------------------------------------------------" << endl; 
  MCSCmatrix<double> m13("coo_matrix.mtx");
  cout << m13 ;
  cout << "---------------------------------------------------------------------------------------" << endl; 
  cout << "---------------------------------------------------------------------------------------" << endl; 
  MCSCmatrix<double> m14("mat002.mtx");
  cout << m14 ;
  cout << "---------------------------------------------------------------------------------------" << endl; 
  MCSCmatrix<double> m15("mat003.mtx");
  cout << m15 ;
  cout << "---------------------------------------------------------------------------------------" << endl; 
  MCSCmatrix<double> m16("mat004.mtx");
  cout << m16 ;
  cout << "---------------------------------------------------------------------------------------" << endl; 
   
 
 return 0;
}


